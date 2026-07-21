#include "6_io/PostDataWriter.h"

#include "0_basic/MPI_WRAPPER.h"
#include "0_basic/BoxOps.h"
#include "0_basic/LayoutTraits.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "1_grid/BlockTraits.h"
#include "2_topology/GlobalIncidence.h"
#include "2_topology/LocalIncidence.h"
#include "3_field/Field.h"
#include "6_io/RunData.h"
#include "7_metric/SingularEdgeRegistry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace POST
{
namespace
{
using TOPO::EntityAxis;
using TOPO::EntityDim;
using TOPO::EntityKey;
using I64 = std::int64_t;

struct EntityRow { I64 gid; EntityKey key; };

struct OperatorTriplet
{
    I64 row = -1;
    I64 column = -1;
    double value = 0.0;
    int destination = 0;
};

enum FaceStorageFlag : std::uint32_t
{
    FaceStorageGhost = 1u << 0,
    FaceStoragePhysical = 1u << 1,
    FaceStorageInterface = 1u << 2,
    FaceStorageCoupling = 1u << 3,
    FaceStorageParallel = 1u << 4
};

struct FaceStorageEntry
{
    EntityKey key;
    I64 column_id = -1;
    I64 quotient_id = -1;
    std::int8_t sign_to_owner = 1;
    std::uint8_t is_owner = 0;
    std::uint32_t flags = 0;
};

struct FaceStorageCatalog
{
    I64 quotient_count = 0;
    I64 global_column_count = 0;
    std::vector<FaceStorageEntry> entries;
    std::unordered_map<EntityKey, std::size_t, EntityKey::Hash> by_key;

    const FaceStorageEntry &at(const EntityKey &key) const
    {
        const auto it = by_key.find(key);
        if (it == by_key.end())
            throw std::runtime_error("PostData Bface storage catalog misses a solver stencil slot");
        return entries[it->second];
    }
};

std::string chunk_name(const char *stem, int rank)
{
    std::ostringstream s;
    s << stem << '_' << std::setw(4) << std::setfill('0') << rank << ".bin";
    return s.str();
}

std::string uuid_string(Uuid u)
{
    std::ostringstream s;
    s << std::hex << std::setfill('0') << std::setw(16) << u.hi << std::setw(16) << u.lo;
    return s.str();
}

Uuid mix_uuid(std::uint64_t seed, std::uint64_t a, std::uint64_t b)
{
    auto mix = [](std::uint64_t x) {
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27; x *= 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    };
    return {mix(seed ^ a), mix((seed << 1) ^ b ^ 0x9e3779b97f4a7c15ULL)};
}

std::size_t node_count(const TOPO::Topology &t) { return t.nodes.rep_to_qid.size(); }
std::size_t edge_count(const TOPO::Topology &t) { return t.edges.qkey_to_qid.size(); }
std::size_t face_count(const TOPO::Topology &t) { return t.faces.qkey_to_qid.size(); }
std::size_t cell_count(const TOPO::Topology &t) { return t.cells.local_to_qid.size(); }

WriteOptions normalize_options(const TOPO::Topology &t, WriteOptions o)
{
    if (o.mesh_uuid.hi == 0 && o.mesh_uuid.lo == 0)
        o.mesh_uuid = mix_uuid(0x4d50434e534d4553ULL,
                               (node_count(t) << 32) ^ edge_count(t),
                               (face_count(t) << 32) ^ cell_count(t));
    if (o.case_uuid.hi == 0 && o.case_uuid.lo == 0)
        o.case_uuid = mix_uuid(0x4d50434e53434153ULL, o.mesh_uuid.hi, o.mesh_uuid.lo);
    return o;
}

template <class T> std::span<const T> span(const std::vector<T> &v) { return {v.data(), v.size()}; }

std::vector<EntityRow> owner_rows(const TOPO::Topology &t, EntityDim dim, int rank)
{
    // Important: EdgeTopology::owner_to_gid/FaceTopology::owner_to_gid use a
    // separate compact ID space for shared-owner synchronization. They are not
    // quotient entity IDs and omit every ordinary unshared entity. Build the
    // output view from local representatives and resolve IDs through id_of().
    std::map<I64, EntityKey> unique;
    if (dim == EntityDim::Node) {
        for (const auto &[key, gid] : t.nodes.rep_to_qid)
            if (key.rank == rank) unique.emplace(gid, key);
    } else if (dim == EntityDim::Edge) {
        for (const auto &[key, qkey] : t.edges.local_to_qkey) {
            (void)qkey;
            if (key.rank == rank && t.is_owner(key))
                unique.emplace(static_cast<I64>(t.id_of(key).id), key);
        }
    } else if (dim == EntityDim::Face) {
        for (const auto &[key, qkey] : t.faces.local_to_qkey) {
            (void)qkey;
            if (key.rank == rank && t.is_owner(key))
                unique.emplace(static_cast<I64>(t.id_of(key).id), key);
        }
    } else {
        for (const auto &[key, gid] : t.cells.local_to_qid)
            if (key.rank == rank) unique.emplace(gid, key);
    }
    std::vector<EntityRow> rows; rows.reserve(unique.size());
    for (const auto &[id, key] : unique) rows.push_back({id, key});
    return rows;
}

void validate_entity_partition(const std::vector<EntityRow> &local_rows,
                               std::size_t global_count, const char *kind)
{
    if (global_count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error(std::string("PostData ") + kind + " count exceeds MPI int range");
    std::vector<int> local(global_count, 0), global(global_count, 0);
    for (const auto &row : local_rows) {
        if (row.gid < 0 || static_cast<std::size_t>(row.gid) >= global_count)
            throw std::runtime_error(std::string("PostData invalid local ") + kind + " output ID");
        ++local[static_cast<std::size_t>(row.gid)];
    }
    MPI_Allreduce(local.data(), global.data(), static_cast<int>(global_count),
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    for (std::size_t id = 0; id < global.size(); ++id)
        if (global[id] != 1)
            throw std::runtime_error(std::string("PostData ") + kind + " quotient ID " +
                                     std::to_string(id) + " has output multiplicity " +
                                     std::to_string(global[id]) + " (expected exactly 1)");
}

std::array<double,3> point(const Grid &g, const EntityKey &n)
{
    auto &b = const_cast<Grid &>(g).grids(n.block);
    return {b.x(n.i,n.j,n.k), b.y(n.i,n.j,n.k), b.z(n.i,n.j,n.k)};
}

std::array<double,3> center(const Grid &g, const EntityKey &e)
{
    if (e.dim == EntityDim::Node) return point(g,e);
    std::vector<EntityKey> nodes;
    if (e.dim == EntityDim::Edge) {
        for (const auto &x : TOPO::boundary_of_edge(e)) nodes.push_back(x.entity);
    } else if (e.dim == EntityDim::Face) {
        const auto c = TOPO::corners(e); nodes.assign(c.begin(), c.end());
    } else {
        for (int di=0;di<2;++di) for (int dj=0;dj<2;++dj) for (int dk=0;dk<2;++dk)
            nodes.push_back(TOPO::make_node(e.rank,e.block,e.i+di,e.j+dj,e.k+dk));
    }
    std::array<double,3> x{};
    for (const auto &n : nodes) { const auto p=point(g,n); for(int a=0;a<3;++a)x[a]+=p[a]; }
    for (double &v:x) v/=static_cast<double>(nodes.size());
    return x;
}

const char *metric_area_name(EntityAxis a)
{
    return a==EntityAxis::Xi ? "JDxi" : a==EntityAxis::Eta ? "JDet" : "JDze";
}
const char *metric_length_name(EntityAxis a)
{
    return a==EntityAxis::Xi ? "dl_xi" : a==EntityAxis::Eta ? "dl_eta" : "dl_zeta";
}
StaggerLocation location(EntityDim d, EntityAxis a)
{
    if(d==EntityDim::Node)return StaggerLocation::Node;
    if(d==EntityDim::Cell)return StaggerLocation::Cell;
    if(d==EntityDim::Edge)return a==EntityAxis::Xi?StaggerLocation::EdgeXi:a==EntityAxis::Eta?StaggerLocation::EdgeEt:StaggerLocation::EdgeZe;
    return a==EntityAxis::Xi?StaggerLocation::FaceXi:a==EntityAxis::Eta?StaggerLocation::FaceEt:StaggerLocation::FaceZe;
}
const char *location_name(StaggerLocation l)
{
    switch(l){case StaggerLocation::Cell:return "cell";case StaggerLocation::Node:return "node";
    case StaggerLocation::FaceXi:return "face_xi";case StaggerLocation::FaceEt:return "face_eta";case StaggerLocation::FaceZe:return "face_zeta";
    case StaggerLocation::EdgeXi:return "edge_xi";case StaggerLocation::EdgeEt:return "edge_eta";case StaggerLocation::EdgeZe:return "edge_zeta";} return "unknown";
}

int restart_location_code(StaggerLocation l)
{
    switch(l){case StaggerLocation::Cell:return 0;case StaggerLocation::FaceXi:return 1;
    case StaggerLocation::FaceEt:return 2;case StaggerLocation::FaceZe:return 3;
    case StaggerLocation::EdgeXi:return 4;case StaggerLocation::EdgeEt:return 5;
    case StaggerLocation::EdgeZe:return 6;case StaggerLocation::Node:return 7;}return -1;
}

I64 gid(const TOPO::Topology &t, const EntityKey &e) { return static_cast<I64>(t.id_of(e).id); }

const char *edge_alpha_name(EntityAxis a)
{
    return a==EntityAxis::Xi ? "Hodge_star_inverse_2form_to_1form_edge_xi_lumped" :
           a==EntityAxis::Eta ? "Hodge_star_inverse_2form_to_1form_edge_eta_lumped" :
                                "Hodge_star_inverse_2form_to_1form_edge_zeta_lumped";
}

const char *face_beta_name(EntityAxis a)
{
    return a==EntityAxis::Xi ? "Hodge_star_2form_to_1form_face_xi_lumped" :
           a==EntityAxis::Eta ? "Hodge_star_2form_to_1form_face_eta_lumped" :
                                "Hodge_star_2form_to_1form_face_zeta_lumped";
}

const char *face_field_name(EntityAxis a)
{
    return a==EntityAxis::Xi ? "B_xi" : a==EntityAxis::Eta ? "B_eta" : "B_zeta";
}

bool starts_with(std::string_view value, std::string_view prefix)
{
    return value.size() >= prefix.size() && value.substr(0, prefix.size()) == prefix;
}

bool contains_point(const Box3 &box, const EntityKey &key)
{
    return BOX::contains_point(box, Int3{key.i,key.j,key.k});
}

std::uint32_t face_storage_flags(const TOPO::Topology &t, const FieldBlock &field,
                                 const EntityKey &key)
{
    if (t.faces.local_to_qkey.find(key) != t.faces.local_to_qkey.end()) return 0;
    std::uint32_t flags=FaceStorageGhost;
    const auto loc=location(EntityDim::Face,key.axis);
    const int ngh=field.descriptor().nghost;
    for(const auto &p:t.inner_patches){
        if(p.this_rank!=key.rank||p.this_block!=key.block)continue;
        const Box3 ghost=LAYOUT::ghost_strip_from_node_box(loc,p.this_box_node,p.direction,ngh);
        if(!contains_point(ghost,key))continue;
        flags|=FaceStorageInterface;
        if(p.kind==TOPO::PatchKind::Parallel)flags|=FaceStorageParallel;
        if(p.is_coupling)flags|=FaceStorageCoupling;
    }
    for(const auto &p:t.parallel_patches){
        if(p.this_rank!=key.rank||p.this_block!=key.block)continue;
        const Box3 ghost=LAYOUT::ghost_strip_from_node_box(loc,p.this_box_node,p.direction,ngh);
        if(!contains_point(ghost,key))continue;
        flags|=FaceStorageInterface|FaceStorageParallel;
        if(p.is_coupling)flags|=FaceStorageCoupling;
    }
    for(const auto &p:t.physical_patches){
        if(p.this_rank!=key.rank||p.this_block!=key.block)continue;
        const Box3 inner=LAYOUT::boundary_inner_slab_one_layer_from_cells(
            GRID_TRAITS::cell_counts(field.get_block()),loc,p.this_box_node,p.direction);
        const Box3 ghost=LAYOUT::ghost_slab_from_inner(inner,p.direction,ngh);
        if(!contains_point(ghost,key))continue;
        if(starts_with(p.bc_name,"Coupled-"))flags|=FaceStorageCoupling|FaceStorageInterface;
        else flags|=FaceStoragePhysical;
    }
    return flags;
}

FaceStorageCatalog build_face_storage_catalog(const TOPO::Topology &t,
                                               const Field &fields,int rank)
{
    FaceStorageCatalog out;
    out.quotient_count=static_cast<I64>(face_count(t));
    I64 local_ghost_count=0,local_storage_count=0;
    for(int ib=0;ib<fields.num_blocks();++ib)for(int a=0;a<3;++a){
        const EntityAxis axis=static_cast<EntityAxis>(a);
        const auto &fb=fields.field(face_field_name(axis),ib);
        if(!fb.is_allocated())continue;
        const Int3 lo=fb.get_lo(),hi=fb.get_hi();
        for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k){
            ++local_storage_count;
            const EntityKey key=TOPO::make_face(rank,ib,i,j,k,axis);
            if(t.faces.local_to_qkey.find(key)==t.faces.local_to_qkey.end())++local_ghost_count;
        }
    }
    I64 ghost_prefix=0,global_ghost_count=0;
    MPI_Exscan(&local_ghost_count,&ghost_prefix,1,MPI_INT64_T,MPI_SUM,MPI_COMM_WORLD);
    if(rank==0)ghost_prefix=0;
    MPI_Allreduce(&local_ghost_count,&global_ghost_count,1,MPI_INT64_T,MPI_SUM,MPI_COMM_WORLD);
    out.global_column_count=out.quotient_count+global_ghost_count;
    out.entries.reserve(static_cast<std::size_t>(local_storage_count));
    I64 next_ghost=out.quotient_count+ghost_prefix;
    for(int ib=0;ib<fields.num_blocks();++ib)for(int a=0;a<3;++a){
        const EntityAxis axis=static_cast<EntityAxis>(a);
        const auto &fb=fields.field(face_field_name(axis),ib);
        if(!fb.is_allocated())continue;
        const Int3 lo=fb.get_lo(),hi=fb.get_hi();
        for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k){
            const EntityKey key=TOPO::make_face(rank,ib,i,j,k,axis);
            FaceStorageEntry entry;entry.key=key;
            if(t.faces.local_to_qkey.find(key)!=t.faces.local_to_qkey.end()){
                entry.quotient_id=gid(t,key);entry.column_id=entry.quotient_id;
                entry.sign_to_owner=static_cast<std::int8_t>(t.sign_to_owner(key));
                entry.is_owner=t.is_owner(key)?1:0;
            }else{
                entry.column_id=next_ghost++;
                entry.flags=face_storage_flags(t,fb,key);
            }
            out.by_key.emplace(key,out.entries.size());out.entries.push_back(entry);
        }
    }
    if(next_ghost!=out.quotient_count+ghost_prefix+local_ghost_count)
        throw std::runtime_error("PostData inconsistent Bface ghost storage count");
    return out;
}

const char *edge_dr_name(EntityAxis a)
{
    return a==EntityAxis::Xi ? "dr_xi" : a==EntityAxis::Eta ? "dr_eta" : "dr_zeta";
}

std::vector<TOPO::IncidenceEntry> adjoint_incident_faces(const EntityKey &e)
{
    // This is the exact local stencil used by CTOperators::CurlAdjFaceToEdge.
    if(e.axis==EntityAxis::Xi) return {
        {TOPO::make_face(e.rank,e.block,e.i,e.j,e.k-1,EntityAxis::Eta),+1},
        {TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Eta),-1},
        {TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Zeta),+1},
        {TOPO::make_face(e.rank,e.block,e.i,e.j-1,e.k,EntityAxis::Zeta),-1}};
    if(e.axis==EntityAxis::Eta) return {
        {TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Xi),+1},
        {TOPO::make_face(e.rank,e.block,e.i,e.j,e.k-1,EntityAxis::Xi),-1},
        {TOPO::make_face(e.rank,e.block,e.i-1,e.j,e.k,EntityAxis::Zeta),+1},
        {TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Zeta),-1}};
    return {
        {TOPO::make_face(e.rank,e.block,e.i,e.j-1,e.k,EntityAxis::Xi),+1},
        {TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Xi),-1},
        {TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Eta),+1},
        {TOPO::make_face(e.rank,e.block,e.i-1,e.j,e.k,EntityAxis::Eta),-1}};
}

std::vector<OperatorTriplet> route_triplets_to_row_owners(
    const std::vector<OperatorTriplet> &local)
{
    int nrank=1; MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    std::vector<int> send_count(nrank,0),recv_count(nrank,0);
    for(const auto &x:local) {
        if(x.destination<0||x.destination>=nrank)
            throw std::runtime_error("PostData operator row has invalid owner rank");
        ++send_count[x.destination];
    }
    MPI_Alltoall(send_count.data(),1,MPI_INT,recv_count.data(),1,MPI_INT,MPI_COMM_WORLD);
    std::vector<int> send_off(nrank,0),recv_off(nrank,0);
    for(int r=1;r<nrank;++r){send_off[r]=send_off[r-1]+send_count[r-1];recv_off[r]=recv_off[r-1]+recv_count[r-1];}
    const int nsend=send_off.back()+send_count.back(), nrecv=recv_off.back()+recv_count.back();
    std::vector<OperatorTriplet> ordered(nsend),received(nrecv);
    std::vector<int> cursor=send_off;
    for(const auto &x:local) ordered[cursor[x.destination]++]=x;
    std::vector<int> sb(nrank),rb(nrank),sd(nrank),rd(nrank);
    for(int r=0;r<nrank;++r){
        sb[r]=send_count[r]*static_cast<int>(sizeof(OperatorTriplet));
        rb[r]=recv_count[r]*static_cast<int>(sizeof(OperatorTriplet));
        sd[r]=send_off[r]*static_cast<int>(sizeof(OperatorTriplet));
        rd[r]=recv_off[r]*static_cast<int>(sizeof(OperatorTriplet));
    }
    MPI_Alltoallv(ordered.data(),sb.data(),sd.data(),MPI_BYTE,
                  received.data(),rb.data(),rd.data(),MPI_BYTE,MPI_COMM_WORLD);
    return received;
}

void validate_global_ids_and_orientation(const TOPO::Topology &t)
{
    auto validate_ids=[](const auto &map,std::size_t count,const char *kind){std::vector<bool>seen(count,false);for(const auto&[key,id]:map){(void)key;if(id<0||static_cast<std::size_t>(id)>=count)throw std::runtime_error(std::string("PostData invalid ")+kind+" global ID");seen[id]=true;}for(bool x:seen)if(!x)throw std::runtime_error(std::string("PostData non-contiguous ")+kind+" global IDs");};
    validate_ids(t.nodes.rep_to_qid,node_count(t),"node");
    validate_ids(t.edges.qkey_to_qid,edge_count(t),"edge");
    validate_ids(t.faces.qkey_to_qid,face_count(t),"face");
    validate_ids(t.cells.local_to_qid,cell_count(t),"cell");
    for(const auto&[key,s]:t.edges.local_to_qsign){(void)key;if(s!=1&&s!=-1)throw std::runtime_error("PostData invalid edge orientation sign");}
    for(const auto&[key,s]:t.faces.local_to_qsign){(void)key;if(s!=1&&s!=-1)throw std::runtime_error("PostData invalid face orientation sign");}
}

void write_geometry(const std::filesystem::path &path, const Grid &grid,
                    const TOPO::Topology &topo, const Field &fields,
                    const FaceStorageCatalog &storage,
                    const METRIC::SingularEdgeRegistry *singular, int rank,
                    const WriteOptions &o)
{
    PostBinaryWriter w(path, FileType::Geometry, o.case_uuid, o.mesh_uuid);
    const auto nodes=owner_rows(topo,EntityDim::Node,rank), edges=owner_rows(topo,EntityDim::Edge,rank);
    const auto faces=owner_rows(topo,EntityDim::Face,rank), cells=owner_rows(topo,EntityDim::Cell,rank);
    std::vector<I64> ids, ends; std::vector<double> xyz, centers, vectors, measures; std::vector<std::uint32_t> flags;
    for(const auto&r:nodes){ids.push_back(r.gid);auto p=point(grid,r.key);xyz.insert(xyz.end(),p.begin(),p.end());}
    w.WriteSection("node_global_id",1,span(ids)); w.WriteSection("node_xyz",3,span(xyz));
    ids.clear();
    for(const auto&r:edges){
        ids.push_back(r.gid); const auto ep=TOPO::endpoints(r.key);
        ends.push_back(gid(topo,ep.first));ends.push_back(gid(topo,ep.second));
        auto c=center(grid,r.key),p0=point(grid,ep.first),p1=point(grid,ep.second);
        centers.insert(centers.end(),c.begin(),c.end()); double l2=0;
        for(int a=0;a<3;++a){const double d=p1[a]-p0[a];vectors.push_back(d);l2+=d*d;}
        const double len=std::sqrt(l2); measures.push_back(len);
        std::uint32_t f=0; if(singular&&singular->find(static_cast<int>(r.gid))) f|=(1u<<1)|(1u<<0);
        flags.push_back(f);
        if(o.validate && (!std::isfinite(len) || (len<=0.0 && (f&(1u<<1))==0)))
            throw std::runtime_error("PostData edge geometry validation failed for gid " + std::to_string(r.gid));
    }
    w.WriteSection("edge_global_id",1,span(ids));w.WriteSection("edge_node_ids",2,span(ends));
    w.WriteSection("edge_center_xyz",3,span(centers));w.WriteSection("edge_directed_dr",3,span(vectors));
    w.WriteSection("edge_length",1,span(measures));w.WriteSection("edge_flags",1,span(flags));
    ids.clear();centers.clear();vectors.clear();measures.clear();flags.clear();
    for(const auto&r:faces){
        ids.push_back(r.gid);auto c=center(grid,r.key);centers.insert(centers.end(),c.begin(),c.end());
        const auto &fb=fields.field(metric_area_name(r.key.axis),r.key.block); double a2=0;
        const int s=topo.sign_to_owner(r.key);
        for(int q=0;q<3;++q){double v=s*fb(r.key.i,r.key.j,r.key.k,q);vectors.push_back(v);a2+=v*v;}
        const double area=std::sqrt(a2);measures.push_back(area);flags.push_back(0);
        if(o.validate&&(!std::isfinite(area)||area<=0))throw std::runtime_error("PostData face area validation failed");
    }
    w.WriteSection("face_global_id",1,span(ids));w.WriteSection("face_center_xyz",3,span(centers));
    w.WriteSection("face_area_vector",3,span(vectors));w.WriteSection("face_area",1,span(measures));w.WriteSection("face_flags",1,span(flags));
    ids.clear();centers.clear();vectors.clear();measures.clear();flags.clear();
    std::vector<I64> ghost_address;
    for(const auto &entry:storage.entries){
        if((entry.flags&FaceStorageGhost)==0)continue;
        ids.push_back(entry.column_id);
        ghost_address.insert(ghost_address.end(),{entry.key.rank,entry.key.block,
            restart_location_code(location(EntityDim::Face,entry.key.axis)),entry.key.i,entry.key.j,entry.key.k});
        const auto c=center(grid,entry.key);centers.insert(centers.end(),c.begin(),c.end());
        const auto &area=fields.field(metric_area_name(entry.key.axis),entry.key.block);double a2=0.0;
        for(int q=0;q<3;++q){const double v=area(entry.key.i,entry.key.j,entry.key.k,q);vectors.push_back(v);a2+=v*v;}
        measures.push_back(std::sqrt(a2));flags.push_back(entry.flags);
    }
    w.WriteSection("Bghost_global_id",1,span(ids));w.WriteSection("Bghost_address",6,span(ghost_address));
    w.WriteSection("Bghost_center_xyz",3,span(centers));w.WriteSection("Bghost_area_vector",3,span(vectors));
    w.WriteSection("Bghost_area",1,span(measures));w.WriteSection("Bghost_flags",1,span(flags));
    ids.clear();centers.clear();measures.clear();flags.clear();
    for(const auto&r:cells){ids.push_back(r.gid);auto c=center(grid,r.key);centers.insert(centers.end(),c.begin(),c.end());
        double v=fields.field("Jac",r.key.block)(r.key.i,r.key.j,r.key.k,0);measures.push_back(v);
        const std::string &physics=const_cast<Grid &>(grid).grids(r.key.block).block_name;
        flags.push_back(physics=="Fluid" ? 1u : physics=="Solid" ? 2u : 4u);
        if(o.validate&&(!std::isfinite(v)||v<=0))throw std::runtime_error("PostData cell volume validation failed");}
    w.WriteSection("cell_global_id",1,span(ids));w.WriteSection("cell_center_xyz",3,span(centers));
    w.WriteSection("cell_volume",1,span(measures));w.WriteSection("cell_flags",1,span(flags));w.Close();
}

template<class Map> void csr_from_map(const Map &rows, std::vector<I64>&offsets,std::vector<I64>&ids)
{
    offsets.assign(1,0);for(const auto &[key,row]:rows){(void)key;ids.insert(ids.end(),row.begin(),row.end());offsets.push_back(static_cast<I64>(ids.size()));}
}

void write_topology(const std::filesystem::path &path,const Grid&grid,const TOPO::Topology&t,
                    const FaceStorageCatalog&storage,int rank,const WriteOptions&o)
{
    PostBinaryWriter w(path,FileType::Topology,o.case_uuid,o.mesh_uuid);
    const auto edges=owner_rows(t,EntityDim::Edge,rank),faces=owner_rows(t,EntityDim::Face,rank),cells=owner_rows(t,EntityDim::Cell,rank);
    std::vector<I64> egid,enode,fgid,foff{0},feid;std::vector<std::int8_t> fsign;
    for(const auto&r:edges){egid.push_back(r.gid);auto ep=TOPO::endpoints(r.key);enode.push_back(gid(t,ep.first));enode.push_back(gid(t,ep.second));}
    for(const auto&r:faces){fgid.push_back(r.gid);for(const auto&e:TOPO::boundary_of_face(r.key)){feid.push_back(gid(t,e.entity));fsign.push_back(static_cast<std::int8_t>(t.sign_to_owner(r.key)*e.sign*t.sign_to_owner(e.entity)));}foff.push_back(feid.size());}
    std::vector<I64> cgid,coff{0},cfid;std::vector<std::int8_t> csign;
    std::map<I64,std::set<I64>> node_cells,edge_cells,face_cells;
    for(const auto&r:cells){cgid.push_back(r.gid);for(int di=0;di<2;++di)for(int dj=0;dj<2;++dj)for(int dk=0;dk<2;++dk)node_cells[gid(t,TOPO::make_node(rank,r.key.block,r.key.i+di,r.key.j+dj,r.key.k+dk))].insert(r.gid);
        for(const auto&f:TOPO::boundary_of_cell(r.key)){const I64 id=gid(t,f.entity);cfid.push_back(id);csign.push_back(static_cast<std::int8_t>(f.sign*t.sign_to_owner(f.entity)));face_cells[id].insert(r.gid);for(const auto&e:TOPO::boundary_of_face(f.entity))edge_cells[gid(t,e.entity)].insert(r.gid);}coff.push_back(cfid.size());}
    auto write_reverse=[&](std::string_view prefix,const auto&map){std::vector<I64> keys,off{0},vals;for(const auto&[k,s]:map){keys.push_back(k);vals.insert(vals.end(),s.begin(),s.end());off.push_back(vals.size());}w.WriteSection(std::string(prefix)+"_global_id",1,span(keys));w.WriteSection(std::string(prefix)+"_cell_offsets",1,span(off));w.WriteSection(std::string(prefix)+"_cell_ids",1,span(vals));};
    w.WriteSection("edge_global_id",1,span(egid));w.WriteSection("edge_node_ids",2,span(enode));
    w.WriteSection("face_global_id",1,span(fgid));w.WriteSection("face_edge_offsets",1,span(foff));w.WriteSection("face_edge_ids",1,span(feid));w.WriteSection("face_edge_signs",1,span(fsign));
    w.WriteSection("cell_global_id",1,span(cgid));w.WriteSection("cell_face_offsets",1,span(coff));w.WriteSection("cell_face_ids",1,span(cfid));w.WriteSection("cell_face_signs",1,span(csign));
    write_reverse("node",node_cells);write_reverse("edge",edge_cells);write_reverse("face",face_cells);
    // One section per block/location keeps each structured local array directly reshapeable.
    for(int ib=0;ib<grid.nblock;++ib){const auto&b=const_cast<Grid &>(grid).grids(ib);for(int ld=0;ld<8;++ld){EntityDim d;EntityAxis a=EntityAxis::None;Int3 sh{};
        if(ld==0){d=EntityDim::Node;sh=Int3{b.mx+1,b.my+1,b.mz+1};}else if(ld==1){d=EntityDim::Cell;sh=Int3{b.mx,b.my,b.mz};}
        else if(ld<5){d=EntityDim::Edge;a=static_cast<EntityAxis>(ld-2);sh=Int3{b.mx+(a!=EntityAxis::Xi),b.my+(a!=EntityAxis::Eta),b.mz+(a!=EntityAxis::Zeta)};}
        else{d=EntityDim::Face;a=static_cast<EntityAxis>(ld-5);sh=Int3{b.mx+(a==EntityAxis::Xi),b.my+(a==EntityAxis::Eta),b.mz+(a==EntityAxis::Zeta)};}
        std::vector<I64> map;std::vector<std::int8_t> signs,owners;map.reserve(sh.i*sh.j*sh.k);
        for(int k=0;k<sh.k;++k)for(int j=0;j<sh.j;++j)for(int i=0;i<sh.i;++i){EntityKey e=d==EntityDim::Node?TOPO::make_node(rank,ib,i,j,k):d==EntityDim::Cell?TOPO::make_cell(rank,ib,i,j,k):d==EntityDim::Edge?TOPO::make_edge(rank,ib,i,j,k,a):TOPO::make_face(rank,ib,i,j,k,a);map.push_back(gid(t,e));signs.push_back(t.sign_to_owner(e));owners.push_back(t.is_owner(e)?1:0);}
        std::ostringstream base;base<<'b'<<std::setw(4)<<std::setfill('0')<<ib<<'_'<<ld;std::vector<I64> shape{ib,ld,sh.i,sh.j,sh.k};
        w.WriteSection(base.str()+"_shape",5,span(shape));w.WriteSection(base.str()+"_gid",1,span(map));w.WriteSection(base.str()+"_sign",1,span(signs));w.WriteSection(base.str()+"_owner",1,span(owners));}}
    std::vector<I64> block_meta;for(int ib=0;ib<grid.nblock;++ib){const std::string&physics=const_cast<Grid &>(grid).grids(ib).block_name;const I64 code=physics=="Fluid"?1:physics=="Solid"?2:0;block_meta.insert(block_meta.end(),{rank,ib,code});}w.WriteSection("block_metadata",3,span(block_meta));
    std::vector<I64> conn;auto add_patch=[&](const TOPO::InterfacePatch&p){conn.insert(conn.end(),{p.this_block,p.nb_block,p.direction,p.nb_direction,p.trans.perm[0],p.trans.perm[1],p.trans.perm[2],p.trans.sign[0],p.trans.sign[1],p.trans.sign[2],p.trans.offset.i,p.trans.offset.j,p.trans.offset.k,p.this_rank,p.nb_rank});};
    for(const auto&p:t.inner_patches)if(p.this_rank==rank)add_patch(p);for(const auto&p:t.parallel_patches)if(p.this_rank==rank)add_patch(p);
    w.WriteSection("block_connections",15,span(conn));

    // Storage topology for the complete restart Bface arrays. Quotient slots
    // retain their canonical Face IDs; every ordinary ghost slot receives a
    // globally unique auxiliary column ID. Entries use the restart writer's
    // exact i,j,k loop order within each (field,block) array.
    std::vector<I64> storage_meta{storage.quotient_count,storage.global_column_count};
    std::vector<I64> storage_ids,storage_qids,storage_address;
    std::vector<std::int8_t> storage_sign;
    std::vector<std::uint8_t> storage_owner;
    std::vector<std::uint32_t> storage_flags;
    storage_ids.reserve(storage.entries.size());storage_qids.reserve(storage.entries.size());
    storage_sign.reserve(storage.entries.size());storage_owner.reserve(storage.entries.size());
    storage_flags.reserve(storage.entries.size());storage_address.reserve(storage.entries.size()*6);
    for(const auto &entry:storage.entries){
        storage_ids.push_back(entry.column_id);storage_qids.push_back(entry.quotient_id);
        storage_address.insert(storage_address.end(),{entry.key.rank,entry.key.block,
            restart_location_code(location(EntityDim::Face,entry.key.axis)),entry.key.i,entry.key.j,entry.key.k});
        storage_sign.push_back(entry.sign_to_owner);storage_owner.push_back(entry.is_owner);storage_flags.push_back(entry.flags);
    }
    w.WriteSection("Bstore_meta",2,span(storage_meta));w.WriteSection("Bstore_global_id",1,span(storage_ids));
    w.WriteSection("Bstore_quotient_id",1,span(storage_qids));w.WriteSection("Bstore_address",6,span(storage_address));
    w.WriteSection("Bstore_sign",1,span(storage_sign));w.WriteSection("Bstore_owner",1,span(storage_owner));
    w.WriteSection("Bstore_flags",1,span(storage_flags));

    std::vector<I64> physical_meta;std::vector<std::uint8_t> physical_names;
    for(const auto &p:t.physical_patches){
        if(p.this_rank!=rank)continue;const I64 name_offset=static_cast<I64>(physical_names.size());
        physical_names.insert(physical_names.end(),p.bc_name.begin(),p.bc_name.end());
        physical_meta.insert(physical_meta.end(),{p.this_rank,p.this_block,p.bc_id,p.direction,
            p.this_box_node.lo.i,p.this_box_node.lo.j,p.this_box_node.lo.k,
            p.this_box_node.hi.i,p.this_box_node.hi.j,p.this_box_node.hi.k,
            name_offset,static_cast<I64>(p.bc_name.size())});
    }
    w.WriteSection("physical_patch_meta",12,span(physical_meta));
    w.WriteSection("physical_patch_names",1,span(physical_names));

    std::vector<I64> interface_meta;
    auto add_interface_v2=[&](const TOPO::InterfacePatch&p){
        if(p.this_rank!=rank)return;
        interface_meta.insert(interface_meta.end(),{p.this_rank,p.nb_rank,p.this_block,p.nb_block,
            p.direction,p.nb_direction,static_cast<I64>(p.kind),p.is_coupling?1:0,
            p.trans.perm[0],p.trans.perm[1],p.trans.perm[2],p.trans.sign[0],p.trans.sign[1],p.trans.sign[2],
            p.trans.offset.i,p.trans.offset.j,p.trans.offset.k,
            p.this_box_node.lo.i,p.this_box_node.lo.j,p.this_box_node.lo.k,
            p.this_box_node.hi.i,p.this_box_node.hi.j,p.this_box_node.hi.k,
            p.nb_box_node.lo.i,p.nb_box_node.lo.j,p.nb_box_node.lo.k,
            p.nb_box_node.hi.i,p.nb_box_node.hi.j,p.nb_box_node.hi.k});
    };
    for(const auto&p:t.inner_patches)add_interface_v2(p);for(const auto&p:t.parallel_patches)add_interface_v2(p);
    w.WriteSection("interface_patch_meta",29,span(interface_meta));w.Close();
}

void write_reconstruction(const std::filesystem::path&path,const Grid&grid,const TOPO::Topology&t,const Field&fields,
                          const FaceStorageCatalog&storage,const METRIC::SingularEdgeRegistry *singular,int rank,const WriteOptions&o)
{
    PostBinaryWriter w(path,FileType::Reconstruction,o.case_uuid,o.mesh_uuid);
    const auto cells=owner_rows(t,EntityDim::Cell,rank);std::vector<I64> cg,off{0},faceids;std::vector<double> weights;
    for(const auto&r:cells){cg.push_back(r.gid);const auto&W=fields.field("Bcell_from_Bface_w",r.key.block);
        const std::array<TOPO::IncidenceEntry,6> stencil{{
            {TOPO::make_face(rank,r.key.block,r.key.i,r.key.j,r.key.k,EntityAxis::Xi),-1},
            {TOPO::make_face(rank,r.key.block,r.key.i+1,r.key.j,r.key.k,EntityAxis::Xi),+1},
            {TOPO::make_face(rank,r.key.block,r.key.i,r.key.j,r.key.k,EntityAxis::Eta),-1},
            {TOPO::make_face(rank,r.key.block,r.key.i,r.key.j+1,r.key.k,EntityAxis::Eta),+1},
            {TOPO::make_face(rank,r.key.block,r.key.i,r.key.j,r.key.k,EntityAxis::Zeta),-1},
            {TOPO::make_face(rank,r.key.block,r.key.i,r.key.j,r.key.k+1,EntityAxis::Zeta),+1}}};
        for(int n=0;n<6;++n){const auto&f=stencil[n];faceids.push_back(gid(t,f.entity));const int transform=f.sign*t.sign_to_owner(f.entity);for(int a=0;a<3;++a){const double value=W(r.key.i,r.key.j,r.key.k,6*a+n)*transform;if(o.validate&&!std::isfinite(value))throw std::runtime_error("PostData non-finite Bcell reconstruction weight");weights.push_back(value);}}off.push_back(faceids.size());}
    w.WriteSection("Bcell_global_id",1,span(cg));w.WriteSection("Bcell_offsets",1,span(off));w.WriteSection("Bcell_face_ids",1,span(faceids));w.WriteSection("Bcell_weights",3,span(weights));
    // Materialize the solver's scalar cell-to-node arithmetic average. Global
    // valences are reduced so rank chunks together form one normalized row.
    const std::size_t nn=node_count(t);std::vector<int> local_count(nn,0),global_count(nn,0);std::map<I64,std::set<I64>> rows;
    for(const auto&r:cells)for(int di=0;di<2;++di)for(int dj=0;dj<2;++dj)for(int dk=0;dk<2;++dk)rows[gid(t,TOPO::make_node(rank,r.key.block,r.key.i+di,r.key.j+dj,r.key.k+dk))].insert(r.gid);
    for(const auto&[n,s]:rows)local_count.at(n)=s.size();MPI_Allreduce(local_count.data(),global_count.data(),static_cast<int>(nn),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    std::vector<I64> ng,noff{0},nc;std::vector<double> nw;for(const auto&[n,s]:rows){ng.push_back(n);for(I64 c:s){nc.push_back(c);nw.push_back(1.0/global_count.at(n));}noff.push_back(nc.size());}
    if(o.validate){std::vector<double>local_sum(nn,0.0),global_sum(nn,0.0);for(const auto&[n,s]:rows)for(I64 c:s){(void)c;local_sum.at(n)+=1.0/global_count.at(n);}MPI_Allreduce(local_sum.data(),global_sum.data(),static_cast<int>(nn),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);for(std::size_t n=0;n<nn;++n)if(global_count[n]>0&&std::abs(global_sum[n]-1.0)>1e-12)throw std::runtime_error("PostData NodeScalar weights do not sum to one");}
    w.WriteSection("NodeScalar_global_id",1,span(ng));w.WriteSection("NodeScalar_offsets",1,span(noff));w.WriteSection("NodeScalar_cell_ids",1,span(nc));w.WriteSection("NodeScalar_weights",1,span(nw));

    // B_face -> J_edge. Build the coefficients from the final Hodge fields,
    // then reproduce the production alias reduction. Contributions are routed
    // to the rank owning the quotient edge so every serialized CSR row is
    // complete and occurs in exactly one rank chunk.
    std::unordered_map<EntityKey,std::size_t,EntityKey::Hash> edge_class_size;
    for(const auto &cls:t.edges.classes)
        for(const auto &m:cls.members) edge_class_size[m.entity]=cls.members.size();
    std::set<I64> singular_rows;
    if(singular) for(const auto &s:singular->entries()) singular_rows.insert(gid(t,s.owner));
    std::vector<OperatorTriplet> local_bj;
    for(const auto &[e,qkey]:t.edges.local_to_qkey){
        (void)qkey;if(e.rank!=rank||singular_rows.count(gid(t,e)))continue;
        const auto size_it=edge_class_size.find(e);
        const double divisor=size_it==edge_class_size.end()?1.0:static_cast<double>(size_it->second);
        const auto &alpha=fields.field(edge_alpha_name(e.axis),e.block);
        const double row_scale=t.sign_to_owner(e)*alpha(e.i,e.j,e.k,0)/divisor;
        const int destination=t.owner_of(e).rank;
        for(const auto &inc:adjoint_incident_faces(e)){
            const auto &input=storage.at(inc.entity);
            const auto &beta=fields.field(face_beta_name(inc.entity.axis),inc.entity.block);
            const double input_transform=input.quotient_id>=0?input.sign_to_owner:1.0;
            const double value=row_scale*inc.sign*input_transform*
                               beta(inc.entity.i,inc.entity.j,inc.entity.k,0);
            local_bj.push_back({gid(t,e),input.column_id,value,destination});
        }
    }
    if(singular) for(const auto &s:singular->entries()){
        int owner_orientation=0;
        for(const auto &alias:s.aliases)if(alias.owner){owner_orientation=alias.orientation;break;}
        if(owner_orientation!=1&&owner_orientation!=-1)
            throw std::runtime_error("PostData singular edge registry has no valid owner orientation");
        for(const auto &inc:s.local_incident_faces){
            int incidence=0;
            for(const auto &x:adjoint_incident_faces(inc.source_alias))
                if(x.entity==inc.entity){incidence=x.sign;break;}
            if(incidence==0)continue;
            const auto &input=storage.at(inc.entity);
            const auto &beta=fields.field(face_beta_name(inc.entity.axis),inc.entity.block);
            // SingularEdgeRegistry writes owner_orientation * global_sum to
            // the owner-local Jedge slot after reducing source-oriented face
            // contributions. Preserve both signs in the serialized row.
            const double value=s.inverse_hodge*owner_orientation*inc.source_orientation*incidence*
                               (input.quotient_id>=0?input.sign_to_owner:1.0)*
                               beta(inc.entity.i,inc.entity.j,inc.entity.k,0);
            local_bj.push_back({gid(t,s.owner),input.column_id,value,s.owner.rank});
        }}
    const auto received_bj=route_triplets_to_row_owners(local_bj);
    std::map<I64,std::map<I64,double>> bj_rows;
    for(const auto &x:received_bj)bj_rows[x.row][x.column]+=x.value;
    // Reproduce the final Pole Jedge physical handler. Tangential/rotational
    // edge values on the active Pole slab are exactly zero. The normal family
    // uses the handler's uncapped dl/Astar formula when its legacy beta_zeta
    // field is present; otherwise the production handler leaves the inner
    // value from Calc_J_Edge unchanged.
    for(const auto &p:t.physical_patches){
        if(p.this_rank!=rank||p.bc_name!="Pole")continue;const int normal=std::abs(p.direction)-1;if(normal<0||normal>1)continue;
        const auto &blk=const_cast<Grid&>(grid).grids(p.this_block);const Int3 nc{blk.mx,blk.my,blk.mz};
        for(int ea=0;ea<3;++ea){const EntityAxis axis=static_cast<EntityAxis>(ea);const Box3 box=LAYOUT::boundary_inner_slab_one_layer_from_cells(nc,location(EntityDim::Edge,axis),p.this_box_node,p.direction);
            for(int i=box.lo.i;i<box.hi.i;++i)for(int j=box.lo.j;j<box.hi.j;++j)for(int k=box.lo.k;k<box.hi.k;++k){const EntityKey e=TOPO::make_edge(rank,p.this_block,i,j,k,axis);if(!t.is_owner(e))continue;
                auto &row=bj_rows[gid(t,e)];if(ea!=normal){row.clear();continue;}
                if(!fields.has_field("beta_zeta"))continue;const auto&beta=fields.field("beta_zeta",p.this_block);const auto&dl=fields.field(metric_length_name(axis),p.this_block);const auto&Astar=fields.field(axis==EntityAxis::Xi?"Astar_xi":"Astar_eta",p.this_block);const double area=Astar(i,j,k,0);row.clear();if(std::abs(area)<=1.e-300)continue;const double scale=dl(i,j,k,0)/area;
                const EntityKey f0=normal==0?TOPO::make_face(rank,p.this_block,i,j,k,EntityAxis::Zeta):TOPO::make_face(rank,p.this_block,i-1,j,k,EntityAxis::Zeta);
                const EntityKey f1=normal==0?TOPO::make_face(rank,p.this_block,i,j-1,k,EntityAxis::Zeta):TOPO::make_face(rank,p.this_block,i,j,k,EntityAxis::Zeta);
                const auto &in0=storage.at(f0),&in1=storage.at(f1);
                row[in0.column_id]+=scale*beta(f0.i,f0.j,f0.k,0)*(in0.quotient_id>=0?in0.sign_to_owner:1.0);
                row[in1.column_id]-=scale*beta(f1.i,f1.j,f1.k,0)*(in1.quotient_id>=0?in1.sign_to_owner:1.0);
            }}
    }
    std::vector<I64> jeg,joff{0},jf;std::vector<double> jw;
    for(const auto &r:owner_rows(t,EntityDim::Edge,rank)){
        jeg.push_back(r.gid);const auto it=bj_rows.find(r.gid);
        if(it!=bj_rows.end())for(const auto &[col,value]:it->second){
            if(o.validate&&!std::isfinite(value))throw std::runtime_error("PostData non-finite B_face_to_J_edge weight");
            if(value!=0.0){jf.push_back(col);jw.push_back(value);}
        }joff.push_back(jf.size());
    }
    w.WriteSection("BfaceJedge_global_id",1,span(jeg));w.WriteSection("BfaceJedge_offsets",1,span(joff));
    w.WriteSection("BfaceJedge_face_ids",1,span(jf));w.WriteSection("BfaceJedge_weights",1,span(jw));

    // J_edge -> cell Cartesian vector. Start with the exact preassembled 36
    // weights and replace Pole rows with the same regularized ring least-
    // squares operator used by MercurySolver::calc_Jcell().
    using Vec3=std::array<double,3>;
    std::map<I64,std::map<I64,Vec3>> jc_rows;
    auto add_jc=[&](std::map<I64,Vec3>&row,const EntityKey&e,const Vec3&v){
        auto &dst=row[gid(t,e)];const double s=t.sign_to_owner(e);
        for(int a=0;a<3;++a)dst[a]+=s*v[a];
    };
    auto cell_edges=[&](const EntityKey&c){return std::array<EntityKey,12>{{
        TOPO::make_edge(rank,c.block,c.i,c.j,c.k,EntityAxis::Xi),TOPO::make_edge(rank,c.block,c.i,c.j+1,c.k,EntityAxis::Xi),
        TOPO::make_edge(rank,c.block,c.i,c.j,c.k+1,EntityAxis::Xi),TOPO::make_edge(rank,c.block,c.i,c.j+1,c.k+1,EntityAxis::Xi),
        TOPO::make_edge(rank,c.block,c.i,c.j,c.k,EntityAxis::Eta),TOPO::make_edge(rank,c.block,c.i+1,c.j,c.k,EntityAxis::Eta),
        TOPO::make_edge(rank,c.block,c.i,c.j,c.k+1,EntityAxis::Eta),TOPO::make_edge(rank,c.block,c.i+1,c.j,c.k+1,EntityAxis::Eta),
        TOPO::make_edge(rank,c.block,c.i,c.j,c.k,EntityAxis::Zeta),TOPO::make_edge(rank,c.block,c.i+1,c.j,c.k,EntityAxis::Zeta),
        TOPO::make_edge(rank,c.block,c.i,c.j+1,c.k,EntityAxis::Zeta),TOPO::make_edge(rank,c.block,c.i+1,c.j+1,c.k,EntityAxis::Zeta)}};};
    for(const auto &r:cells){auto &row=jc_rows[r.gid];const auto edges=cell_edges(r.key);const auto&W=fields.field("Jcell_from_Jedge_w",r.key.block);
        for(int n=0;n<12;++n)add_jc(row,edges[n],{W(r.key.i,r.key.j,r.key.k,n),W(r.key.i,r.key.j,r.key.k,12+n),W(r.key.i,r.key.j,r.key.k,24+n)});}

    auto solve3=[&](double Ain[3][3],const double bin[3],double x[3])->bool{
        double A[3][3];for(int a=0;a<3;++a)for(int c=0;c<3;++c)A[a][c]=Ain[a][c];
        const double reg=1.e-12*std::max(1.0,A[0][0]+A[1][1]+A[2][2]);for(int a=0;a<3;++a)A[a][a]+=reg;
        double M[3][4]={{A[0][0],A[0][1],A[0][2],bin[0]},{A[1][0],A[1][1],A[1][2],bin[1]},{A[2][0],A[2][1],A[2][2],bin[2]}};
        for(int c=0;c<3;++c){int piv=c;double am=std::abs(M[c][c]);for(int rr=c+1;rr<3;++rr)if(std::abs(M[rr][c])>am){am=std::abs(M[rr][c]);piv=rr;}if(am<1.e-300)return false;
            if(piv!=c)for(int q=c;q<4;++q)std::swap(M[c][q],M[piv][q]);const double inv=1.0/M[c][c];for(int q=c;q<4;++q)M[c][q]*=inv;
            for(int rr=0;rr<3;++rr)if(rr!=c){const double fac=M[rr][c];for(int q=c;q<4;++q)M[rr][q]-=fac*M[c][q];}}
        for(int a=0;a<3;++a)x[a]=M[a][3];return true;};
    auto pole_row=[&](int ib,const std::vector<EntityKey>&edges)->std::map<I64,Vec3>{
        double A[3][3]={{0,0,0},{0,0,0},{0,0,0}};struct Term{EntityKey e;double q[3];};std::vector<Term>terms;
        for(const auto&e:edges){const auto&dr=fields.field(edge_dr_name(e.axis),ib);const double dx=dr(e.i,e.j,e.k,0),dy=dr(e.i,e.j,e.k,1),dz=dr(e.i,e.j,e.k,2),l2=dx*dx+dy*dy+dz*dz,L=std::sqrt(l2);if(L<1.e-14)continue;
            const double tau[3]={dx/L,dy/L,dz/L};for(int a=0;a<3;++a)for(int c=0;c<3;++c)A[a][c]+=tau[a]*tau[c];terms.push_back({e,{dx/l2,dy/l2,dz/l2}});}
        std::map<I64,Vec3>row;if(terms.empty())return row;double inv[3][3];for(int c=0;c<3;++c){double b[3]={0,0,0},x[3];b[c]=1;if(!solve3(A,b,x))return {};for(int a=0;a<3;++a)inv[a][c]=x[a];}
        for(const auto&term:terms){Vec3 v{};for(int a=0;a<3;++a)for(int c=0;c<3;++c)v[a]+=inv[a][c]*term.q[c];add_jc(row,term.e,v);}return row;};
    for(const auto&p:t.physical_patches){if(p.this_rank!=rank||p.bc_name!="Pole")continue;const int dir=std::abs(p.direction);if(dir!=1&&dir!=2)continue;const int ib=p.this_block;const bool high=p.direction>0;const Box3&node=p.this_box_node;
        if(dir==1){const int ic=high?node.lo.i-1:node.lo.i,it=high?ic:ic+1;for(int j=node.lo.j;j<node.hi.j-1;++j){std::vector<EntityKey>es;for(int k=node.lo.k;k<node.hi.k-1;++k){
            es.push_back(TOPO::make_edge(rank,ib,ic,j,k,EntityAxis::Xi));es.push_back(TOPO::make_edge(rank,ib,ic,j+1,k,EntityAxis::Xi));es.push_back(TOPO::make_edge(rank,ib,ic,j,k+1,EntityAxis::Xi));es.push_back(TOPO::make_edge(rank,ib,ic,j+1,k+1,EntityAxis::Xi));
            es.push_back(TOPO::make_edge(rank,ib,it,j,k,EntityAxis::Eta));es.push_back(TOPO::make_edge(rank,ib,it,j,k+1,EntityAxis::Eta));es.push_back(TOPO::make_edge(rank,ib,it,j,k,EntityAxis::Zeta));es.push_back(TOPO::make_edge(rank,ib,it,j+1,k,EntityAxis::Zeta));}
            const auto row=pole_row(ib,es);if(!row.empty())for(int k=node.lo.k;k<node.hi.k-1;++k)jc_rows[gid(t,TOPO::make_cell(rank,ib,ic,j,k))]=row;}}
        else{const int jc=high?node.lo.j-1:node.lo.j,jt=high?jc:jc+1;for(int i=node.lo.i;i<node.hi.i-1;++i){std::vector<EntityKey>es;for(int k=node.lo.k;k<node.hi.k-1;++k){
            es.push_back(TOPO::make_edge(rank,ib,i,jc,k,EntityAxis::Eta));es.push_back(TOPO::make_edge(rank,ib,i+1,jc,k,EntityAxis::Eta));es.push_back(TOPO::make_edge(rank,ib,i,jc,k+1,EntityAxis::Eta));es.push_back(TOPO::make_edge(rank,ib,i+1,jc,k+1,EntityAxis::Eta));
            es.push_back(TOPO::make_edge(rank,ib,i,jt,k,EntityAxis::Xi));es.push_back(TOPO::make_edge(rank,ib,i,jt,k+1,EntityAxis::Xi));es.push_back(TOPO::make_edge(rank,ib,i,jt,k,EntityAxis::Zeta));es.push_back(TOPO::make_edge(rank,ib,i+1,jt,k,EntityAxis::Zeta));}
            const auto row=pole_row(ib,es);if(!row.empty())for(int k=node.lo.k;k<node.hi.k-1;++k)jc_rows[gid(t,TOPO::make_cell(rank,ib,i,jc,k))]=row;}}}
    std::vector<I64> jcg,jcoff{0},jce;std::vector<double> jcw;
    for(const auto&r:cells){jcg.push_back(r.gid);for(const auto&[edge,v]:jc_rows[r.gid]){if(o.validate&&(!std::isfinite(v[0])||!std::isfinite(v[1])||!std::isfinite(v[2])))throw std::runtime_error("PostData non-finite J_edge_to_cell_cartesian weight");
        if(v[0]!=0||v[1]!=0||v[2]!=0){jce.push_back(edge);jcw.insert(jcw.end(),v.begin(),v.end());}}jcoff.push_back(jce.size());}
    w.WriteSection("JedgeJcell_global_id",1,span(jcg));w.WriteSection("JedgeJcell_offsets",1,span(jcoff));
    w.WriteSection("JedgeJcell_edge_ids",1,span(jce));w.WriteSection("JedgeJcell_weights",3,span(jcw));

    if(o.validate){
        auto axis_field=[](const char *prefix,EntityAxis a){return std::string(prefix)+(a==EntityAxis::Xi?"xi":a==EntityAxis::Eta?"eta":"zeta");};
        const std::size_t nf=face_count(t),ne=edge_count(t),ncol=static_cast<std::size_t>(storage.global_column_count);if(nf>static_cast<std::size_t>(std::numeric_limits<int>::max())||ne>static_cast<std::size_t>(std::numeric_limits<int>::max())||ncol>static_cast<std::size_t>(std::numeric_limits<int>::max()))throw std::overflow_error("PostData DEC validation exceeds MPI int range");
        std::vector<double> lb(nf,0.0),gb(nf,0.0),lj(ne,0.0),gj(ne,0.0);std::vector<int> lbc(nf,0),gbc(nf,0),ljc(ne,0),gjc(ne,0);
        for(const auto&r:owner_rows(t,EntityDim::Face,rank)){const auto&f=fields.field(axis_field("B_",r.key.axis),r.key.block);lb[r.gid]=t.sign_to_owner(r.key)*f(r.key.i,r.key.j,r.key.k,0);lbc[r.gid]=1;}
        for(const auto&r:owner_rows(t,EntityDim::Edge,rank)){const auto&j=fields.field(axis_field("J_",r.key.axis),r.key.block);lj[r.gid]=t.sign_to_owner(r.key)*j(r.key.i,r.key.j,r.key.k,0);ljc[r.gid]=1;}
        MPI_Allreduce(lb.data(),gb.data(),static_cast<int>(nf),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(lbc.data(),gbc.data(),static_cast<int>(nf),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(lj.data(),gj.data(),static_cast<int>(ne),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(ljc.data(),gjc.data(),static_cast<int>(ne),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        for(std::size_t q=0;q<nf;++q)if(gbc[q]!=1)throw std::runtime_error("PostData B_face_to_J_edge validation missing/duplicate Face owner");
        for(std::size_t q=0;q<ne;++q)if(gjc[q]!=1)throw std::runtime_error("PostData B_face_to_J_edge validation missing/duplicate Edge owner");
        double local_abs=0.0,local_scale=0.0;
        for(const auto&[e,qkey]:t.edges.local_to_qkey){(void)qkey;if(e.rank!=rank)continue;const auto&j=fields.field(axis_field("J_",e.axis),e.block);const double oriented=t.sign_to_owner(e)*j(e.i,e.j,e.k,0),owner=gj[gid(t,e)];local_abs=std::max(local_abs,std::abs(oriented-owner));local_scale=std::max(local_scale,std::max(std::abs(oriented),std::abs(owner)));}
        for(const auto&[f,qkey]:t.faces.local_to_qkey){(void)qkey;if(f.rank!=rank)continue;const auto&b=fields.field(axis_field("B_",f.axis),f.block);const double oriented=t.sign_to_owner(f)*b(f.i,f.j,f.k,0),owner=gb[gid(t,f)];local_abs=std::max(local_abs,std::abs(oriented-owner));local_scale=std::max(local_scale,std::max(std::abs(oriented),std::abs(owner)));}
        double global_abs=0.0,global_scale=0.0;MPI_Allreduce(&local_abs,&global_abs,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);MPI_Allreduce(&local_scale,&global_scale,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        if(global_abs>512.0*std::numeric_limits<double>::epsilon()*std::max(1.0,global_scale))throw std::runtime_error("PostData owner/alias orientation validation failed: max_abs="+std::to_string(global_abs));
        std::vector<double> lbx(ncol,0.0),gbx(ncol,0.0);std::vector<int> lbxc(ncol,0),gbxc(ncol,0);
        for(std::size_t q=0;q<nf;++q){lbx[q]=lb[q];lbxc[q]=lbc[q];}
        for(const auto &entry:storage.entries){
            if((entry.flags&FaceStorageGhost)==0)continue;
            const auto &b=fields.field(face_field_name(entry.key.axis),entry.key.block);
            lbx[entry.column_id]=b(entry.key.i,entry.key.j,entry.key.k,0);lbxc[entry.column_id]=1;
        }
        MPI_Allreduce(lbx.data(),gbx.data(),static_cast<int>(ncol),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(lbxc.data(),gbxc.data(),static_cast<int>(ncol),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        for(std::size_t q=0;q<ncol;++q)if(gbxc[q]!=1)throw std::runtime_error("PostData B_face_to_J_edge validation missing/duplicate extended Bface column");
        local_abs=0.0;local_scale=0.0;I64 local_worst_gid=-1;double local_worst_value=0.0,local_worst_actual=0.0;
        for(const auto&r:owner_rows(t,EntityDim::Edge,rank)){double value=0.0;const auto it=bj_rows.find(r.gid);if(it!=bj_rows.end())for(const auto&[col,weight]:it->second)value+=weight*gbx[col];const double error=std::abs(value-gj[r.gid]);if(error>local_abs){local_abs=error;local_worst_gid=r.gid;local_worst_value=value;local_worst_actual=gj[r.gid];}local_scale=std::max(local_scale,std::max(std::abs(value),std::abs(gj[r.gid])));}
        MPI_Allreduce(&local_abs,&global_abs,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);MPI_Allreduce(&local_scale,&global_scale,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        if(global_abs>2048.0*std::numeric_limits<double>::epsilon()*std::max(1.0,global_scale)){
            struct MaxLoc{double error;int rank;} local_max{local_abs,rank},global_max{};
            MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
            I64 worst_gid=local_worst_gid;double detail[2]={local_worst_value,local_worst_actual};
            MPI_Bcast(&worst_gid,1,MPI_INT64_T,global_max.rank,MPI_COMM_WORLD);
            MPI_Bcast(detail,2,MPI_DOUBLE,global_max.rank,MPI_COMM_WORLD);
            throw std::runtime_error("PostData B_face_to_J_edge does not reproduce synchronized solver Jedge: max_abs="+std::to_string(global_abs)+", edge_gid="+std::to_string(worst_gid)+", recomputed="+std::to_string(detail[0])+", solver="+std::to_string(detail[1])+", singular="+(singular_rows.count(worst_gid)?"true":"false"));
        }
        local_abs=0.0;local_scale=0.0;
        for(const auto&r:cells){Vec3 value{};for(const auto&[edge,weight]:jc_rows[r.gid])for(int a=0;a<3;++a)value[a]+=weight[a]*gj[edge];const auto&j=fields.field("J_cell",r.key.block);for(int a=0;a<3;++a){const double actual=j(r.key.i,r.key.j,r.key.k,a);local_abs=std::max(local_abs,std::abs(value[a]-actual));local_scale=std::max(local_scale,std::max(std::abs(value[a]),std::abs(actual)));}}
        MPI_Allreduce(&local_abs,&global_abs,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);MPI_Allreduce(&local_scale,&global_scale,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        if(global_abs>2048.0*std::numeric_limits<double>::epsilon()*std::max(1.0,global_scale))throw std::runtime_error("PostData J_edge_to_cell_cartesian does not reproduce solver Jcell: max_abs="+std::to_string(global_abs));
    }
    w.Close();
}

void write_fields(const std::filesystem::path&path,FileType type,const Grid&grid,const TOPO::Topology&t,const Field&fields,int rank,const std::vector<std::string>&names,const WriteOptions&o)
{
    PostBinaryWriter w(path,type,o.case_uuid,o.mesh_uuid);std::set<std::string> seen;
    int field_index=0;
    for(const auto&name:names){if(!seen.insert(name).second)continue;if(!fields.has_field(name))throw std::invalid_argument("PostData unknown field: "+name);const auto&d=fields.descriptor(name);
        std::vector<I64> ids;std::vector<double> values;for(int ib=0;ib<grid.nblock;++ib){const auto&fb=fields.field(name,ib);if(!fb.is_allocated())continue;const Int3 lo=fb.inner_lo(),hi=fb.inner_hi();for(int k=lo.k;k<hi.k;++k)for(int j=lo.j;j<hi.j;++j)for(int i=lo.i;i<hi.i;++i){EntityAxis a=EntityAxis::None;EntityDim dim=EntityDim::Cell;
            switch(d.location){case StaggerLocation::Node:dim=EntityDim::Node;break;case StaggerLocation::EdgeXi:dim=EntityDim::Edge;a=EntityAxis::Xi;break;case StaggerLocation::EdgeEt:dim=EntityDim::Edge;a=EntityAxis::Eta;break;case StaggerLocation::EdgeZe:dim=EntityDim::Edge;a=EntityAxis::Zeta;break;case StaggerLocation::FaceXi:dim=EntityDim::Face;a=EntityAxis::Xi;break;case StaggerLocation::FaceEt:dim=EntityDim::Face;a=EntityAxis::Eta;break;case StaggerLocation::FaceZe:dim=EntityDim::Face;a=EntityAxis::Zeta;break;default:break;}
            EntityKey e=dim==EntityDim::Node?TOPO::make_node(rank,ib,i,j,k):dim==EntityDim::Cell?TOPO::make_cell(rank,ib,i,j,k):dim==EntityDim::Edge?TOPO::make_edge(rank,ib,i,j,k,a):TOPO::make_face(rank,ib,i,j,k,a);if((dim==EntityDim::Edge||dim==EntityDim::Face||dim==EntityDim::Node)&&!t.is_owner(e))continue;ids.push_back(gid(t,e));const int os=(dim==EntityDim::Edge||dim==EntityDim::Face)?t.sign_to_owner(e):1;for(int m=0;m<d.ncomp;++m){double v=fb(i,j,k,m)*os;if(o.validate&&!std::isfinite(v))throw std::runtime_error("PostData non-finite field value: "+name);values.push_back(v);}}}
        std::ostringstream key;key<<"field_"<<std::setw(4)<<std::setfill('0')<<field_index++;
        std::vector<I64> meta{static_cast<I64>(d.location),d.ncomp,static_cast<I64>(d.value_kind)};
        w.WriteSection(key.str()+"_meta",3,span(meta));w.WriteSection(key.str()+"_ids",1,span(ids));w.WriteSection(key.str()+"_values",d.ncomp,span(values));}
    w.Close();
}

std::string json_escape(const std::string&s){std::string o;for(char c:s){if(c=='"'||c=='\\')o+='\\';o+=c;}return o;}
void write_manifest(const std::filesystem::path&dir,const Grid&grid,const Field&fields,int size,int total_blocks,const WriteOptions&o)
{
    std::ofstream f(dir/"manifest.json",std::ios::trunc);if(!f)throw std::runtime_error("cannot write manifest.json");
    const auto jr=o.normalization.find("current_density_ref"),lr=o.normalization.find("length_ref");
    const bool has_jedge_ref=jr!=o.normalization.end()&&lr!=o.normalization.end();
    f<<"{\n  \"format_name\": \"MPCNS_PostData\",\n  \"format_version\": 3,\n  \"case_uuid\": \""<<uuid_string(o.case_uuid)<<"\",\n  \"mesh_uuid\": \""<<uuid_string(o.mesh_uuid)<<"\",\n  \"endianness\": \"little\",\n  \"float_type\": \"float64\",\n  \"index_type\": \"int64\",\n  \"dimension\": "<<grid.dimension<<",\n  \"number_of_blocks\": "<<total_blocks<<",\n  \"number_of_ranks\": "<<size<<",\n  \"array_order\": \"C\",\n  \"logical_index_order\": \"i-fastest\",\n  \"linear_index\": \"i + ni * (j + nj * k)\",\n  \"face_magnetic_semantics\": \"oriented_face_2form_flux\",\n  \"Bface_operator_column_space\": {\"kind\":\"quotient_faces_plus_restart_ghost_storage\",\"quotient_id_range\":\"[0,Nface)\",\"ghost_id_range\":\"[Nface,Bstore_global_column_count)\",\"topology_sections\":[\"Bstore_meta\",\"Bstore_global_id\",\"Bstore_quotient_id\",\"Bstore_address\",\"Bstore_sign\",\"Bstore_owner\",\"Bstore_flags\"],\"address_components\":[\"rank\",\"rank_local_block\",\"restart_location_code\",\"i\",\"j\",\"k\"],\"quotient_input_rule\":\"owner_oriented_value; local_value * Bstore_sign\",\"ghost_input_rule\":\"raw solver-native value from the matching restart storage slot; no orientation transform\",\"decomposition_scope\":\"ghost IDs are valid for dynamic files with the same rank/block decomposition and mesh_uuid\"},\n  \"dec_current_semantics\": {\"stored_quantity\":\"J_dot_dr\",\"form_degree\":1,\"variance\":\"covariant\",\"orientation_aware\":true,\"magnetic_source\":\"induced_B_only\",\"excludes_fields\":[\"Badd_xi\",\"Badd_eta\",\"Badd_zeta\"],\"physical_reference_expression\":\"current_density_ref * length_ref\",\"physical_reference_value\":";
    if(has_jedge_ref)f<<std::setprecision(17)<<jr->second*lr->second;else f<<"null";f<<"},\n  \"normalization\": {";
    bool first_item=true;for(const auto&[k,v]:o.normalization){if(!first_item)f<<',';first_item=false;f<<'"'<<json_escape(k)<<"\":"<<std::setprecision(17)<<v;}f<<"},\n  \"physical_constants\": {";
    first_item=true;for(const auto&[k,v]:o.physical_constants){if(!first_item)f<<',';first_item=false;f<<'"'<<json_escape(k)<<"\":"<<std::setprecision(17)<<v;}f<<"},\n  \"species\": [";
    for(std::size_t i=0;i<o.species.size();++i){if(i)f<<',';f<<'"'<<json_escape(o.species[i])<<'"';}f<<"],\n  \"files\": {\n";
    const char*stems[]={"geometry","topology","reconstruction","constant_field"};for(int s=0;s<4;++s){f<<"    \""<<stems[s]<<"\": [";for(int r=0;r<size;++r){if(r)f<<", ";f<<'"'<<chunk_name(stems[s],r)<<'"';}f<<"]"<<(s==3?'\n':',') ;}
    f<<"  },\n  \"existing_dynamic_data\": {\n"
       "    \"writer\": \"IOModule::WriteRestartBinFile\",\n"
       "    \"format_name\": \"MPCNS_Restart\",\n"
       "    \"magic\": \"MPCNSRST\",\n"
       "    \"format_version\": 1,\n"
       "    \"path_pattern\": \""<<json_escape(o.existing_flow_path_pattern)<<"\",\n"
       "    \"number_of_rank_files\": "<<size<<",\n"
       "    \"overwrite_semantics\": \"latest checkpoint per rank\",\n"
       "    \"optional_dec_jedge\": {\"compile_option\":\"OUTPUT_DEC_JEDGE\",\"runtime_parameter\":\"output_dec_jedge\",\"default_enabled\":false,\"cadence\":\"same_as_Bface_restart\",\"snapshot_phase\":\"completed time step after implicit convergence, owner/alias reconciliation, physical boundary handling, and Jedge sync\"},\n"
       "    \"contains_ghost_layers\": true,\n"
       "    \"header_layout\": [\"magic[8]\",\"version:int32\",\"step:int32\",\"time:float64\",\"nblock:int32\",\"nfield:int32\"],\n"
       "    \"field_header_layout\": [\"name_length:int32\",\"name_bytes\",\"location:int32\",\"components:int32\",\"nghost:int32\"],\n"
       "    \"block_header_layout\": [\"lo_i:int32\",\"lo_j:int32\",\"lo_k:int32\",\"hi_i:int32\",\"hi_j:int32\",\"hi_k:int32\",\"active:int32\"],\n"
       "    \"value_type\": \"float64\",\n"
       "    \"value_loop_order\": \"i, j, k, component (component fastest)\",\n"
       "    \"block_index_scope\": \"rank-local; use topology Bstore_address for Bface real/ghost storage and block maps for other global entities\",\n"
       "    \"oriented_field_rule\": \"quotient global owner value = local value * Bstore_sign; auxiliary ghost value is read raw with sign 1\",\n"
       "    \"inactive_semantics\": \"a field is intentionally inactive when its descriptor physics does not match the block physics; do not treat that as missing data for that domain\",\n"
       "    \"fields\": [";
    first_item=true;for(const auto&name:o.existing_flow_fields){if(!fields.has_field(name))continue;const auto&d=fields.descriptor(name);if(!first_item)f<<',';first_item=false;f<<"{\"name\":\""<<json_escape(name)<<"\",\"location\":\""<<location_name(d.location)<<"\",\"location_code\":"<<restart_location_code(d.location)<<",\"components\":"<<d.ncomp<<",\"nghost\":"<<d.nghost<<",\"value_kind\":\""<<field_value_kind_name(d.value_kind)<<"\",\"physics_domain\":\""<<json_escape(d.physics.empty()?"all":d.physics)<<"\"}";}f<<"]\n  },\n"
       "  \"block_physics_codes\": {\"0\":\"unknown\",\"1\":\"Fluid\",\"2\":\"Solid\"},\n"
       "  \"cell_flag_bits\": {\"fluid\":1,\"solid\":2,\"unknown_physics\":4},\n"
       "  \"Bstore_flag_bits\": {\"ghost\":1,\"physical_boundary\":2,\"interface\":4,\"coupling\":8,\"parallel\":16},\n"
       "  \"physical_patch_meta_components\": [\"rank\",\"rank_local_block\",\"bc_id\",\"direction\",\"lo_i\",\"lo_j\",\"lo_k\",\"hi_i\",\"hi_j\",\"hi_k\",\"name_byte_offset\",\"name_byte_length\"],\n"
       "  \"interface_patch_meta_components\": [\"rank\",\"neighbor_rank\",\"rank_local_block\",\"neighbor_rank_local_block\",\"direction\",\"neighbor_direction\",\"patch_kind\",\"is_coupling\",\"perm_0\",\"perm_1\",\"perm_2\",\"sign_0\",\"sign_1\",\"sign_2\",\"offset_i\",\"offset_j\",\"offset_k\",\"this_lo_i\",\"this_lo_j\",\"this_lo_k\",\"this_hi_i\",\"this_hi_j\",\"this_hi_k\",\"neighbor_lo_i\",\"neighbor_lo_j\",\"neighbor_lo_k\",\"neighbor_hi_i\",\"neighbor_hi_j\",\"neighbor_hi_k\"],\n"
       "  \"operators\": [{\"name\":\"B_face_to_cell_cartesian\",\"input_location\":\"face\",\"input_value_kind\":\"face_2form\",\"input_components\":1,\"output_location\":\"cell\",\"output_value_kind\":\"cartesian_vector\",\"output_components\":3,\"operator_version\":1},{\"name\":\"B_face_to_J_edge\",\"input_location\":\"extended_face_storage\",\"input_value_kind\":\"oriented_quotient_face_or_raw_restart_ghost_flux\",\"output_location\":\"edge\",\"output_value_kind\":\"edge_covariant_1form\",\"output_components\":1,\"solver_normalized_formula\":\"M1_inverse * D1_transpose * M2\",\"physical_formula\":\"(1/mu_m) * M1_inverse * D1_transpose * M2\",\"includes_physical_boundary_ghosts\":true,\"includes_interface_and_coupling_ghosts\":true,\"weights_source\":\"solver_final_assembly\",\"operator_version\":2},{\"name\":\"J_edge_to_cell_cartesian\",\"input_location\":\"edge\",\"input_value_kind\":\"edge_covariant_1form\",\"output_location\":\"cell\",\"output_value_kind\":\"cartesian_vector\",\"output_components\":3,\"includes_pole_ring_override\":true,\"weights_source\":\"solver_final_assembly\",\"operator_version\":1},{\"name\":\"cell_scalar_to_node\",\"operator_version\":1}],\n  \"fields\": [";
    bool first=true;int index=0;for(const auto&name:o.constant_fields){if(!fields.has_field(name))continue;const auto&d=fields.descriptor(name);if(!first)f<<',';first=false;std::ostringstream key;key<<"field_"<<std::setw(4)<<std::setfill('0')<<index++;f<<"{\"name\":\""<<json_escape(name)<<"\",\"location\":\""<<location_name(d.location)<<"\",\"components\":"<<d.ncomp<<",\"section_prefix\":\""<<key.str()<<"\"}";}f<<"]\n}\n";
    if(!f)throw std::runtime_error("manifest.json write failed");
}
}

PostDataWriter::PostDataWriter(const Grid&g,const TOPO::Topology&t,const Field&f,const RunData&r,int rank,const METRIC::SingularEdgeRegistry*s)
    :grid_(g),topology_(t),fields_(f),run_data_(r),singular_edges_(s),rank_(rank){}

void PostDataWriter::WriteStaticData(const std::filesystem::path&dir,WriteOptions options) const
{
    options=normalize_options(topology_,std::move(options));if(options.validate)validate_global_ids_and_orientation(topology_);std::filesystem::create_directories(dir);
    if(options.validate){validate_entity_partition(owner_rows(topology_,EntityDim::Node,rank_),node_count(topology_),"node");validate_entity_partition(owner_rows(topology_,EntityDim::Edge,rank_),edge_count(topology_),"edge");validate_entity_partition(owner_rows(topology_,EntityDim::Face,rank_),face_count(topology_),"face");validate_entity_partition(owner_rows(topology_,EntityDim::Cell,rank_),cell_count(topology_),"cell");}
    const auto face_storage=build_face_storage_catalog(topology_,fields_,rank_);
    write_geometry(dir/chunk_name("geometry",rank_),grid_,topology_,fields_,face_storage,singular_edges_,rank_,options);
    write_topology(dir/chunk_name("topology",rank_),grid_,topology_,face_storage,rank_,options);
    write_reconstruction(dir/chunk_name("reconstruction",rank_),grid_,topology_,fields_,face_storage,singular_edges_,rank_,options);
    write_fields(dir/chunk_name("constant_field",rank_),FileType::ConstantField,grid_,topology_,fields_,rank_,options.constant_fields,options);
    PARALLEL::mpi_barrier();int size=1;PARALLEL::mpi_size(&size);int local_blocks=grid_.nblock,total_blocks=0;MPI_Allreduce(&local_blocks,&total_blocks,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);if(rank_==0)write_manifest(dir,grid_,fields_,size,total_blocks,options);PARALLEL::mpi_barrier();
}

} // namespace POST
