#include "6_io/PostDataWriter.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
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

void write_topology(const std::filesystem::path &path,const Grid&grid,const TOPO::Topology&t,int rank,const WriteOptions&o)
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
    w.WriteSection("block_connections",15,span(conn));w.Close();
}

void write_reconstruction(const std::filesystem::path&path,const Grid&grid,const TOPO::Topology&t,const Field&fields,int rank,const WriteOptions&o)
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
    w.WriteSection("NodeScalar_global_id",1,span(ng));w.WriteSection("NodeScalar_offsets",1,span(noff));w.WriteSection("NodeScalar_cell_ids",1,span(nc));w.WriteSection("NodeScalar_weights",1,span(nw));w.Close();
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
    f<<"{\n  \"format_name\": \"MPCNS_PostData\",\n  \"format_version\": 1,\n  \"case_uuid\": \""<<uuid_string(o.case_uuid)<<"\",\n  \"mesh_uuid\": \""<<uuid_string(o.mesh_uuid)<<"\",\n  \"endianness\": \"little\",\n  \"float_type\": \"float64\",\n  \"index_type\": \"int64\",\n  \"dimension\": "<<grid.dimension<<",\n  \"number_of_blocks\": "<<total_blocks<<",\n  \"number_of_ranks\": "<<size<<",\n  \"array_order\": \"C\",\n  \"logical_index_order\": \"i-fastest\",\n  \"linear_index\": \"i + ni * (j + nj * k)\",\n  \"face_magnetic_semantics\": \"oriented_face_2form_flux\",\n  \"normalization\": {";
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
       "    \"contains_ghost_layers\": true,\n"
       "    \"header_layout\": [\"magic[8]\",\"version:int32\",\"step:int32\",\"time:float64\",\"nblock:int32\",\"nfield:int32\"],\n"
       "    \"field_header_layout\": [\"name_length:int32\",\"name_bytes\",\"location:int32\",\"components:int32\",\"nghost:int32\"],\n"
       "    \"block_header_layout\": [\"lo_i:int32\",\"lo_j:int32\",\"lo_k:int32\",\"hi_i:int32\",\"hi_j:int32\",\"hi_k:int32\",\"active:int32\"],\n"
       "    \"value_type\": \"float64\",\n"
       "    \"value_loop_order\": \"i, j, k, component (component fastest)\",\n"
       "    \"block_index_scope\": \"rank-local; use topology chunk block maps for global entities\",\n"
       "    \"oriented_field_rule\": \"global_owner_value = local_value / local_orientation_sign\",\n"
       "    \"inactive_semantics\": \"a field is intentionally inactive when its descriptor physics does not match the block physics; do not treat that as missing data for that domain\",\n"
       "    \"fields\": [";
    first_item=true;for(const auto&name:o.existing_flow_fields){if(!fields.has_field(name))continue;const auto&d=fields.descriptor(name);if(!first_item)f<<',';first_item=false;f<<"{\"name\":\""<<json_escape(name)<<"\",\"location\":\""<<location_name(d.location)<<"\",\"location_code\":"<<restart_location_code(d.location)<<",\"components\":"<<d.ncomp<<",\"nghost\":"<<d.nghost<<",\"value_kind\":\""<<field_value_kind_name(d.value_kind)<<"\",\"physics_domain\":\""<<json_escape(d.physics.empty()?"all":d.physics)<<"\"}";}f<<"]\n  },\n"
       "  \"block_physics_codes\": {\"0\":\"unknown\",\"1\":\"Fluid\",\"2\":\"Solid\"},\n"
       "  \"cell_flag_bits\": {\"fluid\":1,\"solid\":2,\"unknown_physics\":4},\n"
       "  \"operators\": [{\"name\":\"B_face_to_cell_cartesian\",\"input_location\":\"face\",\"input_value_kind\":\"face_2form\",\"input_components\":1,\"output_location\":\"cell\",\"output_value_kind\":\"cartesian_vector\",\"output_components\":3,\"operator_version\":1},{\"name\":\"cell_scalar_to_node\",\"operator_version\":1}],\n  \"fields\": [";
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
    write_geometry(dir/chunk_name("geometry",rank_),grid_,topology_,fields_,singular_edges_,rank_,options);
    write_topology(dir/chunk_name("topology",rank_),grid_,topology_,rank_,options);
    write_reconstruction(dir/chunk_name("reconstruction",rank_),grid_,topology_,fields_,rank_,options);
    write_fields(dir/chunk_name("constant_field",rank_),FileType::ConstantField,grid_,topology_,fields_,rank_,options.constant_fields,options);
    PARALLEL::mpi_barrier();int size=1;PARALLEL::mpi_size(&size);int local_blocks=grid_.nblock,total_blocks=0;MPI_Allreduce(&local_blocks,&total_blocks,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);if(rank_==0)write_manifest(dir,grid_,fields_,size,total_blocks,options);PARALLEL::mpi_barrier();
}

} // namespace POST
