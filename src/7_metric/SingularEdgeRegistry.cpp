#include "7_metric/SingularEdgeRegistry.h"

#include "0_basic/Error.h"
#include "0_basic/MPI_WRAPPER.h"
#include "3_field/Field.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/LocalIncidence.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <set>
#include <sstream>

namespace
{
using TOPO::EntityAxis;
using TOPO::EntityDim;
using TOPO::EntityKey;

bool contains_inner(const FieldBlock &f, int i, int j, int k)
{
    const Int3 lo = f.inner_lo(), hi = f.inner_hi();
    return i >= lo.i && i < hi.i && j >= lo.j && j < hi.j && k >= lo.k && k < hi.k;
}

std::vector<EntityKey> incident_cells(const EntityKey &e)
{
    std::vector<EntityKey> out;
    const int a = static_cast<int>(e.axis);
    for (int s0 : {-1, 0})
        for (int s1 : {-1, 0})
        {
            int i=e.i, j=e.j, k=e.k;
            if (a == 0) { j += s0; k += s1; }
            if (a == 1) { i += s0; k += s1; }
            if (a == 2) { i += s0; j += s1; }
            out.push_back(TOPO::make_cell(e.rank, e.block, i, j, k));
        }
    return out;
}

std::vector<EntityKey> incident_faces(const EntityKey &e)
{
    std::vector<EntityKey> out;
    if (e.axis == EntityAxis::Xi)
    {
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Eta));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j,e.k-1,EntityAxis::Eta));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Zeta));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j-1,e.k,EntityAxis::Zeta));
    }
    else if (e.axis == EntityAxis::Eta)
    {
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Xi));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j,e.k-1,EntityAxis::Xi));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Zeta));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i-1,e.j,e.k,EntityAxis::Zeta));
    }
    else
    {
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Xi));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j-1,e.k,EntityAxis::Xi));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i,e.j,e.k,EntityAxis::Eta));
        out.push_back(TOPO::make_face(e.rank,e.block,e.i-1,e.j,e.k,EntityAxis::Eta));
    }
    return out;
}

const char *axis_suffix(EntityAxis a)
{
    return a == EntityAxis::Xi ? "xi" : (a == EntityAxis::Eta ? "eta" : "zeta");
}
} // namespace

namespace METRIC
{
void SingularEdgeRegistry::build(const TOPO::Topology &topology, Field &fields, Grid &grid, int rank)
{
    rank_ = rank;
    entries_.clear(); curl_entries_.clear(); gid_to_index_.clear();
    const int nglobal = topology.edges.n_global_owner;
    std::vector<double> local_length(nglobal,0.0), local_cell_measure(nglobal,0.0);
    std::vector<double> local_face_measure(nglobal,0.0);
    std::vector<double> local_dr(3*nglobal,0.0);
    std::vector<double> local_center_records;
    std::vector<int> local_length_n(nglobal,0), local_cells(nglobal,0), local_faces(nglobal,0);

    for (const auto &cls : topology.edges.classes)
    {
        // Build the physical-cell polygon for every shared quotient edge.
        // Variable-valence singular edges remain a separate subset for EMF
        // assembly and output, while regular seam edges also use this polygon
        // for the non-orthogonal curl/current operator.
        if (cls.members.size() < 2 || cls.global_id < 0)
            continue;

        SingularPhysicalEdge rec;
        rec.global_id = cls.global_id;
        rec.owner = cls.owner.entity;
        for (const auto &m : cls.members)
            rec.aliases.push_back({m.entity, m.orient_sign, m.is_owner});

        std::set<EntityKey> seen_cells, seen_faces;
        for (const auto &a : rec.aliases)
        {
            const EntityKey &e = a.edge;
            if (e.rank != rank || e.block < 0 || e.block >= fields.num_blocks()) continue;

            FieldBlock &dl = fields.field(std::string("dl_") + axis_suffix(e.axis), e.block);
            if (contains_inner(dl,e.i,e.j,e.k))
            {
                const double L = dl(e.i,e.j,e.k,0);
                if (std::isfinite(L) && L > 0.0) { local_length[rec.global_id] += L; ++local_length_n[rec.global_id]; }
                FieldBlock &dr = fields.field(std::string("dr_") + axis_suffix(e.axis), e.block);
                for (int q=0;q<3;++q)
                    local_dr[3*rec.global_id+q] += a.orientation*dr(e.i,e.j,e.k,q);
            }

            FieldBlock &jac = fields.field("Jac", e.block);
            for (const EntityKey &c : incident_cells(e))
                if (seen_cells.insert(c).second && contains_inner(jac,c.i,c.j,c.k))
                {
                    const double v = std::abs(jac(c.i,c.j,c.k,0));
                    if (std::isfinite(v) && v > 1.e-30)
                    {
                        rec.local_incident_cells.push_back({c,v,0.0,e,a.orientation});
                        local_cell_measure[rec.global_id] += v; ++local_cells[rec.global_id];
                        Block &blk=grid.grids(c.block);
                        local_center_records.push_back(static_cast<double>(rec.global_id));
                        local_center_records.push_back(blk.dual_x(c.i+1,c.j+1,c.k+1));
                        local_center_records.push_back(blk.dual_y(c.i+1,c.j+1,c.k+1));
                        local_center_records.push_back(blk.dual_z(c.i+1,c.j+1,c.k+1));
                    }
                }

            for (const EntityKey &f : incident_faces(e))
                if (seen_faces.insert(f).second)
                {
                    // Count and assemble each quotient (physical) face once.
                    // Block-local aliases of the same interface face are not
                    // distinct sectors around the physical edge.
                    if (topology.faces.local_to_qkey.find(f) ==
                        topology.faces.local_to_qkey.end()) continue;
                    if (!topology.is_owner(f)) continue;
                    FieldBlock &area = fields.field(std::string("Area_") + axis_suffix(f.axis), f.block);
                    if (!contains_inner(area,f.i,f.j,f.k)) continue;
                    const double av = area(f.i,f.j,f.k,0);
                    if (std::isfinite(av) && av > 1.e-30)
                    {
                        int local_incidence=0;
                        for(const auto &term:TOPO::boundary_of_face(f))
                            if(term.entity==e) { local_incidence=term.sign; break; }
                        if(local_incidence==0)
                            ERROR::Abort("SingularEdgeRegistry: incident face does not contain its source edge");
                        METRIC::WeightedIncidentEntity inc{f,av,0.0,e,a.orientation};
                        // Registry canonical vectors/cochains use the quotient
                        // key orientation (the same convention as alias
                        // reductions), not the arbitrary local owner axis.
                        inc.entity_orientation=topology.face_qsign(f);
                        inc.quotient_incidence=inc.entity_orientation*local_incidence*
                                                topology.edge_qsign(e);
                        rec.local_incident_faces.push_back(inc);
                        local_face_measure[rec.global_id]+=av;
                        ++local_faces[rec.global_id];
                    }
                }
        }
        gid_to_index_[rec.global_id] = entries_.size(); entries_.push_back(std::move(rec));
    }

    std::vector<double> gl(nglobal), gv(nglobal);
    std::vector<double> gfm(nglobal);
    std::vector<double> gdr(3*nglobal);
    std::vector<int> gln(nglobal), gc(nglobal), gf(nglobal);
    MPI_Allreduce(local_length.data(),gl.data(),nglobal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_cell_measure.data(),gv.data(),nglobal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_face_measure.data(),gfm.data(),nglobal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_dr.data(),gdr.data(),3*nglobal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_length_n.data(),gln.data(),nglobal,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_cells.data(),gc.data(),nglobal,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_faces.data(),gf.data(),nglobal,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    int local_record_values=static_cast<int>(local_center_records.size());
    int nrank=1;
    MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    std::vector<int> record_counts(nrank,0), record_offsets(nrank,0);
    MPI_Allgather(&local_record_values,1,MPI_INT,record_counts.data(),1,MPI_INT,MPI_COMM_WORLD);
    int total_record_values=0;
    for(int r=0;r<nrank;++r) { record_offsets[r]=total_record_values; total_record_values+=record_counts[r]; }
    std::vector<double> global_center_records(total_record_values);
    MPI_Allgatherv(local_center_records.data(),local_record_values,MPI_DOUBLE,
                   global_center_records.data(),record_counts.data(),record_offsets.data(),
                   MPI_DOUBLE,MPI_COMM_WORLD);
    std::unordered_map<int,std::vector<std::array<double,3>>> incident_centers;
    for(int n=0;n+3<total_record_values;n+=4)
        incident_centers[static_cast<int>(global_center_records[n])].push_back(
            {global_center_records[n+1],global_center_records[n+2],global_center_records[n+3]});

    for (std::size_t entry_index = 0; entry_index < entries_.size(); ++entry_index)
    {
        auto &e = entries_[entry_index];
        e.primal_length = gln[e.global_id] ? gl[e.global_id]/gln[e.global_id] : 0.0;
        // The dual face is the polygon through the centers of the real cells
        // incident on this quotient edge.  Sum(V/L) is a primal cross-section
        // measure and overestimates this polygon at extraordinary valence.
        e.dual_area=0.0;
        auto center_it=incident_centers.find(e.global_id);
        if(center_it!=incident_centers.end() && center_it->second.size()>=3)
        {
            auto points=center_it->second;
            std::array<double,3> center{{0.0,0.0,0.0}};
            for(const auto &p:points) for(int q=0;q<3;++q) center[q]+=p[q];
            for(int q=0;q<3;++q) center[q]/=points.size();
            std::array<double,3> tangent{{gdr[3*e.global_id],gdr[3*e.global_id+1],gdr[3*e.global_id+2]}};
            double tn=std::sqrt(tangent[0]*tangent[0]+tangent[1]*tangent[1]+tangent[2]*tangent[2]);
            if(tn>0.0) for(double &v:tangent) v/=tn;
            std::array<double,3> ref=std::abs(tangent[0])<0.9 ? std::array<double,3>{{1,0,0}} : std::array<double,3>{{0,1,0}};
            std::array<double,3> u{{tangent[1]*ref[2]-tangent[2]*ref[1],
                                    tangent[2]*ref[0]-tangent[0]*ref[2],
                                    tangent[0]*ref[1]-tangent[1]*ref[0]}};
            double un=std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
            if(un>0.0) for(double &v:u) v/=un;
            std::array<double,3> v{{tangent[1]*u[2]-tangent[2]*u[1],
                                    tangent[2]*u[0]-tangent[0]*u[2],
                                    tangent[0]*u[1]-tangent[1]*u[0]}};
            std::sort(points.begin(),points.end(),[&](const auto &a,const auto &b)
            {
                auto angle=[&](const auto &p)
                {
                    double du=0.0,dv=0.0;
                    for(int q=0;q<3;++q) { du+=(p[q]-center[q])*u[q]; dv+=(p[q]-center[q])*v[q]; }
                    return std::atan2(dv,du);
                };
                return angle(a)<angle(b);
            });
            e.ordered_cell_centers=points;
            for(auto &local:e.local_incident_cells)
            {
                Block &blk=grid.grids(local.entity.block);
                const std::array<double,3> p{{
                    blk.dual_x(local.entity.i+1,local.entity.j+1,local.entity.k+1),
                    blk.dual_y(local.entity.i+1,local.entity.j+1,local.entity.k+1),
                    blk.dual_z(local.entity.i+1,local.entity.j+1,local.entity.k+1)}};
                double best=std::numeric_limits<double>::max();
                for(std::size_t n=0;n<points.size();++n)
                {
                    const double dx=p[0]-points[n][0];
                    const double dy=p[1]-points[n][1];
                    const double dz=p[2]-points[n][2];
                    const double d2=dx*dx+dy*dy+dz*dz;
                    if(d2<best) { best=d2; local.sector_index=static_cast<int>(n); }
                }
            }
            std::array<double,3> area_vec{{0.0,0.0,0.0}};
            for(std::size_t n=0;n<points.size();++n)
            {
                const auto &a=points[n], &b=points[(n+1)%points.size()];
                area_vec[0]+=a[1]*b[2]-a[2]*b[1];
                area_vec[1]+=a[2]*b[0]-a[0]*b[2];
                area_vec[2]+=a[0]*b[1]-a[1]*b[0];
            }
            e.dual_area=0.5*std::sqrt(area_vec[0]*area_vec[0]+area_vec[1]*area_vec[1]+area_vec[2]*area_vec[2]);
        }
        if(!(e.dual_area>0.0))
            e.dual_area=e.primal_length>0.0 ? gv[e.global_id]/e.primal_length : 0.0;
        e.inverse_hodge = e.dual_area > 0.0 ? e.primal_length/e.dual_area : 0.0;
        if (gln[e.global_id] > 0)
            for (int q=0;q<3;++q) e.canonical_edge_vector[q]=gdr[3*e.global_id+q]/gln[e.global_id];
        e.global_valid_cell_count=gc[e.global_id]; e.global_valid_face_count=gf[e.global_id];
        for (auto &c : e.local_incident_cells) c.weight = gv[e.global_id] > 0.0 ? c.measure/gv[e.global_id] : 0.0;
        for (auto &f:e.local_incident_faces) f.weight = gfm[e.global_id] > 0.0 ? f.measure/gfm[e.global_id] : 0.0;

        const int physical_faces=gf[e.global_id];
        e.variable_valence=e.aliases.size()>=3 &&
                           (physical_faces==3 || physical_faces>4);
        if(e.variable_valence)
            // Publish the variable-valence scalar Hodge data to every local
            // singular alias. Regular seam metrics remain block-local.
            for (const auto &a:e.aliases)
            {
                const auto &x=a.edge;
                if (x.rank!=rank_ || x.block<0 || x.block>=fields.num_blocks()) continue;
                FieldBlock &astar=fields.field(std::string("Astar_")+axis_suffix(x.axis),x.block);
                FieldBlock &hinv=fields.field(std::string("Hodge_star_inverse_2form_to_1form_edge_")+axis_suffix(x.axis)+"_lumped",x.block);
                if (contains_inner(astar,x.i,x.j,x.k)) astar(x.i,x.j,x.k,0)=e.dual_area;
                if (contains_inner(hinv,x.i,x.j,x.k)) hinv(x.i,x.j,x.k,0)=e.inverse_hodge;
            }
    }

    curl_entries_=entries_;
    curl_entries_.erase(
        std::remove_if(curl_entries_.begin(),curl_entries_.end(),
                       [](const SingularPhysicalEdge &e)
                       {
                           return e.ordered_cell_centers.size()<3 ||
                                  e.global_valid_cell_count!=e.global_valid_face_count ||
                                  !(e.inverse_hodge>0.0) ||
                                  !std::isfinite(e.inverse_hodge);
                       }),
        curl_entries_.end());

    // Singularity is the valence of quotient faces incident on the physical
    // edge, not the number of block-local aliases. Three and greater-than-four
    // sectors use the variable-valence EMF/output path; four is a regular seam
    // edge and remains only in curl_entries_.
    entries_.erase(std::remove_if(entries_.begin(), entries_.end(),
                                  [&](const SingularPhysicalEdge &e)
                                  {
                                      const int physical_faces = gf[e.global_id];
                                      return e.aliases.size()<3 ||
                                             (physical_faces != 3 && physical_faces <= 4);
                                  }),
                   entries_.end());
    gid_to_index_.clear();
    for(std::size_t n=0;n<entries_.size();++n)
        gid_to_index_[entries_[n].global_id]=n;

}

void SingularEdgeRegistry::validate_or_abort() const
{
    for (const auto &e : entries_)
        if (e.aliases.empty() || e.global_id < 0 || !(e.primal_length > 0.0) ||
            !(e.dual_area > 0.0) || !std::isfinite(e.inverse_hodge) || e.global_valid_cell_count < 1)
        {
            std::ostringstream os; os << "SingularEdgeRegistry invalid entry gid=" << e.global_id
                << " aliases=" << e.aliases.size() << " cells=" << e.global_valid_cell_count
                << " length=" << e.primal_length << " dual_area=" << e.dual_area;
            ERROR::Abort(os.str());
        }
}

const SingularPhysicalEdge *SingularEdgeRegistry::find(int gid) const
{
    auto it=gid_to_index_.find(gid); return it==gid_to_index_.end()?nullptr:&entries_[it->second];
}

void SingularEdgeRegistry::assemble_cell_field_to_local_owners(
    Field &fields, const std::string &field_name,
    const CellContribution &contribution) const
{
    if (entries_.empty()) return;
    const int fid = fields.field_id(field_name);
    const StaggerLocation loc = fields.descriptor(fid).location;
    const EntityAxis expected = loc == StaggerLocation::EdgeXi ? EntityAxis::Xi :
                                (loc == StaggerLocation::EdgeEt ? EntityAxis::Eta : EntityAxis::Zeta);
    std::vector<double> lsum(entries_.size(),0.0), gsum(entries_.size(),0.0);
    for (std::size_t n=0; n<entries_.size(); ++n)
        for (const auto &c : entries_[n].local_incident_cells)
            lsum[n] += c.weight*contribution(entries_[n],c);
    MPI_Allreduce(lsum.data(),gsum.data(),static_cast<int>(entries_.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (std::size_t n=0; n<entries_.size(); ++n)
    {
        const auto &e=entries_[n];
        if (e.owner.rank!=rank_ || e.owner.axis!=expected) continue;
        int owner_sign=+1;
        for (const auto &a:e.aliases) if (a.owner) { owner_sign=a.orientation; break; }
        fields.field(fid,e.owner.block)(e.owner.i,e.owner.j,e.owner.k,0)=owner_sign*gsum[n];
    }
}

void SingularEdgeRegistry::assemble_face_triplet_to_local_owners(
    Field &fields, const std::array<std::string,3> &field_names,
    const FaceContribution &contribution) const
{
    if (entries_.empty()) return;
    std::vector<double> lsum(entries_.size(),0.0), gsum(entries_.size(),0.0);
    for (std::size_t n=0; n<entries_.size(); ++n)
        for (const auto &f : entries_[n].local_incident_faces)
            lsum[n] += contribution(entries_[n],f);
    MPI_Allreduce(lsum.data(),gsum.data(),static_cast<int>(entries_.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (std::size_t n=0; n<entries_.size(); ++n)
    {
        const auto &e=entries_[n];
        if (e.owner.rank!=rank_) continue;
        int owner_sign=+1;
        for (const auto &a:e.aliases) if (a.owner) { owner_sign=a.orientation; break; }
        const int fid=fields.field_id(field_names[static_cast<int>(e.owner.axis)]);
        fields.field(fid,e.owner.block)(e.owner.i,e.owner.j,e.owner.k,0)=owner_sign*gsum[n];
    }
}

void SingularEdgeRegistry::assemble_cell_vector_circulation_to_local_owners(
    Field &fields,
    const std::string &cell_vector_field_name,
    const std::array<std::string,3> &edge_field_names,
    double regular_edge_blend) const
{
    if(curl_entries_.empty()) return;

    std::vector<std::size_t> offset(curl_entries_.size()+1,0);
    for(std::size_t n=0;n<curl_entries_.size();++n)
        offset[n+1]=offset[n]+curl_entries_[n].ordered_cell_centers.size();
    const std::size_t nsector=offset.back();
    std::vector<double> local_b(3*nsector,0.0), global_b(3*nsector,0.0);
    std::vector<int> local_count(nsector,0), global_count(nsector,0);
    std::vector<double> local_old(curl_entries_.size(),0.0);
    std::vector<double> global_old(curl_entries_.size(),0.0);

    const int cell_fid=fields.field_id(cell_vector_field_name);
    for(std::size_t n=0;n<curl_entries_.size();++n)
        for(const auto &inc:curl_entries_[n].local_incident_cells)
        {
            if(inc.sector_index<0) continue;
            FieldBlock &b=fields.field(cell_fid,inc.entity.block);
            if(!b.is_allocated()) continue;
            const std::size_t s=offset[n]+static_cast<std::size_t>(inc.sector_index);
            for(int q=0;q<3;++q)
                local_b[3*s+q]+=b(inc.entity.i,inc.entity.j,inc.entity.k,q);
            ++local_count[s];
        }

    for(std::size_t n=0;n<curl_entries_.size();++n)
    {
        const auto &e=curl_entries_[n];
        if(e.owner.rank!=rank_) continue;
        int owner_sign=+1;
        for(const auto &a:e.aliases) if(a.owner) { owner_sign=a.orientation; break; }
        const int fid=fields.field_id(edge_field_names[static_cast<int>(e.owner.axis)]);
        local_old[n]=owner_sign*
            fields.field(fid,e.owner.block)(e.owner.i,e.owner.j,e.owner.k,0);
    }

    MPI_Allreduce(local_b.data(),global_b.data(),static_cast<int>(global_b.size()),
                  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_count.data(),global_count.data(),static_cast<int>(global_count.size()),
                  MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_old.data(),global_old.data(),static_cast<int>(global_old.size()),
                  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for(std::size_t n=0;n<curl_entries_.size();++n)
    {
        const auto &e=curl_entries_[n];
        const std::size_t ns=e.ordered_cell_centers.size();
        if(ns<3) continue;
        double circulation=0.0;
        for(std::size_t m=0;m<ns;++m)
        {
            const std::size_t a=offset[n]+m;
            const std::size_t b=offset[n]+(m+1)%ns;
            if(global_count[a]!=1 || global_count[b]!=1)
                ERROR::Abort("SingularEdgeRegistry: a dual polygon sector does not map to one unique real cell");
            for(int q=0;q<3;++q)
            {
                const double bmid=0.5*(global_b[3*a+q]+global_b[3*b+q]);
                const double dx=e.ordered_cell_centers[(m+1)%ns][q]-
                                e.ordered_cell_centers[m][q];
                circulation+=bmid*dx;
            }
        }
        if(e.owner.rank!=rank_) continue;
        int owner_sign=+1;
        for(const auto &a:e.aliases) if(a.owner) { owner_sign=a.orientation; break; }
        const int fid=fields.field_id(edge_field_names[static_cast<int>(e.owner.axis)]);
        const double consistent=e.inverse_hodge*circulation;
        const double blend=e.variable_valence ? 1.0 : regular_edge_blend;
        const double canonical=global_old[n]+blend*(consistent-global_old[n]);
        fields.field(fid,e.owner.block)(e.owner.i,e.owner.j,e.owner.k,0)=
            owner_sign*canonical;
    }
}

void SingularEdgeRegistry::assemble_consistent_face_coboundary_to_local_owners(
    Field &fields,
    const std::array<std::string,3> &face_field_names,
    const std::array<std::string,3> &edge_field_names) const
{
    if(curl_entries_.empty()) return;

    // M2B is a quotient face 1-cochain.  Assemble d1^T M2B on the unique
    // physical edge first, and apply the physical lumped M1 inverse exactly
    // once.  Averaging already-scaled block-local J aliases is not equivalent
    // on a variable-valence edge.
    std::vector<double> local(curl_entries_.size(),0.0);
    std::vector<double> global(curl_entries_.size(),0.0);
    for(std::size_t n=0;n<curl_entries_.size();++n)
        for(const auto &inc:curl_entries_[n].local_incident_faces)
        {
            if(inc.quotient_incidence==0) continue;
            const int fid=fields.field_id(face_field_names[static_cast<int>(inc.entity.axis)]);
            FieldBlock &h=fields.field(fid,inc.entity.block);
            if(!h.is_allocated()) continue;
            const double canonical_face=inc.entity_orientation*
                h(inc.entity.i,inc.entity.j,inc.entity.k,0);
            local[n]+=inc.quotient_incidence*canonical_face;
        }

    MPI_Allreduce(local.data(),global.data(),static_cast<int>(global.size()),
                  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(std::size_t n=0;n<curl_entries_.size();++n)
    {
        const auto &e=curl_entries_[n];
        if(e.owner.rank!=rank_) continue;
        int owner_sign=+1;
        for(const auto &a:e.aliases) if(a.owner) { owner_sign=a.orientation; break; }
        const int fid=fields.field_id(edge_field_names[static_cast<int>(e.owner.axis)]);
        fields.field(fid,e.owner.block)(e.owner.i,e.owner.j,e.owner.k,0)=
            owner_sign*e.inverse_hodge*global[n];
    }
}

void SingularEdgeRegistry::assemble_cell_vector_affine_curl_to_local_owners(
    Field &fields,
    const std::string &cell_vector_field_name,
    const std::array<std::string,3> &edge_field_names) const
{
    if(entries_.empty()) return;

    // One unique Cartesian B sample per physical sector.  An affine fit in
    // the plane normal to the quotient edge supplies the polynomial-exact
    // target of M1_lumped^{-1} d1^T M2: J_e=(curl B . t_hat)|e|.
    std::vector<std::size_t> offset(entries_.size()+1,0);
    for(std::size_t n=0;n<entries_.size();++n)
        offset[n+1]=offset[n]+entries_[n].ordered_cell_centers.size();
    const std::size_t nsector=offset.back();
    std::vector<double> local_b(3*nsector,0.0),global_b(3*nsector,0.0);
    std::vector<int> local_count(nsector,0),global_count(nsector,0);
    const int cell_fid=fields.field_id(cell_vector_field_name);
    for(std::size_t n=0;n<entries_.size();++n)
        for(const auto &inc:entries_[n].local_incident_cells)
        {
            if(inc.sector_index<0) continue;
            FieldBlock &b=fields.field(cell_fid,inc.entity.block);
            if(!b.is_allocated()) continue;
            const std::size_t s=offset[n]+static_cast<std::size_t>(inc.sector_index);
            for(int q=0;q<3;++q)
                local_b[3*s+q]+=b(inc.entity.i,inc.entity.j,inc.entity.k,q);
            ++local_count[s];
        }
    MPI_Allreduce(local_b.data(),global_b.data(),static_cast<int>(global_b.size()),
                  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_count.data(),global_count.data(),static_cast<int>(global_count.size()),
                  MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    auto solve3=[](double A[3][3],double rhs[3],double x[3])->bool
    {
        double aug[3][4];
        double scale=0.0;
        for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        { aug[i][j]=A[i][j]; scale=std::max(scale,std::abs(A[i][j])); }
        for(int i=0;i<3;++i) aug[i][3]=rhs[i];
        const double tol=std::max(1.e-30,1.e-12*scale);
        for(int col=0;col<3;++col)
        {
            int pivot=col;
            for(int row=col+1;row<3;++row)
                if(std::abs(aug[row][col])>std::abs(aug[pivot][col])) pivot=row;
            if(std::abs(aug[pivot][col])<=tol) return false;
            if(pivot!=col) for(int j=col;j<4;++j) std::swap(aug[col][j],aug[pivot][j]);
            const double d=aug[col][col];
            for(int j=col;j<4;++j) aug[col][j]/=d;
            for(int row=0;row<3;++row) if(row!=col)
            {
                const double f=aug[row][col];
                for(int j=col;j<4;++j) aug[row][j]-=f*aug[col][j];
            }
        }
        for(int i=0;i<3;++i) x[i]=aug[i][3];
        return std::isfinite(x[0])&&std::isfinite(x[1])&&std::isfinite(x[2]);
    };

    for(std::size_t n=0;n<entries_.size();++n)
    {
        const auto &edge=entries_[n];
        const std::size_t ns=edge.ordered_cell_centers.size();
        if(ns<3) continue;
        for(std::size_t m=0;m<ns;++m)
            if(global_count[offset[n]+m]!=1)
                ERROR::Abort("Singular affine patch sector does not map to one unique real cell");

        const double L=std::sqrt(edge.canonical_edge_vector[0]*edge.canonical_edge_vector[0]+
                                 edge.canonical_edge_vector[1]*edge.canonical_edge_vector[1]+
                                 edge.canonical_edge_vector[2]*edge.canonical_edge_vector[2]);
        if(!(L>0.0)) continue;
        const std::array<double,3> t{{edge.canonical_edge_vector[0]/L,
                                      edge.canonical_edge_vector[1]/L,
                                      edge.canonical_edge_vector[2]/L}};
        const std::array<double,3> ref=std::abs(t[0])<0.9
            ? std::array<double,3>{{1.0,0.0,0.0}}
            : std::array<double,3>{{0.0,1.0,0.0}};
        std::array<double,3> u{{t[1]*ref[2]-t[2]*ref[1],
                                t[2]*ref[0]-t[0]*ref[2],
                                t[0]*ref[1]-t[1]*ref[0]}};
        const double un=std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        if(!(un>0.0)) continue;
        for(double &q:u) q/=un;
        const std::array<double,3> v{{t[1]*u[2]-t[2]*u[1],
                                      t[2]*u[0]-t[0]*u[2],
                                      t[0]*u[1]-t[1]*u[0]}};

        std::array<double,3> center{{0.0,0.0,0.0}};
        for(const auto &p:edge.ordered_cell_centers)
            for(int q=0;q<3;++q) center[q]+=p[q]/static_cast<double>(ns);
        std::vector<double> px(ns),py(ns),bu(ns),bv(ns);
        double radius2=0.0;
        for(std::size_t m=0;m<ns;++m)
        {
            std::array<double,3> dr{{edge.ordered_cell_centers[m][0]-center[0],
                                     edge.ordered_cell_centers[m][1]-center[1],
                                     edge.ordered_cell_centers[m][2]-center[2]}};
            px[m]=dr[0]*u[0]+dr[1]*u[1]+dr[2]*u[2];
            py[m]=dr[0]*v[0]+dr[1]*v[1]+dr[2]*v[2];
            radius2+=px[m]*px[m]+py[m]*py[m];
            const double *b=&global_b[3*(offset[n]+m)];
            bu[m]=b[0]*u[0]+b[1]*u[1]+b[2]*u[2];
            bv[m]=b[0]*v[0]+b[1]*v[1]+b[2]*v[2];
        }
        const double radius=std::sqrt(radius2/static_cast<double>(ns));
        if(!(radius>0.0) || !std::isfinite(radius)) continue;

        // Normalize transverse coordinates before forming the normal matrix.
        // This keeps the fixed small solve well-conditioned across radial
        // refinement levels; derivatives are divided by radius afterwards.
        double normal[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        double rhs_u[3]={0,0,0},rhs_v[3]={0,0,0};
        for(std::size_t m=0;m<ns;++m)
        {
            const double p[3]={1.0,px[m]/radius,py[m]/radius};
            for(int i=0;i<3;++i)
            {
                rhs_u[i]+=p[i]*bu[m]; rhs_v[i]+=p[i]*bv[m];
                for(int j=0;j<3;++j) normal[i][j]+=p[i]*p[j];
            }
        }
        double Au[3][3],Av[3][3];
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) Au[i][j]=Av[i][j]=normal[i][j];
        double cu[3],cv[3];
        if(!solve3(Au,rhs_u,cu) || !solve3(Av,rhs_v,cv)) continue;
        const double curl_t=(cv[1]-cu[2])/radius;
        const double canonical=L*curl_t;
        if(!std::isfinite(canonical)) continue;

        if(edge.owner.rank!=rank_) continue;
        int owner_sign=+1;
        for(const auto &a:edge.aliases) if(a.owner){owner_sign=a.orientation;break;}
        const int fid=fields.field_id(edge_field_names[static_cast<int>(edge.owner.axis)]);
        fields.field(fid,edge.owner.block)(edge.owner.i,edge.owner.j,edge.owner.k,0)=
            owner_sign*canonical;
    }
}
} // namespace METRIC
