#include "7_metric/SingularEdgeRegistry.h"

#include "0_basic/Error.h"
#include "0_basic/MPI_WRAPPER.h"
#include "3_field/Field.h"
#include "1_grid/1_MPCNS_Grid.h"

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
    entries_.clear(); gid_to_index_.clear();
    const int nglobal = topology.edges.n_global_owner;
    std::vector<double> local_length(nglobal,0.0), local_cell_measure(nglobal,0.0);
    std::vector<double> local_face_measure(nglobal,0.0);
    std::vector<double> local_dr(3*nglobal,0.0);
    std::vector<double> local_center_records;
    std::vector<int> local_length_n(nglobal,0), local_cells(nglobal,0), local_faces(nglobal,0);

    for (const auto &cls : topology.edges.classes)
    {
        // members counts block-local aliases, not physical sectors.  A
        // conforming block decomposition can give a regular edge three
        // aliases.  Treat this as a candidate only; physical valence is
        // classified below from quotient-face incidence.
        if (cls.members.size() < 3 || cls.global_id < 0)
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
                    { rec.local_incident_faces.push_back({f,av,0.0,e,a.orientation}); local_face_measure[rec.global_id]+=av; ++local_faces[rec.global_id]; }
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

    // Singularity is the valence of quotient faces incident on the physical
    // edge, not the number of block-local edge aliases (or cells).  Three and
    // greater-than-four face sectors use the variable-valence path.  Four is
    // a regular interior edge; two is a regular physical-boundary edge.
    entries_.erase(std::remove_if(entries_.begin(), entries_.end(),
                                  [&](const SingularPhysicalEdge &e)
                                  {
                                      const int physical_faces = gf[e.global_id];
                                      return physical_faces != 3 && physical_faces <= 4;
                                  }),
                   entries_.end());
    gid_to_index_.clear();

    for (std::size_t entry_index = 0; entry_index < entries_.size(); ++entry_index)
    {
        auto &e = entries_[entry_index];
        gid_to_index_[e.global_id] = entry_index;
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

        // Publish the corrected scalar Hodge data to every local alias.  The
        // vector-valued Sstar remains untouched until a dual-polygon vector
        // reconstruction is requested by an operator.
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
} // namespace METRIC
