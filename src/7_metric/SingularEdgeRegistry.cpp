#include "7_metric/SingularEdgeRegistry.h"

#include "0_basic/Error.h"
#include "0_basic/MPI_WRAPPER.h"
#include "3_field/Field.h"

#include <algorithm>
#include <cmath>
#include <iostream>
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
void SingularEdgeRegistry::build(const TOPO::Topology &topology, Field &fields, int rank)
{
    rank_ = rank;
    entries_.clear(); gid_to_index_.clear();
    const int nglobal = topology.edges.n_global_owner;
    std::vector<double> local_length(nglobal,0.0), local_dual(nglobal,0.0), local_cell_measure(nglobal,0.0);
    std::vector<int> local_length_n(nglobal,0), local_cells(nglobal,0), local_faces(nglobal,0);

    for (const auto &cls : topology.edges.classes)
    {
        // A regular multi-block interior edge has four block-local sectors.
        // The registry is intentionally valence based, so future 5+ valence
        // edges enter the same path without another special case.
        if (cls.members.size() == 4 || cls.members.size() < 3 || cls.global_id < 0)
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
                FieldBlock &as = fields.field(std::string("Astar_") + axis_suffix(e.axis), e.block);
                if (std::isfinite(L) && L > 0.0) { local_length[rec.global_id] += L; ++local_length_n[rec.global_id]; }
                const double A = as(e.i,e.j,e.k,0);
                if (std::isfinite(A) && A > 0.0) local_dual[rec.global_id] += A;
            }

            FieldBlock &jac = fields.field("Jac", e.block);
            for (const EntityKey &c : incident_cells(e))
                if (seen_cells.insert(c).second && contains_inner(jac,c.i,c.j,c.k))
                {
                    const double v = std::abs(jac(c.i,c.j,c.k,0));
                    if (std::isfinite(v) && v > 1.e-30)
                    {
                        rec.local_incident_cells.push_back({c,v,0.0});
                        local_cell_measure[rec.global_id] += v; ++local_cells[rec.global_id];
                    }
                }

            for (const EntityKey &f : incident_faces(e))
                if (seen_faces.insert(f).second)
                {
                    FieldBlock &area = fields.field(std::string("Area_") + axis_suffix(f.axis), f.block);
                    if (!contains_inner(area,f.i,f.j,f.k)) continue;
                    const double av = area(f.i,f.j,f.k,0);
                    if (std::isfinite(av) && av > 1.e-30)
                    { rec.local_incident_faces.push_back({f,av,0.0}); ++local_faces[rec.global_id]; }
                }
        }
        gid_to_index_[rec.global_id] = entries_.size(); entries_.push_back(std::move(rec));
    }

    std::vector<double> gl(nglobal), ga(nglobal), gv(nglobal);
    std::vector<int> gln(nglobal), gc(nglobal), gf(nglobal);
    MPI_Allreduce(local_length.data(),gl.data(),nglobal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_dual.data(),ga.data(),nglobal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_cell_measure.data(),gv.data(),nglobal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_length_n.data(),gln.data(),nglobal,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_cells.data(),gc.data(),nglobal,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_faces.data(),gf.data(),nglobal,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    for (auto &e : entries_)
    {
        e.primal_length = gln[e.global_id] ? gl[e.global_id]/gln[e.global_id] : 0.0;
        // Each alias constructs the same complete dual polygon; average it.
        e.dual_area = gln[e.global_id] ? ga[e.global_id]/gln[e.global_id] : 0.0;
        e.inverse_hodge = e.dual_area > 0.0 ? e.primal_length/e.dual_area : 0.0;
        e.global_valid_cell_count=gc[e.global_id]; e.global_valid_face_count=gf[e.global_id];
        for (auto &c : e.local_incident_cells) c.weight = gv[e.global_id] > 0.0 ? c.measure/gv[e.global_id] : 0.0;
        double face_sum=0.0; for (const auto &f:e.local_incident_faces) face_sum += f.measure;
        for (auto &f:e.local_incident_faces) f.weight = face_sum > 0.0 ? f.measure/face_sum : 0.0;
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

void SingularEdgeRegistry::reduce_field_to_local_owners(Field &fields,
                                                         const std::string &field_name) const
{
    if (entries_.empty()) return;
    const int fid = fields.field_id(field_name);
    const StaggerLocation loc = fields.descriptor(fid).location;
    const EntityAxis expected = loc == StaggerLocation::EdgeXi ? EntityAxis::Xi :
                                (loc == StaggerLocation::EdgeEt ? EntityAxis::Eta : EntityAxis::Zeta);
    std::vector<double> lsum(entries_.size(),0.0), gsum(entries_.size(),0.0);
    std::vector<double> ln(entries_.size(),0.0), gn(entries_.size(),0.0);
    for (std::size_t n=0; n<entries_.size(); ++n)
        for (const auto &a : entries_[n].aliases)
        {
            const EntityKey &e=a.edge;
            if (e.rank!=rank_ || e.axis!=expected || e.block<0 || e.block>=fields.num_blocks()) continue;
            FieldBlock &f=fields.field(fid,e.block);
            if (!f.is_allocated() || !contains_inner(f,e.i,e.j,e.k)) continue;
            const double v=f(e.i,e.j,e.k,0);
            if (!std::isfinite(v)) continue;
            lsum[n] += a.orientation*v; ln[n] += 1.0;
        }
    MPI_Allreduce(lsum.data(),gsum.data(),static_cast<int>(entries_.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(ln.data(),gn.data(),static_cast<int>(entries_.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for (std::size_t n=0; n<entries_.size(); ++n)
    {
        const auto &e=entries_[n];
        if (e.owner.rank!=rank_ || e.owner.axis!=expected || gn[n]<=0.0) continue;
        int owner_sign=+1;
        for (const auto &a:e.aliases) if (a.owner) { owner_sign=a.orientation; break; }
        fields.field(fid,e.owner.block)(e.owner.i,e.owner.j,e.owner.k,0)=owner_sign*gsum[n]/gn[n];
    }
}
} // namespace METRIC
