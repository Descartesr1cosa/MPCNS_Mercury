#include "7_metric/SingularEdgeRegistry.h"

#include "0_basic/Error.h"
#include "0_basic/MPI_WRAPPER.h"
#include "3_field/Field.h"
#include "mpcns/compat/legacy_singular_edge_collector.hpp"
#include "mpcns/compat/legacy_singular_edge_metric_publisher.hpp"
#include "mpcns/geometry/singular_edge_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <sstream>

namespace
{
using TOPO::EntityAxis;
} // namespace

namespace METRIC
{
void SingularEdgeRegistry::build(
    const TOPO::Topology &topology,
    const mpcns::compat::LegacySingularEdgeDataSource &source,
    int rank)
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

        const mpcns::compat::LegacySingularEdgeLocalMeasures local =
            mpcns::compat::collect_legacy_singular_edge_measures(
                topology, source, rank, rec.aliases);
        rec.local_incident_cells = local.incident_cells;
        rec.local_incident_faces = local.incident_faces;

        local_length[rec.global_id] += local.primal_length_sum;
        local_cell_measure[rec.global_id] += local.cell_measure_sum;
        local_face_measure[rec.global_id] += local.face_measure_sum;
        for (int component = 0; component < 3; ++component)
            local_dr[3 * rec.global_id + component] +=
                local.tangent_sum[component];
        local_length_n[rec.global_id] += local.primal_length_count;
        local_cells[rec.global_id] += local.valid_cell_count;
        local_faces[rec.global_id] += local.valid_face_count;
        for (const mpcns::geometry::CartesianVector &center :
             local.incident_cell_centers)
        {
            local_center_records.push_back(
                static_cast<double>(rec.global_id));
            local_center_records.push_back(center[0]);
            local_center_records.push_back(center[1]);
            local_center_records.push_back(center[2]);
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
        e.primal_length =
            gln[e.global_id] ? gl[e.global_id] / gln[e.global_id] : 0.0;

        // Topology and MPI collect the physical incident centers above. The
        // quotient-edge dual polygon itself is now a pure geometry closure.
        mpcns::geometry::SingularEdgeGeometryInput closure_input;
        closure_input.primal_length = e.primal_length;
        closure_input.tangent_sum = {{
            gdr[3 * e.global_id],
            gdr[3 * e.global_id + 1],
            gdr[3 * e.global_id + 2]}};
        closure_input.tangent_sample_count = gln[e.global_id];
        const auto center_it = incident_centers.find(e.global_id);
        if (center_it != incident_centers.end())
            closure_input.incident_dual_centers = center_it->second;
        closure_input.fallback_dual_area =
            e.primal_length > 0.0
                ? gv[e.global_id] / e.primal_length
                : 0.0;

        const mpcns::geometry::SingularEdgeGeometry closure =
            mpcns::geometry::close_singular_edge_geometry(closure_input);
        e.dual_area = closure.dual_area;
        e.inverse_hodge = closure.inverse_hodge;
        e.canonical_edge_vector = closure.canonical_edge_vector;
        e.global_valid_cell_count=gc[e.global_id]; e.global_valid_face_count=gf[e.global_id];
        for (auto &c : e.local_incident_cells) c.weight = gv[e.global_id] > 0.0 ? c.measure/gv[e.global_id] : 0.0;
        for (auto &f:e.local_incident_faces) f.weight = gfm[e.global_id] > 0.0 ? f.measure/gfm[e.global_id] : 0.0;

        // Keep the legacy Field catalog at an explicit compatibility seam.
        mpcns::compat::publish_legacy_singular_edge_hodge(source.fields, e, rank_);
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
} // namespace METRIC
