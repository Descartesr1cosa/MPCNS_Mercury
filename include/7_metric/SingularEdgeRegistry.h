#pragma once

#include "2_topology/Topology.h"

#include <cstddef>
#include <unordered_map>
#include <vector>
#include <string>
#include <array>
#include <functional>

class Field;
class Grid;

namespace METRIC
{
struct SingularEdgeAlias
{
    TOPO::EntityKey edge;
    int orientation = +1;
    bool owner = false;
};

struct WeightedIncidentEntity
{
    TOPO::EntityKey entity;
    double measure = 0.0;
    double weight = 0.0;
    TOPO::EntityKey source_alias{TOPO::EntityDim::Edge,0,0,0,0,0,TOPO::EntityAxis::Xi};
    int source_orientation = +1;
    int sector_index = -1;
    // Orientation of this entity to its quotient owner and the coefficient
    // in the canonical quotient coboundary row.  The latter is meaningful
    // for incident faces of an edge.
    int entity_orientation = +1;
    int quotient_incidence = 0;
};

// One record represents one edge of the quotient (physical) mesh.  Records
// are variable-valence: no four-cell/four-block assumption is made here.
struct SingularPhysicalEdge
{
    int global_id = -1;
    TOPO::EntityKey owner;
    std::vector<SingularEdgeAlias> aliases;
    std::vector<WeightedIncidentEntity> local_incident_cells;
    std::vector<WeightedIncidentEntity> local_incident_faces;

    double primal_length = 0.0;
    double dual_area = 0.0;
    double inverse_hodge = 0.0; // |edge| / |dual face|
    std::array<double,3> canonical_edge_vector{{0.0,0.0,0.0}};
    // Unique real-cell centers ordered counter-clockwise when viewed along
    // canonical_edge_vector.  They form the boundary of the dual face.
    std::vector<std::array<double,3>> ordered_cell_centers;
    bool variable_valence = false;
    int global_valid_cell_count = 0;
    int global_valid_face_count = 0;
};

class SingularEdgeRegistry
{
public:
    void build(const TOPO::Topology &topology, Field &fields, Grid &grid, int rank);
    void validate_or_abort() const;

    const std::vector<SingularPhysicalEdge> &entries() const { return entries_; }
    const SingularPhysicalEdge *find(int global_edge_id) const;
    using CellContribution = std::function<double(const SingularPhysicalEdge &,
                                                   const WeightedIncidentEntity &)>;
    using FaceContribution = std::function<double(const SingularPhysicalEdge &,
                                                   const WeightedIncidentEntity &)>;
    void assemble_cell_field_to_local_owners(Field &fields,
                                             const std::string &field_name,
                                             const CellContribution &contribution) const;
    void assemble_face_triplet_to_local_owners(Field &fields,
                                               const std::array<std::string,3> &field_names,
                                               const FaceContribution &contribution) const;
    void assemble_cell_vector_circulation_to_local_owners(
        Field &fields,
        const std::string &cell_vector_field_name,
        const std::array<std::string,3> &edge_field_names,
        double regular_edge_blend = 1.0) const;
    void assemble_consistent_face_coboundary_to_local_owners(
        Field &fields,
        const std::array<std::string,3> &face_field_names,
        const std::array<std::string,3> &edge_field_names) const;
    bool empty() const { return entries_.empty(); }
    std::size_t size() const { return entries_.size(); }

private:
    std::vector<SingularPhysicalEdge> entries_;
    // All shared quotient edges with a closed real-cell dual polygon.  This
    // includes ordinary two-panel seam edges as well as variable-valence
    // singular edges, but is used only by the current/curl operator.
    std::vector<SingularPhysicalEdge> curl_entries_;
    std::unordered_map<int, std::size_t> gid_to_index_;
    int rank_ = 0;
};
} // namespace METRIC
