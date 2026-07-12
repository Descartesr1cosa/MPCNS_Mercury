#pragma once

#include "2_topology/Topology.h"

#include <cstddef>
#include <unordered_map>
#include <vector>
#include <string>
#include <array>
#include <functional>

class Field;

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
    int global_valid_cell_count = 0;
    int global_valid_face_count = 0;
};

class SingularEdgeRegistry
{
public:
    void build(const TOPO::Topology &topology, Field &fields, int rank);
    void validate_or_abort() const;

    const std::vector<SingularPhysicalEdge> &entries() const { return entries_; }
    const SingularPhysicalEdge *find(int global_edge_id) const;
    using CellContribution = std::function<double(const SingularPhysicalEdge &,
                                                   const WeightedIncidentEntity &)>;
    void assemble_cell_field_to_local_owners(Field &fields,
                                             const std::string &field_name,
                                             const CellContribution &contribution) const;
    bool empty() const { return entries_.empty(); }
    std::size_t size() const { return entries_.size(); }

private:
    std::vector<SingularPhysicalEdge> entries_;
    std::unordered_map<int, std::size_t> gid_to_index_;
    int rank_ = 0;
};
} // namespace METRIC
