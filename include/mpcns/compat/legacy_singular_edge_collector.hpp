#pragma once

#include "2_topology/Topology.h"
#include "mpcns/geometry/geometry_primitives.hpp"
#include "mpcns/mesh/singular_edge_records.hpp"

#include <array>
#include <vector>

class Field;
class Grid;

namespace mpcns::compat
{

// Explicit legacy data source for the singular-edge compatibility boundary.
// It keeps Grid/Field access out of SingularEdgeRegistry's public contract.
struct LegacySingularEdgeDataSource
{
    Field &fields;
    Grid &grid;
};

// Local, pre-MPI data for one quotient-edge candidate. The registry owns
// global reduction, physical-valence classification, geometry closure, and
// normalized weights; this collector owns only legacy Field/Grid reads.
struct LegacySingularEdgeLocalMeasures
{
    double primal_length_sum = 0.0;
    std::array<double, 3> tangent_sum{{0.0, 0.0, 0.0}};
    double cell_measure_sum = 0.0;
    double face_measure_sum = 0.0;
    int primal_length_count = 0;
    int valid_cell_count = 0;
    int valid_face_count = 0;

    std::vector<mesh::WeightedIncidentEntity> incident_cells;
    std::vector<mesh::WeightedIncidentEntity> incident_faces;
    std::vector<geometry::CartesianVector> incident_cell_centers;
};

LegacySingularEdgeLocalMeasures collect_legacy_singular_edge_measures(
    const TOPO::Topology &topology,
    const LegacySingularEdgeDataSource &source,
    int local_rank,
    const std::vector<mesh::SingularEdgeAlias> &aliases);

} // namespace mpcns::compat
