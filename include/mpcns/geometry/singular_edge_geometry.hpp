#pragma once

#include "mpcns/geometry/geometry_primitives.hpp"

#include <vector>

namespace mpcns::geometry
{

// Geometry-only input for the closure around a quotient/extraordinary edge.
// The registry is responsible for topology, ownership, MPI reduction, and
// collecting these values; this type has no dependency on any of them.
struct SingularEdgeGeometryInput
{
    double primal_length = 0.0;
    CartesianVector tangent_sum{{0.0, 0.0, 0.0}};
    int tangent_sample_count = 0;
    std::vector<CartesianVector> incident_dual_centers;
    double fallback_dual_area = 0.0;
};

struct SingularEdgeGeometry
{
    double dual_area = 0.0;
    double inverse_hodge = 0.0;
    CartesianVector canonical_edge_vector{{0.0, 0.0, 0.0}};
    bool used_incident_polygon = false;
};

// Close the dual polygon around a singular edge. If the polygon is unusable,
// retain the established volume-over-length fallback. This is a pure geometry
// operation and preserves the legacy ordering of the polygon arithmetic.
SingularEdgeGeometry close_singular_edge_geometry(
    const SingularEdgeGeometryInput &input);

} // namespace mpcns::geometry
