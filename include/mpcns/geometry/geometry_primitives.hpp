#pragma once

#include <array>

namespace mpcns::geometry
{

using CartesianVector = std::array<double, 3>;

// Oriented area vector of a possibly non-planar quadrilateral.  The points
// are ordered as p00, p10, p01, p11, matching the two-triangle stencil used
// by the legacy structured-grid metrics.
CartesianVector quadrilateral_area_vector(const CartesianVector &p00,
                                          const CartesianVector &p10,
                                          const CartesianVector &p01,
                                          const CartesianVector &p11);

// Euclidean norm of a Cartesian geometric vector.
double euclidean_norm(const CartesianVector &vector);

// Area vector and centroid of one cell face, oriented outward from the cell.
struct OrientedFaceGeometry
{
    CartesianVector area_vector;
    CartesianVector centroid;
};

// Volume from the divergence-theorem face quadrature used by the legacy
// structured-grid metric.  Faces must be ordered deterministically by the
// caller so the floating-point accumulation remains reproducible.
double cell_volume_from_oriented_faces(
    const CartesianVector &reference_point,
    const std::array<OrientedFaceGeometry, 6> &faces);

} // namespace mpcns::geometry
