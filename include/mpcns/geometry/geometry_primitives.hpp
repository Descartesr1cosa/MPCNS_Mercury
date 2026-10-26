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

} // namespace mpcns::geometry
