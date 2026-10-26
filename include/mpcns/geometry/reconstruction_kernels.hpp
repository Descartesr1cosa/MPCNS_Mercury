#pragma once

#include "mpcns/geometry/geometry_primitives.hpp"

#include <array>

namespace mpcns::geometry
{

using Matrix3 = std::array<double, 9>;

struct CellCovariantBasis
{
    CartesianVector xi;
    CartesianVector eta;
    CartesianVector zeta;
};

// The corners are ordered r000, r100, r010, r110, r001, r101, r011, r111.
// This preserves the legacy four-edge averages and floating-point grouping.
CellCovariantBasis cell_covariant_basis_from_corners(
    const std::array<CartesianVector, 8> &corners);

// Copy an existing row-major C-array into the stable matrix value type.
Matrix3 matrix3_from_row_major(const double values[9]);

// Tikhonov pseudo-inverse used by the legacy reconstruction cache:
// (X X^T + reg I)^-1 X, with the same scale-aware regularization,
// eight retry attempts, determinant cutoff, and zero fallback as production.
Matrix3 regularized_pseudoinverse_transpose(const Matrix3 &matrix);

} // namespace mpcns::geometry
