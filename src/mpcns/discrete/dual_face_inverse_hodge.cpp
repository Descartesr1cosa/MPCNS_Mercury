#include "mpcns/discrete/dual_face_inverse_hodge.hpp"

#include <algorithm>
#include <cmath>

namespace mpcns::discrete
{
namespace
{

using CartesianPoint = std::array<double, 3>;

CartesianPoint subtract(const CartesianPoint &lhs,
                        const CartesianPoint &rhs)
{
    return {lhs[0] - rhs[0],
            lhs[1] - rhs[1],
            lhs[2] - rhs[2]};
}

double euclidean_norm(const CartesianPoint &vector)
{
    return std::sqrt(vector[0] * vector[0] +
                     vector[1] * vector[1] +
                     vector[2] * vector[2]);
}

} // namespace

AxisSafeInverseStar1Policy::AxisSafeInverseStar1Policy(
    AxisSafeInverseStar1Parameters parameters)
    : parameters_(parameters)
{
}

double AxisSafeInverseStar1Policy::apply(
    double primal_edge_length,
    const DualFaceQuadrilateral &dual_face,
    double dual_face_area,
    double zero_area_tolerance) const
{
    const double inverse_star_raw =
        (dual_face_area < zero_area_tolerance)
            ? 0.0
            : (primal_edge_length / dual_face_area);

    const double d00_10 = euclidean_norm(subtract(dual_face.p10, dual_face.p00));
    const double d00_01 = euclidean_norm(subtract(dual_face.p01, dual_face.p00));
    const double d10_11 = euclidean_norm(subtract(dual_face.p11, dual_face.p10));
    const double d01_11 = euclidean_norm(subtract(dual_face.p11, dual_face.p01));

    const double d00_11 = euclidean_norm(subtract(dual_face.p11, dual_face.p00));
    const double d10_01 = euclidean_norm(subtract(dual_face.p01, dual_face.p10));

    const double href = std::max(
        std::max(std::max(d00_10, d00_01), std::max(d10_11, d01_11)),
        std::max(std::max(d00_11, d10_01),
                 std::max(primal_edge_length, zero_area_tolerance)));

    const double side_min = std::min(
        std::min(d00_10, d00_01),
        std::min(d10_11, d01_11));
    const double diagonal_min = std::min(d00_11, d10_01);

    const bool dual_touches_axis =
        (primal_edge_length < parameters_.primal_edge_absolute_tolerance) ||
        (side_min < parameters_.dual_collapse_ratio * href) ||
        (diagonal_min < parameters_.dual_collapse_ratio * href) ||
        (dual_face_area < parameters_.dual_collapse_ratio * href * href);
    if (dual_touches_axis)
        return 0.0;

    return inverse_star_raw;
}

} // namespace mpcns::discrete
