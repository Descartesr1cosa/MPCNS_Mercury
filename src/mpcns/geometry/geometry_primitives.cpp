#include "mpcns/geometry/geometry_primitives.hpp"

namespace mpcns::geometry
{
namespace
{

CartesianVector subtract(const CartesianVector &lhs,
                         const CartesianVector &rhs)
{
    return {lhs[0] - rhs[0],
            lhs[1] - rhs[1],
            lhs[2] - rhs[2]};
}

CartesianVector cross(const CartesianVector &lhs,
                      const CartesianVector &rhs)
{
    return {lhs[1] * rhs[2] - lhs[2] * rhs[1],
            lhs[2] * rhs[0] - lhs[0] * rhs[2],
            lhs[0] * rhs[1] - lhs[1] * rhs[0]};
}

double dot(const CartesianVector &lhs,
           const CartesianVector &rhs)
{
    return lhs[0] * rhs[0] +
           lhs[1] * rhs[1] +
           lhs[2] * rhs[2];
}

} // namespace

CartesianVector quadrilateral_area_vector(const CartesianVector &p00,
                                          const CartesianVector &p10,
                                          const CartesianVector &p01,
                                          const CartesianVector &p11)
{
    const CartesianVector first_triangle =
        cross(subtract(p10, p00), subtract(p01, p00));
    const CartesianVector second_triangle =
        cross(subtract(p01, p11), subtract(p10, p11));
    return {0.5 * (first_triangle[0] + second_triangle[0]),
            0.5 * (first_triangle[1] + second_triangle[1]),
            0.5 * (first_triangle[2] + second_triangle[2])};
}

double cell_volume_from_oriented_faces(
    const CartesianVector &reference_point,
    const std::array<OrientedFaceGeometry, 6> &faces)
{
    double signed_volume = 0.0;
    for (const OrientedFaceGeometry &face : faces)
    {
        signed_volume += dot(face.area_vector,
                             subtract(face.centroid, reference_point));
    }
    return signed_volume / 3.0;
}

} // namespace mpcns::geometry
