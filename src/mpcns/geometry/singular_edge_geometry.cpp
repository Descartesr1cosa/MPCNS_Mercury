#include "mpcns/geometry/singular_edge_geometry.hpp"

#include <algorithm>
#include <cmath>

namespace mpcns::geometry
{
namespace
{

double norm(const CartesianVector &vector)
{
    return std::sqrt(vector[0] * vector[0] +
                     vector[1] * vector[1] +
                     vector[2] * vector[2]);
}

} // namespace

SingularEdgeGeometry close_singular_edge_geometry(
    const SingularEdgeGeometryInput &input)
{
    SingularEdgeGeometry result;

    if (input.tangent_sample_count > 0)
        for (int component = 0; component < 3; ++component)
            result.canonical_edge_vector[component] =
                input.tangent_sum[component] / input.tangent_sample_count;

    if (input.incident_dual_centers.size() >= 3)
    {
        std::vector<CartesianVector> points = input.incident_dual_centers;
        CartesianVector center{{0.0, 0.0, 0.0}};
        for (const CartesianVector &point : points)
            for (int component = 0; component < 3; ++component)
                center[component] += point[component];
        for (int component = 0; component < 3; ++component)
            center[component] /= points.size();

        CartesianVector tangent = input.tangent_sum;
        const double tangent_norm = norm(tangent);
        if (tangent_norm > 0.0)
            for (double &component : tangent)
                component /= tangent_norm;

        const CartesianVector reference =
            std::abs(tangent[0]) < 0.9
                ? CartesianVector{{1.0, 0.0, 0.0}}
                : CartesianVector{{0.0, 1.0, 0.0}};
        CartesianVector u{{tangent[1] * reference[2] -
                               tangent[2] * reference[1],
                           tangent[2] * reference[0] -
                               tangent[0] * reference[2],
                           tangent[0] * reference[1] -
                               tangent[1] * reference[0]}};
        const double u_norm = norm(u);
        if (u_norm > 0.0)
            for (double &component : u)
                component /= u_norm;
        const CartesianVector v{{tangent[1] * u[2] - tangent[2] * u[1],
                                 tangent[2] * u[0] - tangent[0] * u[2],
                                 tangent[0] * u[1] - tangent[1] * u[0]}};

        std::sort(points.begin(), points.end(),
                  [&](const CartesianVector &lhs, const CartesianVector &rhs)
                  {
                      const auto angle = [&](const CartesianVector &point)
                      {
                          double du = 0.0;
                          double dv = 0.0;
                          for (int component = 0; component < 3; ++component)
                          {
                              du += (point[component] - center[component]) *
                                    u[component];
                              dv += (point[component] - center[component]) *
                                    v[component];
                          }
                          return std::atan2(dv, du);
                      };
                      return angle(lhs) < angle(rhs);
                  });

        CartesianVector area_vector{{0.0, 0.0, 0.0}};
        for (std::size_t index = 0; index < points.size(); ++index)
        {
            const CartesianVector &lhs = points[index];
            const CartesianVector &rhs = points[(index + 1) % points.size()];
            area_vector[0] += lhs[1] * rhs[2] - lhs[2] * rhs[1];
            area_vector[1] += lhs[2] * rhs[0] - lhs[0] * rhs[2];
            area_vector[2] += lhs[0] * rhs[1] - lhs[1] * rhs[0];
        }
        result.dual_area = 0.5 * norm(area_vector);
        result.used_incident_polygon = result.dual_area > 0.0;
    }

    if (!(result.dual_area > 0.0))
        result.dual_area = input.fallback_dual_area;
    result.inverse_hodge =
        result.dual_area > 0.0 ? input.primal_length / result.dual_area : 0.0;
    return result;
}

} // namespace mpcns::geometry
