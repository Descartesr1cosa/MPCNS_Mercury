#pragma once

#include <array>

namespace mpcns::discrete
{

// Ordered corners of the dual quadrilateral surrounding one primal edge.
// The order p00, p10, p01, p11 matches the legacy structured-grid stencil.
struct DualFaceQuadrilateral
{
    std::array<double, 3> p00;
    std::array<double, 3> p10;
    std::array<double, 3> p01;
    std::array<double, 3> p11;
};

// Safety thresholds for the diagonal inverse star-1 coefficient near a
// collapsed dual face. Defaults preserve the established Mercury treatment.
struct AxisSafeInverseStar1Parameters
{
    double dual_collapse_ratio = 1.0e-3;
    double primal_edge_absolute_tolerance = 1.0e-12;
};

// Computes |e| / |*e| unless the dual quadrilateral is axis-touching or
// degenerate. This is a pure, Field-independent policy so its thresholds and
// geometric conditions can be reviewed separately from metric storage.
class AxisSafeInverseStar1Policy
{
public:
    explicit AxisSafeInverseStar1Policy(
        AxisSafeInverseStar1Parameters parameters = {});

    double apply(double primal_edge_length,
                 const DualFaceQuadrilateral &dual_face,
                 double dual_face_area,
                 double zero_area_tolerance) const;

private:
    AxisSafeInverseStar1Parameters parameters_;
};

} // namespace mpcns::discrete
