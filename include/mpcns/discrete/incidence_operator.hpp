#pragma once

#include "3_field/FieldBlock.h"

namespace mpcns::discrete
{

// Transitional typed views over the existing three directional FieldBlock
// storage.  They preserve the current layout while preventing generic DEC
// kernels from depending on application field names or IdTriplet.
struct PrimalEdgeFormView
{
    FieldBlock &xi;
    FieldBlock &eta;
    FieldBlock &zeta;
    int component = 0;
};

struct PrimalFaceFormView
{
    FieldBlock &xi;
    FieldBlock &eta;
    FieldBlock &zeta;
    int component = 0;
};

namespace detail
{

inline double d1_xi_at(const PrimalEdgeFormView &edge, int i, int j, int k)
{
    return (edge.eta(i, j, k, edge.component) -
            edge.eta(i, j, k + 1, edge.component)) +
           (edge.zeta(i, j + 1, k, edge.component) -
            edge.zeta(i, j, k, edge.component));
}

inline double d1_eta_at(const PrimalEdgeFormView &edge, int i, int j, int k)
{
    return (edge.xi(i, j, k + 1, edge.component) -
            edge.xi(i, j, k, edge.component)) +
           (edge.zeta(i, j, k, edge.component) -
            edge.zeta(i + 1, j, k, edge.component));
}

inline double d1_zeta_at(const PrimalEdgeFormView &edge, int i, int j, int k)
{
    return (edge.xi(i, j, k, edge.component) -
            edge.xi(i, j + 1, k, edge.component)) +
           (edge.eta(i + 1, j, k, edge.component) -
            edge.eta(i, j, k, edge.component));
}

} // namespace detail

// Apply the oriented primal incidence operator d1: C^1 -> C^2.  The optional
// scale preserves the legacy CT multiplier without introducing a second curl
// stencil.
void apply_d1(const PrimalEdgeFormView &edge,
              PrimalFaceFormView face,
              double scale = 1.0);

// Subtract a scaled oriented primal incidence operator from an existing
// 2-form.  This preserves the established Faraday update order,
// face -= dt * d1(edge).
void subtract_d1(const PrimalEdgeFormView &edge,
                 PrimalFaceFormView face,
                 double scale);

} // namespace mpcns::discrete
