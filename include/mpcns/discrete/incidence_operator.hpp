#pragma once

#include "mpcns/state/cochain_view.hpp"

namespace mpcns::discrete
{

namespace detail
{

inline double d1_xi_at(const state::PrimalEdgeFormView &edge, int i, int j, int k)
{
    return (edge.eta(i, j, k, edge.component) -
            edge.eta(i, j, k + 1, edge.component)) +
           (edge.zeta(i, j + 1, k, edge.component) -
            edge.zeta(i, j, k, edge.component));
}

inline double d1_eta_at(const state::PrimalEdgeFormView &edge, int i, int j, int k)
{
    return (edge.xi(i, j, k + 1, edge.component) -
            edge.xi(i, j, k, edge.component)) +
           (edge.zeta(i, j, k, edge.component) -
            edge.zeta(i + 1, j, k, edge.component));
}

inline double d1_zeta_at(const state::PrimalEdgeFormView &edge, int i, int j, int k)
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
void apply_d1(const state::PrimalEdgeFormView &edge,
              state::PrimalFaceFormView face,
              double scale = 1.0);

// Subtract a scaled oriented primal incidence operator from an existing
// 2-form.  This preserves the established Faraday update order,
// face -= dt * d1(edge).
void subtract_d1(const state::PrimalEdgeFormView &edge,
                 state::PrimalFaceFormView face,
                 double scale);

// Apply the transpose incidence operator d1^T: dual 1-form -> dual 2-form.
// The caller supplies the Hodge-transformed dual 1-form; no geometric scale is
// embedded in this purely topological operator.
void apply_d1_transpose(const state::DualEdgeFormView &dual_edge,
                        state::DualFaceFormView dual_face,
                        double scale = 1.0);

} // namespace mpcns::discrete
