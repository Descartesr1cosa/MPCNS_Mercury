#pragma once

#include "mpcns/state/cochain_view.hpp"
#include "mpcns/compat/legacy_cochain_view.hpp"

namespace mpcns::discrete
{

// Oriented primal incidence operators. They contain no metric, Field catalog,
// MPI, ownership, or boundary decision; callers synchronize input cochains
// before application and reconcile output cochains through their Halo policy.
void apply_d0(const state::PrimalZeroFormView &node,
              state::PrimalEdgeFormView edge);

void apply_d1(const state::PrimalEdgeFormView &edge,
              state::PrimalFaceFormView face,
              double scale = 1.0);

void apply_d1(const state::PrimalEdgeFormStateView &edge,
              state::PrimalFaceFormStateView face,
              double scale = 1.0);

void apply_d2(const state::PrimalFaceFormView &face,
              state::PrimalThreeFormView cell);

// Subtract a scaled d1 from an existing 2-form. This is the Faraday update
// primitive, face -= dt * d1(edge).
void subtract_d1(const state::PrimalEdgeFormView &edge,
                 state::PrimalFaceFormView face,
                 double scale);

// Apply d1 transpose. The input must already have the metric/Hodge action
// appropriate for the requested codifferential.
void apply_d1_transpose(const state::DualEdgeFormView &dual_edge,
                        state::DualFaceFormView dual_face,
                        double scale = 1.0);

} // namespace mpcns::discrete
