#pragma once

#include "mpcns/state/state_view.hpp"

namespace mpcns::state
{

// Storage-independent cochain views for new operators.  The legacy
// FieldBlock-backed definitions are isolated in compat/ while callers migrate.
struct PrimalZeroFormStateView { StateView node; int component = 0; };
struct PrimalEdgeFormStateView { StateView xi, eta, zeta; int component = 0; };
struct PrimalFaceFormStateView { StateView xi, eta, zeta; int component = 0; };
struct PrimalThreeFormStateView { StateView cell; int component = 0; };
struct DualEdgeFormStateView { StateView xi, eta, zeta; int component = 0; };
struct DualFaceFormStateView { StateView xi, eta, zeta; int component = 0; };
struct CellScalarStateView { StateView cell; int component = 0; };
struct CellVectorStateView { StateView cell; };

} // namespace mpcns::state
