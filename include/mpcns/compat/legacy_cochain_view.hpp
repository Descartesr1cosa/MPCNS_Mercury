#pragma once

#include "3_field/FieldBlock.h"

namespace mpcns::state
{

// Compatibility-only FieldBlock-backed cochains. New public state APIs use
// the storage-independent views in mpcns/state/cochain_view.hpp.
struct PrimalZeroFormView { FieldBlock &node; int component = 0; };
struct PrimalEdgeFormView { FieldBlock &xi, &eta, &zeta; int component = 0; };
struct PrimalFaceFormView { FieldBlock &xi, &eta, &zeta; int component = 0; };
struct PrimalThreeFormView { FieldBlock &cell; int component = 0; };
struct DualEdgeFormView { FieldBlock &xi, &eta, &zeta; int component = 0; };
struct DualFaceFormView { FieldBlock &xi, &eta, &zeta; int component = 0; };
struct PrimalFaceScalarView { FieldBlock &xi, &eta, &zeta; int component = 0; };
struct PrimalEdgeScalarView { FieldBlock &xi, &eta, &zeta; int component = 0; };

} // namespace mpcns::state
