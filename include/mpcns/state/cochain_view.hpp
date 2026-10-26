#pragma once

#include "3_field/FieldBlock.h"
#include "mpcns/state/state_view.hpp"

namespace mpcns::state
{

// Transitional typed views over the existing FieldBlock storage. A future
// StateLayout may provide the same views without exposing FieldBlock.
struct PrimalZeroFormView
{
    FieldBlock &node;
    int component = 0;
};

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

struct PrimalThreeFormView
{
    FieldBlock &cell;
    int component = 0;
};

// The dual complex is stored on complementary primal locations during the
// transition. Distinct types preserve mathematical meaning even when storage
// shape is currently identical.
struct DualEdgeFormView
{
    FieldBlock &xi;
    FieldBlock &eta;
    FieldBlock &zeta;
    int component = 0;
};

struct DualFaceFormView
{
    FieldBlock &xi;
    FieldBlock &eta;
    FieldBlock &zeta;
    int component = 0;
};

// Scalar geometric data live at staggered locations but are not cochains.
struct PrimalFaceScalarView
{
    FieldBlock &xi;
    FieldBlock &eta;
    FieldBlock &zeta;
    int component = 0;
};

struct PrimalEdgeScalarView
{
    FieldBlock &xi;
    FieldBlock &eta;
    FieldBlock &zeta;
    int component = 0;
};

// Storage-independent cochain views for new operators.  The legacy
// FieldBlock-backed views above remain transitional adapters while callers
// migrate one production path at a time.
struct PrimalZeroFormStateView { StateView node; int component = 0; };
struct PrimalEdgeFormStateView { StateView xi, eta, zeta; int component = 0; };
struct PrimalFaceFormStateView { StateView xi, eta, zeta; int component = 0; };
struct PrimalThreeFormStateView { StateView cell; int component = 0; };
struct DualEdgeFormStateView { StateView xi, eta, zeta; int component = 0; };
struct DualFaceFormStateView { StateView xi, eta, zeta; int component = 0; };
struct CellScalarStateView { StateView cell; int component = 0; };
struct CellVectorStateView { StateView cell; };

} // namespace mpcns::state
