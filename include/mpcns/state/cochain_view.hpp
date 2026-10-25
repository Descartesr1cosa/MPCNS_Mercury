#pragma once

#include "3_field/FieldBlock.h"

namespace mpcns::state
{

// Transitional typed views over the existing three directional FieldBlock
// storage.  A future StateLayout will provide the same views without exposing
// FieldBlock to numerical operators.
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

// The dual complex is stored on the complementary primal locations during the
// transition: dual edges live on primal faces and dual faces live on primal
// edges.  The distinct types preserve the mathematical meaning even though
// their current backend storage has the same shape.
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

// Scalar geometric data live at a staggered location but are not cochains.
// These types prevent a diagonal Hodge scale from being confused with a field
// value at a call site.
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

} // namespace mpcns::state
