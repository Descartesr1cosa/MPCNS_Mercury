#pragma once

#include "3_field/FieldBlock.h"

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

} // namespace mpcns::state
