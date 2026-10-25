#pragma once

#include <cstdint>

namespace mpcns::core
{

// A logical structured-grid index, independent of any particular legacy
// array or field implementation.
struct LogicalIndex
{
    std::int32_t i = 0;
    std::int32_t j = 0;
    std::int32_t k = 0;
};

// Half-open logical extent [lo, hi).  It is deliberately separate from the
// legacy Box3 so new low-level APIs do not acquire a 0_basic dependency.
struct LogicalBox
{
    LogicalIndex lo;
    LogicalIndex hi;

    bool contains(LogicalIndex index) const
    {
        return index.i >= lo.i && index.i < hi.i &&
               index.j >= lo.j && index.j < hi.j &&
               index.k >= lo.k && index.k < hi.k;
    }
};

} // namespace mpcns::core
