#pragma once

#include "1_grid/1_MPCNS_Grid.h"

namespace GRID_TRAITS
{
    inline Int3 cell_counts(const Block &blk)
    {
        return {blk.mx, blk.my, blk.mz};
    }
}
