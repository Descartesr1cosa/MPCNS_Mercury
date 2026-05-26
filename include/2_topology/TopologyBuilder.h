// 我和谁连，在什么地方，怎么对齐（拓扑/连接）
// 把 1_grid 里“块之间的几何连接关系”抽出来，变成一组统一的接口描述（patch），供以后任何 field / halo 使用。
#pragma once

#include "2_topology/Topology.h"

#include "1_grid/1_MPCNS_Grid.h"
#include "1_grid/Grid_Boundary.h"

namespace TOPO
{
    // Stable construction entry: build patches and entity equivalence in one call.
    Topology build_topology(Grid &grid, int my_rank, int dimension);
} // namespace TOPO
