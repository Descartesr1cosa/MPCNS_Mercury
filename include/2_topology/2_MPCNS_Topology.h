// 我和谁连，在什么地方，怎么对齐（拓扑/连接）
// 把 1_grid 里“块之间的几何连接关系”抽出来，变成一组统一的接口描述（patch），供以后任何 field / halo 使用。
#pragma once

#include <string>

#include "2_topology/TopologyTypes.h"

#include "1_grid/1_MPCNS_Grid.h"
#include "1_grid/Grid_Boundary.h"

namespace TOPO
{
    // 从 Grid 构造 topology，每个 rank 本地调用一次
    Topology build_topology(Grid &grid, int my_rank, int dimension);

    // 旧接口保留，内部已转调 TopologyOps。
    void node_box_from_subsup(const int sub[3], const int sup[3], Box3 &box);

    // 从现有 inner/parallel/physical patch 自动生成 edge / vertex patch
    void build_edge_patches(Grid &grid, Topology &topo, int dimension);
    void build_vertex_patches(Grid &grid, Topology &topo, int dimension);

    // 在 build_edge_patches / build_vertex_patches 之后调用：
    // 把 coupling 的 InterfacePatch 追加为 PhysicalPatch，
    // bc_name = prefix + nb_block_name
    void append_coupling_faces_as_physical_patches(Grid &grid,
                                                   Topology &topo,
                                                   int dimension,
                                                   const std::string &prefix = "Coupled-");
} // namespace TOPO