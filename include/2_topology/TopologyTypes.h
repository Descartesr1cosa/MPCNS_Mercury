#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "0_basic/TYPES.h"

// 为了暂时保留 PhysicalPatch::raw，不在这里 include Grid_Boundary.h。
// 后续建议彻底删除 raw 指针。
class Physical_Boundary;

namespace TOPO
{
    // 接口类型：内部同 rank、跨 rank、物理外边界
    enum class PatchKind
    {
        Inner,
        Parallel,
        Physical
    };

    // point_local -> point_nb 的索引变换：
    // nb[ perm[a] ] = sign[a] * local[a] + offset[a]
    struct IndexTransform
    {
        int perm[3]; // {0,1,2} 的排列，表示本地哪个坐标对应到邻居的 i/j/k
        int sign[3]; // +1 或 -1；二维时第三维可为 0
        Int3 offset; // 整数偏移
    };

    // 块-块接口：Inner + Parallel 共用
    struct InterfacePatch
    {
        PatchKind kind;

        int this_rank;
        int nb_rank;

        int this_block;
        int nb_block;

        std::string this_block_name;
        std::string nb_block_name;

        // node-space half-open box: [lo, hi)
        Box3 this_box_node;
        Box3 nb_box_node;

        // 本块接口面的方向：
        //   -1 XMinus, +1 XPlus
        //   -2 YMinus, +2 YPlus
        //   -3 ZMinus, +3 ZPlus
        int direction = 0;

        // 邻居块对应接口面的方向
        int nb_direction = 0;

        // this -> neighbor 的 node index 映射
        IndexTransform trans;

        // Parallel 时生效
        int32_t send_flag = 0;
        int32_t recv_flag = 0;

        // 是否是多物理场耦合面
        bool is_coupling = false;
    };

    // 物理边界 patch：只在本块/本 rank 的一侧
    struct PhysicalPatch
    {
        int this_rank = 0;
        int this_block = 0;
        std::string this_block_name;

        int bc_id = 0;
        std::string bc_name;

        Box3 this_box_node;
        int direction = 0;

        // 临时兼容旧逻辑。
        // 后续建议删除 raw，使 TopologyTypes 完全不引用 Grid_Boundary 概念。
        const Physical_Boundary *raw = nullptr;
    };

    // 二维角区 / edge strip
    struct EdgePatch
    {
        PatchKind kind;

        int this_rank = 0;
        int nb_rank = 0;

        int this_block = 0;
        int nb_block = 0;

        std::string this_block_name;
        std::string nb_block_name;

        Box3 this_box_node;
        Box3 nb_box_node;

        IndexTransform trans;

        // 该 edge 贴着两个出界方向，例如 XMinus(-1) + YPlus(+2)
        int dir1 = 0;
        int dir2 = 0;

        int32_t send_flag = 0;
        int32_t recv_flag = 0;

        bool is_coupling = false;
    };

    // 三维角区 / vertex patch
    struct VertexPatch
    {
        PatchKind kind;

        int this_rank = 0;
        int nb_rank = 0;

        int this_block = 0;
        int nb_block = 0;

        std::string this_block_name;
        std::string nb_block_name;

        Box3 this_box_node;
        Box3 nb_box_node;

        IndexTransform trans;

        int dir1 = 0;
        int dir2 = 0;
        int dir3 = 0;

        int32_t send_flag = 0;
        int32_t recv_flag = 0;

        bool is_coupling = false;
    };

    // 汇总 topology
    struct Topology
    {
        std::vector<InterfacePatch> inner_patches;
        std::vector<InterfacePatch> parallel_patches;
        std::vector<PhysicalPatch> physical_patches;

        std::vector<EdgePatch> inner_edge_patches;
        std::vector<EdgePatch> parallel_edge_patches;
        std::vector<EdgePatch> physical_edge_patches;

        std::vector<VertexPatch> inner_vertex_patches;
        std::vector<VertexPatch> parallel_vertex_patches;
        std::vector<VertexPatch> physical_vertex_patches;
    };

} // namespace TOPO