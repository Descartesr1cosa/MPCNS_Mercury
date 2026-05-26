//=======================================================================================
// TOPOLOGY 模块：定义数据结构
//=======================================================================================
// TopologyTypes.h
//     数据层。
//     定义 patch、transform、Topology 容器。
//     不构建、不通信、不访问 Field。

// TopologyBuilder.h
//     构建入口层。
//     从 Grid 生成 Topology。

// TopologyOps.h
//     内部工具层。
//     做 box transform、direction ordering、patch priority 等纯操作。

// TopologyView.h
//     查询 / 适配层。
//     给 Halo、Boundary、Coupling 提供统一 patch view。

// Topology.h
//     等价类 / owner-alias 层。
//     处理 node/edge/face 物理同一 DOF、owner、alias、orientation sign。

// TopologyDebug.h
//     诊断层。
//     打印、检查、abort，不参与正常数据流。
//=======================================================================================
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
    // =========================================================================
    // File role
    // ---------
    // TopologyTypes.h defines the lightweight data records used by the topology
    // layer.
    //
    // This file answers:
    //   - Which block is connected to which neighbor block?
    //   - On which face / edge / vertex region are they connected?
    //   - How are node indices mapped from this block to the neighbor block?
    //   - Is a patch an inner, parallel, physical, or coupling interface?
    //
    // This file does NOT:
    //   - build topology from Grid;
    //   - perform MPI communication;
    //   - allocate Field data;
    //   - apply boundary conditions;
    //   - perform Halo pack/unpack.
    //
    // Important viewpoint convention:
    //   Every patch is stored from the viewpoint of `this_block`.
    //
    // Coupling convention:
    //   For InterfacePatch / EdgePatch / VertexPatch with is_coupling == true,
    //   the directed coupling pair is:
    //
    //       nb_block_name -> this_block_name
    //
    //   In other words, `this_block` is the local destination / receiver side
    //   of the coupling buffer.
    //
    // Halo redesign contract:
    //   These patches describe topology-side adjacency only.  A future halo
    //   pattern builder may combine codim-1 / codim-2 / codim-3 patches into
    //   one precomputed send/receive plan, with StorageAddress or HaloAddress
    //   describing actual field storage.  Topology does not create EntityId
    //   values for ordinary ghost slots, decide which fields need ghosts, or
    //   pack/unpack communication buffers.
    // =========================================================================

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

    // Codimension-1 face adjacency: a block interface available to a future
    // face communication pattern builder.  This is connectivity, not a halo
    // storage address and not an MPI operation.
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

    // Codimension-2 adjacency: edge / edge-strip / corner-strip connectivity.
    // A halo builder may use it directly for edge-region ghost items instead
    // of relying on face-to-edge propagation.
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

    // Codimension-3 adjacency: vertex / corner connectivity.  A future halo
    // plan may map it directly to corner StorageAddress/HaloAddress items.
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

} // namespace TOPO
