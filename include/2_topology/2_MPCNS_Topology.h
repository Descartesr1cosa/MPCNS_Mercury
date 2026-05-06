// 我和谁连，在什么地方，怎么对齐（拓扑/连接）
// 把 1_grid 里“块之间的几何连接关系”抽出来，变成一组统一的接口描述（patch），供以后任何 field / halo 使用。
#pragma once
#include <vector>
#include <string>

#include "0_basic/TYPES.h"        // Int3, Box3
#include "1_grid/1_MPCNS_Grid.h"  // Grid, Block
#include "1_grid/Grid_Boundary.h" // Inner/Parallel/Physical_Boundary

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
    // nb[ perm[a] ] = sign[a] * local[ a ] + offset[a]
    struct IndexTransform
    {
        int perm[3]; // {0,1,2} 的排列，表示本地哪个坐标对应到邻居的 i/j/k
        int sign[3]; // +1 或 -1（二维时可以把第三维设成 0）
        Int3 offset; // 整数偏移（大多数情况为 0）
    };

    // 块-块接口（Inner + Parallel 都用这个）
    struct InterfacePatch
    {
        PatchKind kind;

        int this_rank; // this myid（0-based）
        int nb_rank;   // neighbor myid（Inner: 同 rank；Parallel: tar_myid）

        int this_block; // 0-based: Block index in this rank
        int nb_block;   // Inner: 对方 block index；Parallel: 如不知道可先设为 -1
        std::string this_block_name;
        std::string nb_block_name;

        // 在“节点 index 空间里的 box”，半开区间 [lo, hi)
        Box3 this_box_node; // 本块 interface 上的 node 区域
        Box3 nb_box_node;   // 对方块对应区域（Inner 可填，Parallel 暂时可留空或等于 this_box）

        // 本块接口面的方向，约定：
        //   -1 XMinus, +1 XPlus
        //   -2 YMinus, +2 YPlus
        //   -3 ZMinus, +3 ZPlus
        int direction = 0;

        // 邻居块对应接口面的方向，同样使用 ±1/±2/±3。
        // 对 Inner/Parallel 都应该填。
        // 对将来 coupling / validation 很有用。
        int nb_direction = 0;

        // 从本块 (i,j,k) 到邻居块 (i',j',k') 的映射
        IndexTransform trans;

        int32_t send_flag, recv_flag; // kind == PatchKind::Parallel才生效

        bool is_coupling; // 记录是否为耦合边界面
    };

    // 物理边界 patch：只在本块/本 rank 的一侧
    struct PhysicalPatch
    {
        int this_rank;
        int this_block;
        std::string this_block_name;

        int bc_id;           // boundary_num
        std::string bc_name; // boundary_name

        Box3 this_box_node; // 外边界在 node 空间中的区域 [lo,hi)
        int direction;      // ±1,±2,±3（和 Physical_Boundary 相同）

        const Physical_Boundary *raw = nullptr; // 回指原始结构（可选）
    };

    // 二维角区（edge strip）：至少两个坐标在 inner 外
    // 说明：这里不用 Direction，而是用 int dir1/dir2，约定：
    //   ±1 -> X±, ±2 -> Y±, ±3 -> Z±
    // 与 PhysicalPatch::direction 的约定保持一致，避免头文件循环依赖。
    struct EdgePatch
    {
        PatchKind kind; // Inner / Parallel / Physical

        int this_rank;
        int nb_rank;

        int this_block;
        int nb_block;
        std::string this_block_name;
        std::string nb_block_name;

        // 角区在 node 空间里的交集区域（沿着一条棱的 node strip）[lo, hi)
        Box3 this_box_node;
        Box3 nb_box_node;

        IndexTransform trans; // this -> nb

        // 该 edge 贴着两个出界方向，例如 XMinus(-1) + YPlus(+2)
        int dir1; // ±1, ±2, ±3
        int dir2; // ±1, ±2, ±3

        // 只有 Parallel 时才真正使用，用于 MPI 通信
        int32_t send_flag = 0;
        int32_t recv_flag = 0;

        bool is_coupling; // 记录是否为耦合边界面
    };

    // 三维角区（vertex 区）：三个坐标都在 inner 外
    struct VertexPatch
    {
        PatchKind kind; // Inner / Parallel / Physical

        int this_rank;
        int nb_rank;

        int this_block;
        int nb_block;
        std::string this_block_name;
        std::string nb_block_name;

        // 顶点附近的 node 区域（通常是一小块 [i,i+1)×[j,j+1)×[k,k+1)）
        Box3 this_box_node;
        Box3 nb_box_node;

        IndexTransform trans;

        // 三个出界方向，例如 XPlus(+1), YPlus(+2), ZMinus(-3)
        int dir1; // ±1, ±2, ±3
        int dir2; // ±1, ±2, ±3
        int dir3; // ±1, ±2, ±3

        int32_t send_flag = 0;
        int32_t recv_flag = 0;

        bool is_coupling; // 记录是否为耦合边界面
    };

    // 汇总：以后 3_field 只拿 Topology 这一个对象
    struct Topology
    {
        std::vector<InterfacePatch> inner_patches;
        std::vector<InterfacePatch> parallel_patches;
        std::vector<PhysicalPatch> physical_patches;

        // 二维角区（edge strip）
        std::vector<EdgePatch> inner_edge_patches;
        std::vector<EdgePatch> parallel_edge_patches;
        std::vector<EdgePatch> physical_edge_patches;

        // 三维角区（vertex 区）
        std::vector<VertexPatch> inner_vertex_patches;
        std::vector<VertexPatch> parallel_vertex_patches;
        std::vector<VertexPatch> physical_vertex_patches;
    };

    // 从 Grid 构造 topology（每个 rank 本地调用一次）
    Topology build_topology(Grid &grid, int my_rank, int dimension);

    void node_box_from_subsup(const int sub[3], const int sup[3], Box3 &box);

    // 从现有 inner/parallel/physical patch 自动生成 edge / vertex patch
    void build_edge_patches(Grid &grid, Topology &topo, int dimension);
    void build_vertex_patches(Grid &grid, Topology &topo, int dimension);

    // 在 build_edge_patches / build_vertex_patches 之后调用：
    // 把 coupling 的 InterfacePatch 追加为 PhysicalPatch，bc_name = prefix + nb_block_name
    void append_coupling_faces_as_physical_patches(Grid &grid,
                                                   Topology &topo,
                                                   int dimension,
                                                   const std::string &prefix = "Coupled-");
} // namespace TOPO
