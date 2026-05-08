#include "2_topology/TopologyBuilder.h"
#include "2_topology/TopologyOps.h"

#include "0_basic/BoxOps.h"
#include "0_basic/Direction.h"
#include "0_basic/Error.h"

#include <map>
#include <algorithm>
#include <iostream>

namespace TOPO
{

    // 统一 inner / parallel / physical 三种面的信息（只在本 .cpp 里用）
    struct FaceOnBlock
    {
        PatchKind kind; // Inner / Parallel / Physical

        int this_rank;
        int nb_rank; // physical 时可以设成 this_rank 或 -1

        int this_block;
        int nb_block; // physical 时设 -1

        std::string this_block_name;
        std::string nb_block_name;

        Box3 this_box_node; // 在本块 Node 空间里的面区域
        Box3 nb_box_node;   // 对应到 nb_block 的 Node 面区域（physical 可以留空）

        IndexTransform trans; // this -> nb 的索引变换（physical 可以 identity）

        // 方向编码：±1/±2/±3 = X± / Y± / Z±
        int dir_code = 0;

        int32_t send_flag = 0;
        int32_t recv_flag = 0;

        bool is_coupling = false;
    };

    void build_edge_patches(Grid &grid, Topology &topo, int dimension)
    {
        topo.inner_edge_patches.clear();
        topo.parallel_edge_patches.clear();
        topo.physical_edge_patches.clear();

        if (dimension < 2)
            return;

        //==============================
        // 1. 按 this_block 把所有 face patch 收集成 FaceOnBlock
        //==============================

        std::map<int, std::vector<FaceOnBlock>> faces_of_block;

        // 1.1 inner + parallel 面
        auto append_interface_list =
            [&](const std::vector<InterfacePatch> &lst)
        {
            for (const auto &p : lst)
            {
                FaceOnBlock f;

                f.kind = p.kind; // PatchKind::Inner / Parallel
                f.this_rank = p.this_rank;
                f.nb_rank = p.nb_rank;
                f.this_block = p.this_block;
                f.nb_block = p.nb_block;
                f.this_block_name = p.this_block_name;
                f.nb_block_name = p.nb_block_name;
                f.this_box_node = p.this_box_node;
                f.nb_box_node = p.nb_box_node;
                f.trans = p.trans;
                f.send_flag = p.send_flag;
                f.recv_flag = p.recv_flag;

                f.is_coupling = p.is_coupling;

                f.dir_code = p.direction;

                if (!DIR::is_valid(f.dir_code))
                {
                    ERROR::Abort("[build_edge_patches] InterfacePatch has invalid direction");
                }

                faces_of_block[f.this_block].push_back(f);
            }
        };

        append_interface_list(topo.inner_patches);
        append_interface_list(topo.parallel_patches);

        // 1.2 physical 面
        auto append_physical_list =
            [&](const std::vector<PhysicalPatch> &lst)
        {
            for (const auto &p : lst)
            {
                FaceOnBlock f;

                f.kind = PatchKind::Physical;
                f.this_rank = p.this_rank;
                f.nb_rank = -1;
                f.this_block = p.this_block;
                f.nb_block = -1;
                f.this_block_name = p.this_block_name;
                f.nb_block_name = "";
                f.this_box_node = p.this_box_node;
                f.nb_box_node = Box3{};
                f.trans = IndexTransform{}; // identity
                f.send_flag = 0;
                f.recv_flag = 0;

                f.is_coupling = false; // 物理边界不是“跨物理耦合面”
                // 假设 PhysicalPatch 里有 int direction 字段（±1/±2/±3）
                f.dir_code = p.direction;

                if (!DIR::is_valid(f.dir_code))
                    ERROR::Abort("[build_edge_patches] PhysicalPatch has invalid direction");

                faces_of_block[f.this_block].push_back(f);
            }
        };

        append_physical_list(topo.physical_patches);

        //==============================
        // 2. 对每个 block，枚举“正交的面对” → 生成 EdgePatch
        //==============================

        for (auto &kv : faces_of_block)
        {
            int b = kv.first;
            auto &flist = kv.second;

            for (size_t i = 0; i < flist.size(); ++i)
            {
                for (size_t j = i + 1; j < flist.size(); ++j)
                {
                    const FaceOnBlock &f1 = flist[i];
                    const FaceOnBlock &f2 = flist[j];

                    if (f1.this_block != b || f2.this_block != b)
                        continue;

                    // 只要轴不同即可视为正交组合（X vs Y, Y vs Z, Z vs X）
                    if (!DIR::distinct_axes(f1.dir_code, f2.dir_code))
                        continue;

                    // 在本块 node 空间里的交集，就是这一条棱的 node strip
                    Box3 node_edge = BOX::intersect(f1.this_box_node, f2.this_box_node);
                    if (BOX::empty(node_edge))
                        continue;

                    //==============================
                    // 2.1 只用这“一对面”决定 edge 的 owner
                    //==============================
                    const FaceOnBlock *owner = &f1;
                    if (patch_priority(f2.kind, f2.is_coupling) >
                        patch_priority(f1.kind, f1.is_coupling))
                        owner = &f2;

                    //==============================
                    // 2.2 构造 EdgePatch
                    //==============================
                    EdgePatch ep;
                    ep.kind = owner->kind;
                    ep.this_rank = owner->this_rank;
                    ep.nb_rank = owner->nb_rank;
                    ep.this_block = owner->this_block;
                    ep.nb_block = owner->nb_block;
                    ep.this_block_name = owner->this_block_name;
                    ep.nb_block_name = owner->nb_block_name;
                    ep.is_coupling = owner->is_coupling;

                    // 本块上的棱 node 区域
                    ep.this_box_node = node_edge;

                    if (owner->kind == PatchKind::Inner ||
                        owner->kind == PatchKind::Parallel)
                    {
                        ep.nb_box_node = transform_node_box(node_edge, owner->trans);
                        ep.trans = owner->trans;
                        ep.send_flag = owner->send_flag;
                        ep.recv_flag = owner->recv_flag;
                    }
                    else
                    {
                        ep.nb_box_node = Box3{};
                        ep.trans = identity_transform(); // 物理边界 patch 的 transform 不会是未初始化语义
                        ep.send_flag = 0;
                        ep.recv_flag = 0;
                    }

                    // 方向：这条棱是由两个面 f1,f2 交出来的
                    ep.dir1 = f1.dir_code;
                    ep.dir2 = f2.dir_code;
                    // 需要保证，法向一定为dir1方向，便于发送者处理发送范围
                    order_edge_dirs_by_owner(owner->dir_code, ep.dir1, ep.dir2);

                    //==============================
                    // 2.3 按 kind 分类存到不同的 vector
                    //==============================
                    if (ep.kind == PatchKind::Inner)
                        topo.inner_edge_patches.push_back(ep);
                    else if (ep.kind == PatchKind::Parallel)
                        topo.parallel_edge_patches.push_back(ep);
                    else // Physical
                        topo.physical_edge_patches.push_back(ep);
                }
            }
        }
    }

    void build_vertex_patches(Grid &grid, Topology &topo, int dimension)
    {
        topo.inner_vertex_patches.clear();
        topo.parallel_vertex_patches.clear();
        topo.physical_vertex_patches.clear();

        // 只有 3D 情况才有真正的“三维角区”
        if (dimension < 3)
            return;

        //==============================
        // 1. 按 this_block 收集所有 edge
        //==============================

        std::map<int, std::vector<const EdgePatch *>> edges_of_block;

        for (const auto &ep : topo.inner_edge_patches)
            edges_of_block[ep.this_block].push_back(&ep);
        for (const auto &ep : topo.parallel_edge_patches)
            edges_of_block[ep.this_block].push_back(&ep);
        for (const auto &ep : topo.physical_edge_patches)
            edges_of_block[ep.this_block].push_back(&ep);

        // 避免同一 block、同一顶点重复建多个 VertexPatch
        std::map<int, std::vector<Int3>> used_vertex_node;

        //==============================
        // 2. 对每个 block，枚举 edge 与 edge 的交点 → Vertex
        //==============================

        for (auto &kv : edges_of_block)
        {
            int b = kv.first;
            auto &elist = kv.second;

            if (elist.size() < 2)
                continue;

            for (size_t i = 0; i < elist.size(); ++i)
            {
                for (size_t j = i + 1; j < elist.size(); ++j)
                {
                    const EdgePatch &e1 = *elist[i];
                    const EdgePatch &e2 = *elist[j];

                    if (e1.this_block != b || e2.this_block != b)
                        continue;

                    // 同一方向上的两段 edge（比如 0..10 和 10..mz 的 Z-edge）不构成顶点
                    int ax1 = edge_axis1_from_face_dirs(e1.dir1, e1.dir2);
                    int ax2 = edge_axis1_from_face_dirs(e2.dir1, e2.dir2);
                    if (ax1 == ax2)
                        continue;

                    // 在本块 node 空间里的交集
                    Box3 node_vert = BOX::intersect(e1.this_box_node, e2.this_box_node);
                    if (BOX::empty(node_vert))
                        continue;

                    int len_i = node_vert.hi.i - node_vert.lo.i;
                    int len_j = node_vert.hi.j - node_vert.lo.j;
                    int len_k = node_vert.hi.k - node_vert.lo.k;
                    if (!(len_i == 1 && len_j == 1 && len_k == 1))
                        continue; // 必须是一个 node

                    // 去重：同一 block、同一 node 只生成一个 vertex
                    auto &used_list = used_vertex_node[b];
                    bool duplicate = false;
                    for (const auto &v : used_list)
                    {
                        if (v.i == node_vert.lo.i &&
                            v.j == node_vert.lo.j &&
                            v.k == node_vert.lo.k)
                        {
                            duplicate = true;
                            break;
                        }
                    }
                    if (duplicate)
                        continue;
                    used_list.push_back(node_vert.lo);

                    //==============================
                    // 2.1 收集所有覆盖到这个顶点的 edge（候选）
                    //==============================
                    std::vector<const EdgePatch *> candidates;
                    candidates.push_back(&e1);
                    candidates.push_back(&e2);

                    for (size_t k = 0; k < elist.size(); ++k)
                    {
                        const EdgePatch &ek = *elist[k];
                        if (ek.this_block != b)
                            continue;
                        if (&ek == &e1 || &ek == &e2)
                            continue;

                        Box3 inter3 = BOX::intersect(ek.this_box_node, node_vert);
                        if (!BOX::empty(inter3))
                            candidates.push_back(&ek);
                    }

                    //==============================
                    // 2.2 按优先级 Inner > Parallel > Physical 选 owner edge
                    //==============================
                    const EdgePatch *owner = candidates[0];
                    for (size_t k = 1; k < candidates.size(); ++k)
                    {
                        if (patch_priority(candidates[k]->kind, candidates[k]->is_coupling) >
                            patch_priority(owner->kind, owner->is_coupling))
                            owner = candidates[k];
                    }

                    //==============================
                    // 2.3 构造 VertexPatch（基本信息）
                    //==============================
                    VertexPatch vp;
                    vp.kind = owner->kind;
                    vp.this_rank = owner->this_rank;
                    vp.nb_rank = owner->nb_rank;
                    vp.this_block = owner->this_block;
                    vp.nb_block = owner->nb_block;
                    vp.this_block_name = owner->this_block_name;
                    vp.nb_block_name = owner->nb_block_name;
                    vp.is_coupling = owner->is_coupling;

                    vp.this_box_node = node_vert;

                    if (owner->kind == PatchKind::Inner ||
                        owner->kind == PatchKind::Parallel)
                    {
                        vp.nb_box_node = transform_node_box(node_vert, owner->trans);
                        vp.trans = owner->trans;
                        vp.send_flag = owner->send_flag;
                        vp.recv_flag = owner->recv_flag;
                    }
                    else
                    {
                        vp.nb_box_node = Box3{};
                        vp.trans = identity_transform();
                        vp.send_flag = 0;
                        vp.recv_flag = 0;
                    }

                    //==============================
                    // 2.4 从 candidates 里统计 3 个面方向 dir1/dir2/dir3
                    //==============================
                    int dir_by_axis[4] = {0, 0, 0, 0}; // 1..3 用来存 X/Y/Z 方向

                    auto try_add_dir = [&](int dir)
                    {
                        if (!DIR::is_valid(dir))
                            return;

                        const int ax = DIR::axis1(dir);

                        if (dir_by_axis[ax] == 0)
                            dir_by_axis[ax] = dir;
                    };

                    for (const EdgePatch *ep : candidates)
                    {
                        try_add_dir(ep->dir1);
                        try_add_dir(ep->dir2);
                    }

                    // 填到 VertexPatch
                    vp.dir1 = dir_by_axis[1]; // X±
                    vp.dir2 = dir_by_axis[2]; // Y±
                    vp.dir3 = dir_by_axis[3]; // Z±

                    // 需要保证，法向一定为dir1方向，便于发送者处理发送范围
                    order_vertex_dirs_by_owner(owner->dir1, vp.dir1, vp.dir2, vp.dir3);

                    //==============================
                    // 2.5 按 kind 分类存储
                    //==============================
                    if (vp.kind == PatchKind::Inner)
                        topo.inner_vertex_patches.push_back(vp);
                    else if (vp.kind == PatchKind::Parallel)
                        topo.parallel_vertex_patches.push_back(vp);
                    else
                        topo.physical_vertex_patches.push_back(vp);
                }
            }
        }
    }

} // namespace TOPO
