#include "0_basic/Error.h"
#include "2_topology/2_MPCNS_Topology.h"
#include <unordered_set>
#include <sstream>
#include <stdexcept>
namespace TOPO
{
    // // 推断一个 face patch 的 direction：±1/±2/±3
    // // 依据：this_box_node 在某一维厚度为 1，且贴在该维的 0 或 max（max=imax/jmax/kmax）
    // static int infer_face_direction_from_node_box(const Block &blk,
    //                                               const Box3 &b,
    //                                               int dimension)
    // {
    //     auto extent_i = b.hi.i - b.lo.i;
    //     auto extent_j = b.hi.j - b.lo.j;
    //     auto extent_k = b.hi.k - b.lo.k;

    //     // 只在 active 维度里找法向：2D 只考虑 i/j；3D 考虑 i/j/k
    //     if (dimension >= 1 && extent_i == 1)
    //     {
    //         if (b.lo.i == 0)
    //             return -1;
    //         if (b.lo.i == blk.mx)
    //             return +1;
    //     }
    //     if (dimension >= 2 && extent_j == 1)
    //     {
    //         if (b.lo.j == 0)
    //             return -2;
    //         if (b.lo.j == blk.my)
    //             return +2;
    //     }
    //     if (dimension >= 3 && extent_k == 1)
    //     {
    //         if (b.lo.k == 0)
    //             return -3;
    //         if (b.lo.k == blk.mz)
    //             return +3;
    //     }

    //     std::ostringstream oss;
    //     oss << "[append_coupling_faces_as_physical_patches] cannot infer direction. "
    //         << "dim=" << dimension
    //         << " box=[(" << b.lo.i << "," << b.lo.j << "," << b.lo.k << ")->("
    //         << b.hi.i << "," << b.hi.j << "," << b.hi.k << ")] "
    //         << "blk(max)=(" << blk.imax << "," << blk.jmax << "," << blk.kmax << ")";
    //     ERROR::Abort(oss.str());
    // }

    // 生成一个去重 key：block + dir + bc_name + box
    static std::string make_physical_key(const PhysicalPatch &p)
    {
        std::ostringstream oss;
        oss << p.this_block << "|"
            << p.direction << "|"
            << p.bc_name << "|"
            << p.this_box_node.lo.i << "," << p.this_box_node.lo.j << "," << p.this_box_node.lo.k << "|"
            << p.this_box_node.hi.i << "," << p.this_box_node.hi.j << "," << p.this_box_node.hi.k;
        return oss.str();
    }

    void append_coupling_faces_as_physical_patches(Grid &grid,
                                                   Topology &topo,
                                                   int dimension,
                                                   const std::string &prefix)
    {
        // 记录已有 physical patch，避免重复插入
        std::unordered_set<std::string> existed;
        existed.reserve(topo.physical_patches.size() * 2 + 64);
        for (const auto &p : topo.physical_patches)
            existed.insert(make_physical_key(p));

        int new_bc_id = -1000000; // 给新增 Coupled-* 一个负 id（仅用于标识/调试）

        auto append_one = [&](const InterfacePatch &iface)
        {
            if (!iface.is_coupling)
                return;

            if (iface.nb_block_name.empty())
                ERROR::Abort("[append_coupling_faces_as_physical_patches] nb_block_name is empty");

            const Block &blk = grid.grids(iface.this_block);

            PhysicalPatch p;
            p.this_rank = iface.this_rank;
            p.this_block = iface.this_block;
            p.this_block_name = iface.this_block_name;

            p.bc_id = new_bc_id--;
            p.bc_name = prefix + iface.nb_block_name;

            p.this_box_node = iface.this_box_node;
            p.direction = iface.direction;

            if (p.direction == 0)
                ERROR::Abort("[append_coupling_faces_as_physical_patches] iface.direction is zero");

            p.raw = nullptr; // 不修改 PhysicalPatch 类型：Coupled-* 没有对应 Physical_Boundary*

            const auto key = make_physical_key(p);
            if (existed.insert(key).second)
                topo.physical_patches.push_back(std::move(p));
        };

        // 从 inner/parallel 接口里把 coupling 面转成 PhysicalPatch
        for (const auto &iface : topo.inner_patches)
            append_one(iface);
        for (const auto &iface : topo.parallel_patches)
            append_one(iface);

        // 注意：这里不重排 topo.physical_patches，以免改变你原有物理边界覆盖顺序。
        // 如果你确实需要 Coupled-* 参与 Priority 排序，你可以在调用方再跑一次 stable_sort。
    }
}