#include "0_basic/Error.h"
#include "2_topology/TopologyBuildDetail.h"
#include "1_grid/1_MPCNS_Grid.h"
#include <unordered_set>
#include <sstream>
#include <stdexcept>
namespace TOPO
{
    namespace detail
    {
    // Coupling interfaces reuse their established face direction and box.
    // Build a key only to avoid inserting an identical physical patch twice.
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
}
