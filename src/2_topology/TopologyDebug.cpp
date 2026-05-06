
#include "2_topology/TopologyDebug.h"

#include "0_basic/BoxOps.h"
#include "0_basic/Direction.h"
#include "0_basic/Error.h"
#include "0_basic/LayoutTraits.h"
#include "1_grid/BlockTraits.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <iostream>
#include <sstream>

namespace
{
    Box3 node_domain_of_block(const Block &blk)
    {
        const Int3 ncells = GRID_TRAITS::cell_counts(blk);
        const Int3 nnodes = LAYOUT::node_size_from_cells(ncells);

        return Box3{{0, 0, 0}, {nnodes.i, nnodes.j, nnodes.k}};
    }

    std::string patch_where(const char *kind, int block, const Box3 &box)
    {
        std::ostringstream oss;
        oss << kind << " block=" << block
            << " box=" << BOX::to_string(box);
        return oss.str();
    }

    void validate_block_id(int ib, int nblock, const char *where)
    {
        if (ib < 0 || ib >= nblock)
        {
            std::ostringstream oss;
            oss << where << ": invalid block id " << ib
                << ", nblock=" << nblock;
            ERROR::Abort(oss.str());
        }
    }

    void validate_node_box_inside_block(const Block &blk,
                                        const Box3 &box,
                                        const char *where)
    {
        const Box3 domain = node_domain_of_block(blk);
        if (!BOX::contains(domain, box))
        {
            std::ostringstream oss;
            oss << where
                << ": node box outside block node domain. "
                << "box=" << BOX::to_string(box)
                << " domain=" << BOX::to_string(domain);
            ERROR::Abort(oss.str());
        }

        BOX::assert_nonempty(box, where);
    }

    void validate_transform(const TOPO::IndexTransform &tr,
                            int dimension,
                            const char *where)
    {
        bool used[3] = {false, false, false};

        for (int d = 0; d < 3; ++d)
        {
            if (tr.perm[d] < 0 || tr.perm[d] > 2)
            {
                std::ostringstream oss;
                oss << where << ": invalid perm[" << d << "]=" << tr.perm[d];
                ERROR::Abort(oss.str());
            }

            if (used[tr.perm[d]])
            {
                std::ostringstream oss;
                oss << where << ": duplicated perm value " << tr.perm[d];
                ERROR::Abort(oss.str());
            }

            used[tr.perm[d]] = true;
        }

        for (int d = 0; d < 3; ++d)
        {
            const int s = tr.sign[d];

            if (d < dimension)
            {
                if (s != -1 && s != +1)
                {
                    std::ostringstream oss;
                    oss << where << ": invalid active sign[" << d << "]=" << s;
                    ERROR::Abort(oss.str());
                }
            }
            else
            {
                // 2D 情况下第三维可能为 0，也允许 ±1。
                if (s != -1 && s != 0 && s != +1)
                {
                    std::ostringstream oss;
                    oss << where << ": invalid inactive sign[" << d << "]=" << s;
                    ERROR::Abort(oss.str());
                }
            }
        }
    }

    void validate_interface_patch(const TOPO::InterfacePatch &p,
                                  Grid &grid,
                                  int dimension,
                                  const char *list_name,
                                  int index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "]";

        validate_block_id(p.this_block, grid.nblock, where.str().c_str());

        const Block &blk = grid.grids(p.this_block);

        validate_node_box_inside_block(
            blk,
            p.this_box_node,
            patch_where(where.str().c_str(), p.this_block, p.this_box_node).c_str());

        if (!DIR::is_valid(p.direction))
        {
            std::ostringstream oss;
            oss << where.str() << ": invalid direction=" << p.direction;
            ERROR::Abort(oss.str());
        }

        if (!DIR::is_valid(p.nb_direction))
        {
            std::ostringstream oss;
            oss << where.str() << ": invalid nb_direction=" << p.nb_direction;
            ERROR::Abort(oss.str());
        }

        validate_transform(p.trans, dimension, where.str().c_str());
    }

    void validate_physical_patch(const TOPO::PhysicalPatch &p,
                                 Grid &grid,
                                 const char *list_name,
                                 int index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "]";

        validate_block_id(p.this_block, grid.nblock, where.str().c_str());

        const Block &blk = grid.grids(p.this_block);

        validate_node_box_inside_block(
            blk,
            p.this_box_node,
            patch_where(where.str().c_str(), p.this_block, p.this_box_node).c_str());

        if (!DIR::is_valid(p.direction))
        {
            std::ostringstream oss;
            oss << where.str() << ": invalid direction=" << p.direction;
            ERROR::Abort(oss.str());
        }
    }

    void validate_edge_patch(const TOPO::EdgePatch &p,
                             Grid &grid,
                             const char *list_name,
                             int index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "]";

        validate_block_id(p.this_block, grid.nblock, where.str().c_str());

        const Block &blk = grid.grids(p.this_block);

        validate_node_box_inside_block(
            blk,
            p.this_box_node,
            patch_where(where.str().c_str(), p.this_block, p.this_box_node).c_str());

        if (!DIR::distinct_axes(p.dir1, p.dir2))
        {
            std::ostringstream oss;
            oss << where.str() << ": edge directions share the same axis: "
                << p.dir1 << ", " << p.dir2;
            ERROR::Abort(oss.str());
        }
    }

    void validate_vertex_patch(const TOPO::VertexPatch &p,
                               Grid &grid,
                               const char *list_name,
                               int index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "]";

        validate_block_id(p.this_block, grid.nblock, where.str().c_str());

        const Block &blk = grid.grids(p.this_block);

        validate_node_box_inside_block(
            blk,
            p.this_box_node,
            patch_where(where.str().c_str(), p.this_block, p.this_box_node).c_str());

        if (!DIR::distinct_axes(p.dir1, p.dir2, p.dir3))
        {
            std::ostringstream oss;
            oss << where.str() << ": vertex directions are not three distinct axes: "
                << p.dir1 << ", " << p.dir2 << ", " << p.dir3;
            ERROR::Abort(oss.str());
        }
    }
}

namespace TOPO_DEBUG
{
    void dump_topology_summary(const TOPO::Topology &topo, int my_rank)
    {
        if (my_rank != 0)
            return;

        std::cout << "\n========== Topology summary ==========\n";
        std::cout << "faces:\n";
        std::cout << "  inner    = " << topo.inner_patches.size() << "\n";
        std::cout << "  parallel = " << topo.parallel_patches.size() << "\n";
        std::cout << "  physical = " << topo.physical_patches.size() << "\n";

        std::cout << "edges:\n";
        std::cout << "  inner    = " << topo.inner_edge_patches.size() << "\n";
        std::cout << "  parallel = " << topo.parallel_edge_patches.size() << "\n";
        std::cout << "  physical = " << topo.physical_edge_patches.size() << "\n";

        std::cout << "vertices:\n";
        std::cout << "  inner    = " << topo.inner_vertex_patches.size() << "\n";
        std::cout << "  parallel = " << topo.parallel_vertex_patches.size() << "\n";
        std::cout << "  physical = " << topo.physical_vertex_patches.size() << "\n";
        std::cout << "======================================\n\n";
    }

    void validate_topology_or_abort(TOPO::Topology &topo,
                                    Grid &grid,
                                    int my_rank,
                                    int dimension)
    {
        (void)my_rank;

        for (int n = 0; n < static_cast<int>(topo.inner_patches.size()); ++n)
            validate_interface_patch(topo.inner_patches[n], grid, dimension, "inner_patches", n);

        for (int n = 0; n < static_cast<int>(topo.parallel_patches.size()); ++n)
            validate_interface_patch(topo.parallel_patches[n], grid, dimension, "parallel_patches", n);

        for (int n = 0; n < static_cast<int>(topo.physical_patches.size()); ++n)
            validate_physical_patch(topo.physical_patches[n], grid, "physical_patches", n);

        for (int n = 0; n < static_cast<int>(topo.inner_edge_patches.size()); ++n)
            validate_edge_patch(topo.inner_edge_patches[n], grid, "inner_edge_patches", n);

        for (int n = 0; n < static_cast<int>(topo.parallel_edge_patches.size()); ++n)
            validate_edge_patch(topo.parallel_edge_patches[n], grid, "parallel_edge_patches", n);

        for (int n = 0; n < static_cast<int>(topo.physical_edge_patches.size()); ++n)
            validate_edge_patch(topo.physical_edge_patches[n], grid, "physical_edge_patches", n);

        for (int n = 0; n < static_cast<int>(topo.inner_vertex_patches.size()); ++n)
            validate_vertex_patch(topo.inner_vertex_patches[n], grid, "inner_vertex_patches", n);

        for (int n = 0; n < static_cast<int>(topo.parallel_vertex_patches.size()); ++n)
            validate_vertex_patch(topo.parallel_vertex_patches[n], grid, "parallel_vertex_patches", n);

        for (int n = 0; n < static_cast<int>(topo.physical_vertex_patches.size()); ++n)
            validate_vertex_patch(topo.physical_vertex_patches[n], grid, "physical_vertex_patches", n);
    }
}