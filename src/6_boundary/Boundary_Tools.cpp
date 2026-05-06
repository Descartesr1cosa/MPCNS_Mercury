#include "6_boundary/Boundary.h"
#include "0_basic/LayoutTraits.h"
#include "1_grid/BlockTraits.h"

// ------------------------------------------------------------
// Physical patch cache
// ------------------------------------------------------------
void BoundaryCore::BuildPhysicalPatterns()
{
    phy_patterns_.clear();

    for (auto loc : enabled_locs_)
    {
        BOUND::PhysicalPattern pat;
        pat.location = loc;
        phy_patterns_[loc] = std::move(pat);
    }

    for (const auto &p : topo_->physical_patches)
    {
        const int ib = p.this_block;
        const Block &blk = grd_->grids(ib); // 你工程里取 block 的接口按实际改

        const Box3 face_node_box = p.this_box_node; // topo 给的 node patch box（半开区间）
        const int dir = p.direction;

        for (auto loc : enabled_locs_)
        {
            BOUND::PhysicalRegion r;
            r.this_block = p.this_block;
            r.this_block_name = p.this_block_name;
            r.bc_id = p.bc_id;
            r.bc_name = p.bc_name;
            r.direction = dir;
            r.raw = p.raw;

            // cycle 可选：也可完全由 direction 推导
            if (p.raw)
            {
                r.cycle.i = p.raw->cycle[0];
                r.cycle.j = p.raw->cycle[1];
                r.cycle.k = p.raw->cycle[2];
            }

            // 关键：Build 阶段推导“域内贴边1层”的 inner_slab（loc 坐标）
            r.inner_slab = MakeInnerSlabBox_OneLayer(blk, loc, face_node_box, dir);

            phy_patterns_[loc].regions.push_back(std::move(r));
        }
    }
}

// ------------------------------------------------------------
// Geometry helpers (same style as coupling buffer build)
// ------------------------------------------------------------
Int3 BoundaryCore::LocDelta(StaggerLocation loc)
{
    return LAYOUT::dof_delta(loc);
}

Int3 BoundaryCore::LocInnerHi(const Block &blk, StaggerLocation loc)
{
    return LAYOUT::owned_box_from_cells(GRID_TRAITS::cell_counts(blk), loc).hi;
}

void BoundaryCore::ConvertTangent(int lo_n, int hi_n, int delta, int &lo, int &hi)
{
    lo = lo_n;
    hi = (delta == 0) ? hi_n : (hi_n - 1);
}

Box3 BoundaryCore::MakeInnerSlabBox_OneLayer(const Block &blk,
                                             StaggerLocation loc,
                                             const Box3 &face_node_box,
                                             int direction)
{
    return LAYOUT::boundary_inner_slab_one_layer_from_cells(
        GRID_TRAITS::cell_counts(blk),
        loc,
        face_node_box,
        direction);
}

Box3 BoundaryCore::MakeGhostSlabFromInner(const Box3 &inner_slab,
                                          int direction,
                                          int nghost)
{
    return LAYOUT::ghost_slab_from_inner(inner_slab, direction, nghost);
}
