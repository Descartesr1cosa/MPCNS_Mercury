#include "6_boundary/Boundary.h"

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
    switch (loc)
    {
    case StaggerLocation::Cell:
        return {1, 1, 1};
    case StaggerLocation::Node:
        return {0, 0, 0};
    case StaggerLocation::FaceXi:
        return {0, 1, 1};
    case StaggerLocation::FaceEt:
        return {1, 0, 1};
    case StaggerLocation::FaceZe:
        return {1, 1, 0};
    case StaggerLocation::EdgeXi:
        return {1, 0, 0};
    case StaggerLocation::EdgeEt:
        return {0, 1, 0};
    case StaggerLocation::EdgeZe:
        return {0, 0, 1};
    default:
        return {0, 0, 0};
    }
}

Int3 BoundaryCore::LocInnerHi(const Block &blk, StaggerLocation loc)
{
    // blk.mx/my/mz 是 cell counts (Ni,Nj,Nk)
    const Int3 nodes = {blk.mx + 1, blk.my + 1, blk.mz + 1};
    const Int3 d = LocDelta(loc);
    return {nodes.i - d.i, nodes.j - d.j, nodes.k - d.k}; // half-open [0,hi)
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
    const int ax = std::abs(direction); // 1/2/3
    const int sgn = (direction > 0) ? +1 : -1;

    const Int3 hi_in = LocInnerHi(blk, loc);
    const Int3 d = LocDelta(loc);

    // 法向轴以外的两条切向轴
    int t1, t2;
    if (ax == 1)
    {
        t1 = 2;
        t2 = 3;
    }
    else if (ax == 2)
    {
        t1 = 1;
        t2 = 3;
    }
    else
    {
        t1 = 1;
        t2 = 2;
    }

    Box3 b{};

    // 1) 先填切向范围：node_box -> loc_box
    auto set_tangent = [&](int t)
    {
        int lo, hi;
        if (t == 1)
            ConvertTangent(face_node_box.lo.i, face_node_box.hi.i, d.i, lo, hi);
        if (t == 2)
            ConvertTangent(face_node_box.lo.j, face_node_box.hi.j, d.j, lo, hi);
        if (t == 3)
            ConvertTangent(face_node_box.lo.k, face_node_box.hi.k, d.k, lo, hi);

        if (t == 1)
        {
            b.lo.i = lo;
            b.hi.i = hi;
        }
        if (t == 2)
        {
            b.lo.j = lo;
            b.hi.j = hi;
        }
        if (t == 3)
        {
            b.lo.k = lo;
            b.hi.k = hi;
        }
    };
    set_tangent(t1);
    set_tangent(t2);

    // 2) 再填法向方向：域内贴边 1 层 slab
    if (ax == 1)
    {
        if (sgn < 0)
        {
            b.lo.i = 0;
            b.hi.i = 1;
        }
        else
        {
            b.lo.i = hi_in.i - 1;
            b.hi.i = hi_in.i;
        }
    }
    else if (ax == 2)
    {
        if (sgn < 0)
        {
            b.lo.j = 0;
            b.hi.j = 1;
        }
        else
        {
            b.lo.j = hi_in.j - 1;
            b.hi.j = hi_in.j;
        }
    }
    else
    {
        if (sgn < 0)
        {
            b.lo.k = 0;
            b.hi.k = 1;
        }
        else
        {
            b.lo.k = hi_in.k - 1;
            b.hi.k = hi_in.k;
        }
    }

    return b;
}

Box3 BoundaryCore::MakeGhostSlabFromInner(const Box3 &inner_slab,
                                          int direction,
                                          int nghost)
{
    const int ax = std::abs(direction);
    const int sgn = (direction > 0) ? +1 : -1;

    Box3 g = inner_slab; // 切向范围直接继承

    if (ax == 1)
    {
        if (sgn < 0)
        {
            g.hi.i = inner_slab.lo.i;
            g.lo.i = inner_slab.lo.i - nghost;
        }
        else
        {
            g.lo.i = inner_slab.hi.i;
            g.hi.i = inner_slab.hi.i + nghost;
        }
    }
    else if (ax == 2)
    {
        if (sgn < 0)
        {
            g.hi.j = inner_slab.lo.j;
            g.lo.j = inner_slab.lo.j - nghost;
        }
        else
        {
            g.lo.j = inner_slab.hi.j;
            g.hi.j = inner_slab.hi.j + nghost;
        }
    }
    else
    {
        if (sgn < 0)
        {
            g.hi.k = inner_slab.lo.k;
            g.lo.k = inner_slab.lo.k - nghost;
        }
        else
        {
            g.lo.k = inner_slab.hi.k;
            g.hi.k = inner_slab.hi.k + nghost;
        }
    }

    return g;
}