#include "5_boundary/Boundary.h"

#include "0_basic/Error.h"
#include "2_topology/Topology.h"
#include "4_halo/detail/HaloBuildBoxMakers.h"

#include <algorithm>
#include <cstdlib>

namespace
{
    Direction physical_direction_from_int(int d)
    {
        switch (d)
        {
        case -1:
            return Direction::XMinus;
        case 1:
            return Direction::XPlus;
        case -2:
            return Direction::YMinus;
        case 2:
            return Direction::YPlus;
        case -3:
            return Direction::ZMinus;
        case 3:
            return Direction::ZPlus;
        default:
            ERROR::Abort("[BoundaryCore] physical corner default: invalid direction.");
            return Direction::XMinus;
        }
    }

    int clamp_physical_index(int v, int lo, int hi)
    {
        return std::max(lo, std::min(v, hi - 1));
    }

    void fill_physical_corner_ghost(FieldBlock &U,
                                    const Int3 &hi_in,
                                    const Box3 &ghost_box)
    {
        if (!U.is_allocated())
            return;

        const int ncomp = U.descriptor().ncomp;

        for (int i = ghost_box.lo.i; i < ghost_box.hi.i; ++i)
            for (int j = ghost_box.lo.j; j < ghost_box.hi.j; ++j)
                for (int k = ghost_box.lo.k; k < ghost_box.hi.k; ++k)
                {
                    const int ii = clamp_physical_index(i, 0, hi_in.i);
                    const int jj = clamp_physical_index(j, 0, hi_in.j);
                    const int kk = clamp_physical_index(k, 0, hi_in.k);

                    for (int m = 0; m < ncomp; ++m)
                        U(i, j, k, m) = U(ii, jj, kk, m);
                }
    }
}

void BoundaryCore::RegisterPhysical(StaggerLocation loc,
                                    const std::string &field_name,
                                    const std::string &bc_name,
                                    BOUND::PhysicalHandler h)
{
    phy_reg_[BOUND::PhysicalKey{loc, field_name, bc_name}] = std::move(h);
}

void BoundaryCore::RegisterPhysical(const std::string &field_name,
                                    const std::string &bc_name,
                                    BOUND::PhysicalHandler h)
{
    if (!fld_)
        ERROR::Abort("[BoundaryCore] RegisterPhysical called before SetUp");

    const int fid = fld_->field_id(field_name);
    const auto &desc = fld_->descriptor(fid);
    const StaggerLocation loc = desc.location;

    if (enabled_locs_.find(loc) == enabled_locs_.end())
        ERROR::Abort("[BoundaryCore] RegisterPhysical: location not enabled by SetUp(field_ids).");

    RegisterPhysical(loc, field_name, bc_name, std::move(h));
}

void BoundaryCore::CheckPhysicalHandlers(const std::vector<std::string> &field_names) const
{
    if (!fld_ || !topo_ || !grd_)
        ERROR::Abort("[BoundaryCore] CheckPhysicalHandlers called before SetUp");

    std::map<StaggerLocation, std::map<std::string, const BOUND::PhysicalRegion *>> bcset;
    for (const auto &[loc, pat] : phy_patterns_)
    {
        for (const auto &r : pat.regions)
        {
            if (!bcset[loc].count(r.bc_name))
                bcset[loc][r.bc_name] = &r;
        }
    }

    struct Missing
    {
        std::string field;
        StaggerLocation loc{};
        std::string bc;
        int example_block = -1;
        int bc_id = -1;
        int direction = 0;
    };

    std::vector<Missing> missing;

    for (const auto &field_name : field_names)
    {
        const int fid = fld_->field_id(field_name);
        const auto &desc = fld_->descriptor(fid);
        const StaggerLocation loc = desc.location;

        if (phy_patterns_.find(loc) == phy_patterns_.end())
            ERROR::Abort("[BoundaryCore] No physical pattern built for location of field: " + field_name);

        auto bit = bcset.find(loc);
        if (bit == bcset.end())
            continue;

        for (const auto &[bc_name, sample] : bit->second)
        {
            BOUND::PhysicalKey key{loc, field_name, bc_name};
            if (phy_reg_.find(key) == phy_reg_.end())
            {
                missing.push_back(Missing{
                    field_name, loc, bc_name,
                    sample ? sample->this_block : -1,
                    sample ? sample->bc_id : -1,
                    sample ? sample->direction : 0});
            }
        }
    }

    if (missing.empty())
        return;

    std::string msg = "[BoundaryCore] Missing Physical BC handlers:\n";
    for (const auto &m : missing)
    {
        msg += "  field=" + m.field +
               " loc=" + std::to_string((int)m.loc) +
               " bc_name=" + m.bc +
               " (example: block=" + std::to_string(m.example_block) +
               " bc_id=" + std::to_string(m.bc_id) +
               " dir=" + std::to_string(m.direction) + ")\n";
    }
    ERROR::Abort(msg);
    std::abort();
}

void BoundaryCore::ApplyPhysical(const std::string &field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    const StaggerLocation loc = desc.location;
    const int nghost = desc.nghost;

    auto pit = phy_patterns_.find(loc);
    if (pit == phy_patterns_.end())
        ERROR::Abort("[BoundaryCore] ApplyPhysical: pattern not built for this field/location: " + field_name);

    for (const auto &cached : pit->second.regions)
    {
        FieldBlock &U = fld_->field(fid, cached.this_block);
        if (!U.is_allocated())
            continue;

        BOUND::PhysicalRegion work = cached;
        work.box = MakeGhostSlabFromInner(cached.inner_slab, cached.direction, nghost);

        auto h = ResolvePhysical(loc, field_name, cached.bc_name);
        if (!h)
            ERROR::Abort("[BoundaryCore] ApplyPhysical: missing handler for field=" + field_name + " bc=" + cached.bc_name);
        h(U, fld_, work, nghost);
    }
}

void BoundaryCore::ApplyPhysical(const std::vector<std::string> &field_names)
{
    for (const auto &fn : field_names)
        ApplyPhysical(fn);
}

void BoundaryCore::ApplyPhysical(const std::string &field_name, HaloLevel stage)
{
    if (stage == HaloLevel::Corner1D)
        ApplyPhysical(field_name);
    else if (stage == HaloLevel::Corner2D)
        ApplyPhysicalEdgeDefault(field_name);
    else if (stage == HaloLevel::Corner3D)
        ApplyPhysicalVertexDefault(field_name);
}

void BoundaryCore::ApplyPhysical(const std::vector<std::string> &field_names, HaloLevel stage)
{
    for (const auto &fn : field_names)
        ApplyPhysical(fn, stage);
}

void BoundaryCore::ApplyPhysicalEdgeDefault(const std::string &field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    const StaggerLocation loc = desc.location;
    const int nghost = desc.nghost;

    for (const auto &ep : topo_->physical_edge_patches)
    {
        FieldBlock &U = fld_->field(fid, ep.this_block);
        const Block &blk = grd_->grids(ep.this_block);
        const Box3 g = HALO_BOX::make_2DCorner_ghost_box(
            loc,
            ep.this_box_node,
            physical_direction_from_int(ep.dir1),
            physical_direction_from_int(ep.dir2),
            nghost);
        fill_physical_corner_ghost(U, LocInnerHi(blk, loc), g);
    }
}

void BoundaryCore::ApplyPhysicalEdgeDefault(const std::vector<std::string> &field_names)
{
    for (const auto &fn : field_names)
        ApplyPhysicalEdgeDefault(fn);
}

void BoundaryCore::ApplyPhysicalVertexDefault(const std::string &field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    const StaggerLocation loc = desc.location;
    const int nghost = desc.nghost;

    for (const auto &vp : topo_->physical_vertex_patches)
    {
        FieldBlock &U = fld_->field(fid, vp.this_block);
        const Block &blk = grd_->grids(vp.this_block);
        const Box3 g = HALO_BOX::make_3DCorner_ghost_box(
            loc,
            vp.this_box_node,
            physical_direction_from_int(vp.dir1),
            physical_direction_from_int(vp.dir2),
            physical_direction_from_int(vp.dir3),
            nghost);
        fill_physical_corner_ghost(U, LocInnerHi(blk, loc), g);
    }
}

void BoundaryCore::ApplyPhysicalVertexDefault(const std::vector<std::string> &field_names)
{
    for (const auto &fn : field_names)
        ApplyPhysicalVertexDefault(fn);
}

void BoundaryCore::ApplyPhysicalCornerDefault(const std::string &field_name)
{
    ApplyPhysicalEdgeDefault(field_name);
    ApplyPhysicalVertexDefault(field_name);
}

void BoundaryCore::ApplyPhysicalCornerDefault(const std::vector<std::string> &field_names)
{
    for (const auto &fn : field_names)
        ApplyPhysicalCornerDefault(fn);
}

BOUND::PhysicalHandler BoundaryCore::ResolvePhysical(StaggerLocation loc,
                                                     const std::string &field_name,
                                                     const std::string &bc_name) const
{
    BOUND::PhysicalKey key{loc, field_name, bc_name};
    auto it = phy_reg_.find(key);
    if (it != phy_reg_.end())
        return it->second;
    return nullptr;
}

void BoundaryCore::DefaultPhysicalCopy(FieldBlock &U, Field * /*fld*/, const BOUND::PhysicalRegion &r, int nghost)
{
    if (!U.is_allocated())
        return;

    const Box3 &inner = r.inner_slab;
    const int ax = std::abs(r.direction);
    const int sgn = (r.direction > 0) ? +1 : -1;

    const int i_ref = (ax == 1) ? ((sgn < 0) ? inner.lo.i : (inner.hi.i - 1)) : 0;
    const int j_ref = (ax == 2) ? ((sgn < 0) ? inner.lo.j : (inner.hi.j - 1)) : 0;
    const int k_ref = (ax == 3) ? ((sgn < 0) ? inner.lo.k : (inner.hi.k - 1)) : 0;

    const int ncomp = U.descriptor().ncomp;
    const Box3 g = MakeGhostSlabFromInner(inner, r.direction, nghost);

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                const int ii = (ax == 1) ? i_ref : i;
                const int jj = (ax == 2) ? j_ref : j;
                const int kk = (ax == 3) ? k_ref : k;

                for (int m = 0; m < ncomp; ++m)
                    U(i, j, k, m) = U(ii, jj, kk, m);
            }
}
