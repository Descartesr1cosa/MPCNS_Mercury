#include "1_Boundary.h"
#include "0_basic/Error.h"

void MercuryBoundary::Setup(Grid *grd, Field *fld, TOPO::Topology *topo, Halo *halo, Param *par, const std::vector<std::string> &boundary_fields)
{
    grd_ = grd;
    fld_ = fld;
    topo_ = topo;
    halo_ = halo;
    par_ = par;
    boundary_fields_ = boundary_fields;

    if (!grd_ || !fld_ || !topo_ || !halo_ || !par_)
        ERROR::Abort("MercuryBoundary::Setup: null pointer(s)");

    // 1) 初始化 BoundaryCore：Build 阶段会按 location 缓存每个 patch 的 inner_slab（法向1层）
    bound_.SetUp(grd_, fld_, topo_, par_, boundary_fields_);

    InitBCStateFromParam_();

    InstallHandlers();

    InstallDefaultGroups();

    built_ = false;

    Build(true);
}

void MercuryBoundary::Sync(const std::string &group_name)
{
    if (!built_)
        ERROR::Abort("MercuryBoundary::Sync: call Build() first");

    auto it = groups_.find(group_name);
    if (it == groups_.end())
        ERROR::Abort(("MercuryBoundary::Sync: unknown group: " + group_name).c_str());

    Sync_(it->second);
}
