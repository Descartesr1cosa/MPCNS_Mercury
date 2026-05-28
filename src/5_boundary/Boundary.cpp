#include "5_boundary/Boundary.h"
#include "0_basic/Error.h"

void BoundaryCore::SetUp(Grid *grd, Field *fld, TOPO::Topology *topo, Param *par, const std::vector<std::string> &field_names)
{
    grd_ = grd;
    fld_ = fld;
    topo_ = topo;
    par_ = par;

    if (!grd_ || !fld_ || !topo_ || !par_)
        ERROR::Abort("[BoundaryCore] SetUp got null pointer");

    // 1) 从 field_ids 推导需要的 locations
    enabled_locs_.clear();
    for (auto field : field_names)
    {
        int fid = fld_->field_id(field);
        const auto &desc = fld_->descriptor(fid);
        enabled_locs_.insert(desc.location);
    }

    // 2) 只为这些 locations 构建 Pattern（inner_slab，不扩 ghost）
    BuildPhysicalPatterns();
}
