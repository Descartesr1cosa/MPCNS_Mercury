// 2_Fields.cpp
#include "3_field/2_MPCNS_Field.h"

void Field::set_blocks(Grid *grd)
{
    blocks_.resize(grd->nblock);
    for (int i = 0; i < grd->nblock; i++)
        blocks_[i] = &(grd->grids(i));

    // 构建blocks_by_name_，e.g. 本进程中"Fluids"有哪些块
    blocks_by_name_.clear();
    for (int ib = 0; ib < (int)blocks_.size(); ++ib)
    {
        blocks_by_name_[blocks_[ib]->block_name].push_back(ib);
    }

    // 如果之前已经注册了场，就重新分配一次, for safety
    if (catalog_.size() > 0)
    {
        // 保证 field_blocks_ 维度正确
        field_blocks_.resize(catalog_.size());
        for (int32_t fid = 0; fid < catalog_.size(); ++fid)
        {
            field_blocks_[fid].resize(blocks_.size());
            allocate(fid);
        }
    }
}

void Field::register_field(const FieldDescriptor &desc)
{
    const int32_t old_size = catalog_.size();
    const int32_t fid = catalog_.add_or_get(desc);
    if (fid < old_size)
        return;

    // 保持 field_blocks_ 与 catalog_ 同步
    if (field_blocks_.size() != static_cast<size_t>(catalog_.size()))
        field_blocks_.resize(catalog_.size());

    // 如果已经有 blocks_，可以只为这个 field 分配一遍；
    if (!blocks_.empty())
    {
        field_blocks_[fid].resize(blocks_.size());
        allocate(fid);
    }

    return;
}

void Field::allocate(int32_t fieldID)
{
    const int nb = static_cast<int>(blocks_.size());
    const auto &desc = catalog_.descriptor(fieldID);

    std::vector<FieldBlock> &tmp = field_blocks_[fieldID];
    if (tmp.size() != nb)
        tmp.resize(nb);

    for (int b = 0; b < nb; ++b)
    {
        const Block &blk = *blocks_[b];
        const bool active = (desc.physics.empty() || blk.block_name == desc.physics);

        if (active)
            tmp[b].allocate(blk, desc);
        else
            tmp[b].bind_inactive(blk, desc);
    }
}

std::vector<std::string> Field::boundary_field_names() const
{
    return catalog_.boundary_field_names();
}

std::vector<std::string> Field::coupled_field_names() const
{
    return catalog_.coupled_field_names();
}

std::vector<FieldHaloRequest> Field::halo_requests() const
{
    return catalog_.halo_requests();
}
