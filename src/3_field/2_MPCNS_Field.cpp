// 2_Fields.cpp
#include "3_field/2_MPCNS_Field.h"
#include "0_basic/Error.h"

int Field::field_id(const std::string &field_name) const
{
    auto it = name_to_id_.find(field_name);
    if (it == name_to_id_.end())
        ERROR::Abort("Field::field_id: unknown field name: " + field_name);
    return it->second;
}

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
    if (!field_descs_.empty())
    {
        // 保证 field_blocks_ 维度正确
        field_blocks_.resize(field_descs_.size());
        for (int32_t fid = 0; fid < (int32_t)field_descs_.size(); ++fid)
        {
            field_blocks_[fid].resize(blocks_.size());
            allocate(fid);
        }
    }
}

void Field::register_field(const FieldDescriptor &desc)
{
    auto existing = name_to_id_.find(desc.name);
    if (existing != name_to_id_.end())
    {
        const FieldDescriptor &old = field_descs_[existing->second];
        const bool same = old.location == desc.location &&
                          old.ncomp == desc.ncomp &&
                          old.nghost == desc.nghost &&
                          old.physics == desc.physics;
        if (same)
            return;

        ERROR::Abort("Field::register_field: duplicate field has different descriptor: " + desc.name);
    }

    int32_t fid = static_cast<int32_t>(field_descs_.size());
    field_descs_.push_back(desc);
    name_to_id_[desc.name] = fid;

    // 保持 field_blocks_ 与 field_descs_ 同步
    if (field_blocks_.size() != field_descs_.size())
        field_blocks_.resize(field_descs_.size());

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
    const auto &desc = field_descs_[fieldID];

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
