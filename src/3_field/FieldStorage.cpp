#include "3_field/FieldStorage.h"

void FieldStorage::bind_grid(Grid *grid)
{
    grid_ = grid;

    blocks_.clear();
    blocks_by_name_.clear();

    if (!grid_)
    {
        field_blocks_.clear();
        return;
    }

    blocks_.resize(grid_->nblock);

    for (int ib = 0; ib < grid_->nblock; ++ib)
    {
        blocks_[ib] = &(grid_->grids(ib));
        blocks_by_name_[blocks_[ib]->block_name].push_back(ib);
    }
}

bool FieldStorage::has_grid() const
{
    return grid_ != nullptr;
}

int FieldStorage::num_blocks() const
{
    return static_cast<int>(blocks_.size());
}

Block &FieldStorage::block(int iblock)
{
    return *blocks_[iblock];
}

const Block &FieldStorage::block(int iblock) const
{
    return *blocks_[iblock];
}

void FieldStorage::resize_for_catalog(const FieldCatalog &catalog)
{
    if (field_blocks_.size() != static_cast<size_t>(catalog.size()))
        field_blocks_.resize(catalog.size());
}

void FieldStorage::allocate_field(int32_t fid, const FieldCatalog &catalog)
{
    resize_for_catalog(catalog);

    const int nb = static_cast<int>(blocks_.size());
    const FieldDescriptor &desc = catalog.descriptor(fid);

    std::vector<FieldBlock> &field_blocks_for_fid = field_blocks_[fid];

    if (field_blocks_for_fid.size() != static_cast<size_t>(nb))
        field_blocks_for_fid.resize(nb);

    for (int ib = 0; ib < nb; ++ib)
    {
        const Block &blk = *blocks_[ib];

        const bool active =
            desc.physics.empty() ||
            blk.block_name == desc.physics;

        if (active)
            field_blocks_for_fid[ib].allocate(blk, desc);
        else
            field_blocks_for_fid[ib].bind_inactive(blk, desc);
    }
}

void FieldStorage::allocate_all(const FieldCatalog &catalog)
{
    resize_for_catalog(catalog);
    for (int32_t fid = 0; fid < catalog.size(); ++fid)
        allocate_field(fid, catalog);
}

std::vector<FieldBlock> &FieldStorage::field(int32_t fid)
{
    return field_blocks_[fid];
}

const std::vector<FieldBlock> &FieldStorage::field(int32_t fid) const
{
    return field_blocks_[fid];
}

FieldBlock &FieldStorage::field(int32_t fid, int iblock)
{
    return field_blocks_[fid][iblock];
}

const FieldBlock &FieldStorage::field(int32_t fid, int iblock) const
{
    return field_blocks_[fid][iblock];
}

bool FieldStorage::field_active(int32_t fid, int iblock) const
{
    return field_blocks_[fid][iblock].is_allocated();
}

const std::unordered_map<std::string, std::vector<int>> &FieldStorage::blocks_by_name() const
{
    return blocks_by_name_;
}
