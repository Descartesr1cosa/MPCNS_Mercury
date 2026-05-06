#include "3_field/FieldStorage.h"

void FieldStorage::bind_grid(Grid *grid)
{
    grid_ = grid;

    blocks_.resize(grid_->nblock);
    for (int i = 0; i < grid_->nblock; ++i)
        blocks_[i] = &(grid_->grids(i));

    blocks_by_name_.clear();
    for (int ib = 0; ib < static_cast<int>(blocks_.size()); ++ib)
        blocks_by_name_[blocks_[ib]->block_name].push_back(ib);
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

    std::vector<FieldBlock> &blocks = field_blocks_[fid];
    if (blocks.size() != static_cast<size_t>(nb))
        blocks.resize(nb);

    for (int b = 0; b < nb; ++b)
    {
        const Block &blk = *blocks_[b];
        const bool active = (desc.physics.empty() || blk.block_name == desc.physics);

        if (active)
            blocks[b].allocate(blk, desc);
        else
            blocks[b].bind_inactive(blk, desc);
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
