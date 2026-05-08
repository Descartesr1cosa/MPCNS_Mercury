#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/FieldBlock.h"
#include "3_field/FieldCatalog.h"

class FieldStorage
{
public:
    void bind_grid(Grid *grid);
    bool has_grid() const;
    int num_blocks() const;

    Block &block(int iblock);
    const Block &block(int iblock) const;

    void resize_for_catalog(const FieldCatalog &catalog);
    void allocate_field(int32_t fid, const FieldCatalog &catalog);
    void allocate_all(const FieldCatalog &catalog);

    std::vector<FieldBlock> &field(int32_t fid);
    const std::vector<FieldBlock> &field(int32_t fid) const;
    FieldBlock &field(int32_t fid, int iblock);
    const FieldBlock &field(int32_t fid, int iblock) const;

    bool field_active(int32_t fid, int iblock) const;

    const std::unordered_map<std::string, std::vector<int>> &blocks_by_name() const;

private:
    Grid *grid_ = nullptr;
    std::vector<Block *> blocks_;
    std::unordered_map<std::string, std::vector<int>> blocks_by_name_;
    std::vector<std::vector<FieldBlock>> field_blocks_;
};
