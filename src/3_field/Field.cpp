#include "3_field/Field.h"

Field::Field(Grid *grid, Param *param, int geometry_ghost)
    : grd(grid), par(param)
{
    bind_grid(grid);
    build_geometry(geometry_ghost);
}

void Field::bind_grid(Grid *grid)
{
    storage_.bind_grid(grid);
    storage_.allocate_all(catalog_);
}

void Field::register_field(const FieldDescriptor &desc)
{
    const int32_t old_size = catalog_.size();
    const int32_t fid = catalog_.add_or_get(desc);
    if (fid < old_size)
        return;

    if (storage_.has_grid())
        storage_.allocate_field(fid, catalog_);

    return;
}

void Field::allocate(int32_t fieldID)
{
    storage_.allocate_field(fieldID, catalog_);
}

int Field::field_id(const std::string &field_name) const
{
    return catalog_.field_id(field_name);
}

bool Field::has_field(const std::string &name) const
{
    return catalog_.has_field(name);
}

const FieldDescriptor &Field::descriptor(int32_t fid) const
{
    return catalog_.descriptor(fid);
}

const FieldDescriptor &Field::descriptor(const std::string &field_name) const
{
    return catalog_.descriptor(field_name);
}

const std::vector<FieldDescriptor> &Field::descriptors() const
{
    return catalog_.descriptors();
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

std::vector<FieldBlock> &Field::field(int32_t fid)
{
    return storage_.field(fid);
}

const std::vector<FieldBlock> &Field::field(int32_t fid) const
{
    return storage_.field(fid);
}

FieldBlock &Field::field(int32_t fid, int iblock)
{
    return storage_.field(fid, iblock);
}

const FieldBlock &Field::field(int32_t fid, int iblock) const
{
    return storage_.field(fid, iblock);
}

std::vector<FieldBlock> &Field::field(const std::string &name)
{
    return field(field_id(name));
}

const std::vector<FieldBlock> &Field::field(const std::string &name) const
{
    return field(field_id(name));
}

FieldBlock &Field::field(const std::string &name, int iblock)
{
    return field(field_id(name), iblock);
}

const FieldBlock &Field::field(const std::string &name, int iblock) const
{
    return field(field_id(name), iblock);
}
