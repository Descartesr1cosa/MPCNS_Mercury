// 2_Fields.cpp
#include "3_field/Field.h"

void Field::set_blocks(Grid *grd)
{
    storage_.bind_grid(grd);
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
