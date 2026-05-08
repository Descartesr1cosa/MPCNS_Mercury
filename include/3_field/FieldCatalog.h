#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "3_field/FieldDescriptor.h"

class FieldCatalog
{
public:
    bool has_field(const std::string &name) const;
    int32_t field_id(const std::string &name) const;
    int32_t size() const;

    const FieldDescriptor &descriptor(int32_t fid) const;
    const FieldDescriptor &descriptor(const std::string &field_name) const;
    const std::vector<FieldDescriptor> &descriptors() const;

    int32_t add_or_get(const FieldDescriptor &desc);

    std::vector<std::string> boundary_field_names() const;
    std::vector<std::string> coupled_field_names() const;
    std::vector<FieldHaloRequest> halo_requests() const;

private:
    std::vector<FieldDescriptor> descs_;
    std::unordered_map<std::string, int32_t> name_to_id_;
};
