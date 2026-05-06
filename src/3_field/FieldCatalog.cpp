#include "3_field/FieldCatalog.h"

#include "0_basic/Error.h"

namespace
{
    bool same_descriptor(const FieldDescriptor &old, const FieldDescriptor &desc)
    {
        return old.location == desc.location &&
               old.ncomp == desc.ncomp &&
               old.nghost == desc.nghost &&
               old.physics == desc.physics &&
               old.value_kind == desc.value_kind &&
               old.sync.group == desc.sync.group &&
               old.sync.do_coupling == desc.sync.do_coupling &&
               old.sync.do_physical == desc.sync.do_physical &&
               old.sync.do_halo == desc.sync.do_halo &&
               old.sync.halo_level == desc.sync.halo_level &&
               old.sync.owner_sync == desc.sync.owner_sync &&
               old.sync.orientation_aware == desc.sync.orientation_aware;
    }
}

bool FieldCatalog::has_field(const std::string &name) const
{
    return name_to_id_.find(name) != name_to_id_.end();
}

int32_t FieldCatalog::field_id(const std::string &name) const
{
    auto it = name_to_id_.find(name);
    if (it == name_to_id_.end())
        ERROR::Abort("FieldCatalog::field_id: unknown field name: " + name);
    return it->second;
}

int32_t FieldCatalog::size() const
{
    return static_cast<int32_t>(descs_.size());
}

const FieldDescriptor &FieldCatalog::descriptor(int32_t fid) const
{
    return descs_[fid];
}

const FieldDescriptor &FieldCatalog::descriptor(const std::string &field_name) const
{
    return descriptor(field_id(field_name));
}

const std::vector<FieldDescriptor> &FieldCatalog::descriptors() const
{
    return descs_;
}

int32_t FieldCatalog::add_or_get(const FieldDescriptor &desc)
{
    auto existing = name_to_id_.find(desc.name);
    if (existing != name_to_id_.end())
    {
        const FieldDescriptor &old = descs_[existing->second];
        if (same_descriptor(old, desc))
            return existing->second;

        ERROR::Abort("FieldCatalog::add_or_get: duplicate field has different descriptor: " + desc.name);
    }

    const int32_t fid = static_cast<int32_t>(descs_.size());
    descs_.push_back(desc);
    name_to_id_[desc.name] = fid;
    return fid;
}

std::vector<std::string> FieldCatalog::boundary_field_names() const
{
    std::vector<std::string> names;
    for (const auto &desc : descs_)
    {
        if (desc.sync.do_physical)
            names.push_back(desc.name);
    }
    return names;
}

std::vector<std::string> FieldCatalog::coupled_field_names() const
{
    std::vector<std::string> names;
    for (const auto &desc : descs_)
    {
        if (desc.sync.do_coupling)
            names.push_back(desc.name);
    }
    return names;
}

std::vector<FieldHaloRequest> FieldCatalog::halo_requests() const
{
    std::vector<FieldHaloRequest> requests;
    for (const auto &desc : descs_)
    {
        if (!desc.sync.do_halo)
            continue;

        FieldHaloRequest req;
        req.field_name = desc.name;
        req.sync_group = desc.sync.group.empty() ? desc.name : desc.sync.group;
        req.location = desc.location;
        req.value_kind = desc.value_kind;
        req.ncomp = desc.ncomp;
        req.nghost = desc.nghost;
        req.level = desc.sync.halo_level;
        req.owner_sync = desc.sync.owner_sync;
        req.orientation_aware = desc.sync.orientation_aware;
        requests.push_back(req);
    }
    return requests;
}
