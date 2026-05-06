#include "4_halo/1_MPCNS_Halo.h"
#include "0_basic/Error.h"

#include <vector>

const FieldHaloRequest &Halo::halo_request_(const std::string &field_name) const
{
    auto it = halo_registry_.find(field_name);
    if (it == halo_registry_.end())
        ERROR::Abort("[Halo] halo_request_: field is not registered: " + field_name);

    return it->second;
}

void Halo::sync_component_copy_field_(const std::string &field_name)
{
    const FieldHaloRequest &req = halo_request_(field_name);

    exchange_inner(field_name);
    exchange_parallel(field_name);

    if (static_cast<int>(req.level) >= static_cast<int>(HaloLevel::Edge))
    {
        exchange_inner_edge(field_name);
        exchange_parallel_edge(field_name);
    }

    if (static_cast<int>(req.level) >= static_cast<int>(HaloLevel::Vertex))
    {
        exchange_inner_vertex(field_name);
        exchange_parallel_vertex(field_name);
    }
}

void Halo::sync_component_copy_registered_()
{
    for (const std::string &field_name : component_copy_fields_)
        sync_component_copy_field_(field_name);
}

void Halo::sync_edge_1form_triplet_(const HaloTripletRequest &tri)
{
    if (tri.xi.empty() || tri.eta.empty() || tri.zeta.empty())
        ERROR::Abort("[Halo] sync_edge_1form_triplet_: incomplete triplet group: " + tri.group_name);

    const std::vector<std::string> fields{tri.xi, tri.eta, tri.zeta};
    data_trans_edge_1form_triplet(fields, tri.level);
}

void Halo::sync_edge_1form_triplets_registered_()
{
    for (const auto &kv : edge_1form_triplets_)
        sync_edge_1form_triplet_(kv.second);
}

void Halo::sync_face_2form_triplet_(const HaloTripletRequest &tri)
{
    if (tri.xi.empty() || tri.eta.empty() || tri.zeta.empty())
        ERROR::Abort("[Halo] sync_face_2form_triplet_: incomplete triplet group: " + tri.group_name);

    ERROR::Abort("[Halo] Face2FormTriplet sync is registered but not implemented yet: " + tri.group_name);
}

void Halo::sync_face_2form_triplets_registered_()
{
    for (const auto &kv : face_2form_triplets_)
        sync_face_2form_triplet_(kv.second);
}

void Halo::sync_owner_alias_request_(const HaloOwnerRequest &req)
{
    if (req.policy == OwnerSyncPolicy::None)
        return;

    ERROR::Abort("[Halo] OwnerAliasSync is registered but not implemented yet for field: " + req.field_name);
}

void Halo::sync_owner_alias_registered_()
{
    for (const auto &req : owner_sync_requests_)
        sync_owner_alias_request_(req);
}

bool Halo::field_is_component_copy_(const std::string &field_name) const
{
    for (const auto &name : component_copy_fields_)
    {
        if (name == field_name)
            return true;
    }

    return false;
}

void Halo::sync_registered()
{
    sync_component_copy_registered_();

    sync_edge_1form_triplets_registered_();

    sync_face_2form_triplets_registered_();

    sync_owner_alias_registered_();
}

void Halo::sync_field(const std::string &field_name)
{
    const FieldHaloRequest &req = halo_request_(field_name);

    const std::string group =
        req.sync_group.empty() ? req.field_name : req.sync_group;

    const HaloSyncSemantics sem = sync_semantics_(req);

    if (sem == HaloSyncSemantics::ComponentCopy)
    {
        sync_component_copy_field_(field_name);
    }
    else if (sem == HaloSyncSemantics::Edge1FormTriplet)
    {
        auto it = edge_1form_triplets_.find(group);
        if (it == edge_1form_triplets_.end())
            ERROR::Abort("[Halo] sync_field: missing edge 1-form group for field: " + field_name);

        sync_edge_1form_triplet_(it->second);
    }
    else if (sem == HaloSyncSemantics::Face2FormTriplet)
    {
        auto it = face_2form_triplets_.find(group);
        if (it == face_2form_triplets_.end())
            ERROR::Abort("[Halo] sync_field: missing face 2-form group for field: " + field_name);

        sync_face_2form_triplet_(it->second);
    }

    for (const auto &own : owner_sync_requests_)
    {
        if (own.field_name == field_name)
            sync_owner_alias_request_(own);
    }
}

void Halo::sync_group(const std::string &group_name)
{
    bool handled = false;

    auto eit = edge_1form_triplets_.find(group_name);
    if (eit != edge_1form_triplets_.end())
    {
        sync_edge_1form_triplet_(eit->second);
        handled = true;
    }

    auto fit = face_2form_triplets_.find(group_name);
    if (fit != face_2form_triplets_.end())
    {
        sync_face_2form_triplet_(fit->second);
        handled = true;
    }

    for (const auto &field_name : component_copy_fields_)
    {
        auto git = field_to_group_.find(field_name);

        if (git != field_to_group_.end() && git->second == group_name)
        {
            sync_component_copy_field_(field_name);
            handled = true;
        }
    }

    for (const auto &own : owner_sync_requests_)
    {
        if (own.sync_group == group_name)
        {
            sync_owner_alias_request_(own);
            handled = true;
        }
    }

    if (!handled)
        ERROR::Abort("[Halo] sync_group: unknown or empty group: " + group_name);
}
