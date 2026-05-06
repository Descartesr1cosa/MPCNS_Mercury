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

bool Halo::halo_level_includes_edge_(HaloLevel level) const
{
    return level == HaloLevel::Edge ||
           level == HaloLevel::Vertex;
}

bool Halo::halo_level_includes_vertex_(HaloLevel level) const
{
    return level == HaloLevel::Vertex;
}

void Halo::sync_component_copy_field_(const std::string &field_name)
{
    const FieldHaloRequest &req = halo_request_(field_name);

    exchange_inner(field_name);
    exchange_parallel(field_name);

    if (halo_level_includes_edge_(req.level))
    {
        exchange_inner_edge(field_name);
        exchange_parallel_edge(field_name);
    }

    if (halo_level_includes_vertex_(req.level))
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

    if (tri.value_kind != FieldValueKind::EdgeCovariant1Form)
        ERROR::Abort("[Halo] sync_edge_1form_triplet_: group is not EdgeCovariant1Form: " + tri.group_name);

    if (!tri.orientation_aware)
        ERROR::Abort("[Halo] sync_edge_1form_triplet_: group is not orientation-aware: " + tri.group_name);

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

    if (tri.value_kind != FieldValueKind::FaceContravariant2Form)
        ERROR::Abort("[Halo] sync_face_2form_triplet_: group is not FaceContravariant2Form: " + tri.group_name);

    if (!tri.orientation_aware)
        ERROR::Abort("[Halo] sync_face_2form_triplet_: group is not orientation-aware: " + tri.group_name);

    for (const std::string *field_name : {&tri.xi, &tri.eta, &tri.zeta})
    {
        const FieldHaloRequest &req = halo_request_(*field_name);
        const std::string group = req.sync_group.empty() ? req.field_name : req.sync_group;
        if (group != tri.group_name ||
            req.nghost != tri.nghost ||
            req.level != tri.level ||
            req.value_kind != FieldValueKind::FaceContravariant2Form ||
            !req.orientation_aware)
        {
            ERROR::Abort("[Halo] sync_face_2form_triplet_: inconsistent face 2-form triplet request: " + tri.group_name);
        }
    }

    sync_face_2form_triplet_face_level_(tri);
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

    require_owner_equiv_available_(req.policy, req.field_name);

    const int fid = fld_->field_id(req.field_name);
    const TOPO::EquivDofKind kind = owner_policy_to_equiv_kind_(req.policy);
    const auto &classes = equiv_->classes(kind);

    int my_rank = 0;
    PARALLEL::mpi_rank(&my_rank);

    for (const auto &cls : classes)
    {
        const TOPO::EquivMember &owner = cls.owner;

        if (!owner_member_matches_field_(req, owner))
            continue;

        bool has_local_alias = false;
        for (const auto &member : cls.members)
        {
            if (member.rank == my_rank &&
                !member.is_owner &&
                owner_member_matches_field_(req, member))
            {
                has_local_alias = true;
                break;
            }
        }

        if (!has_local_alias)
            continue;

        if (owner.rank != my_rank)
        {
            ERROR::Abort("[Halo] cross-rank OwnerAliasSync is not implemented yet for field: " +
                         req.field_name);
        }

        for (const auto &alias : cls.members)
        {
            if (alias.rank != my_rank)
                continue;

            if (alias.is_owner)
                continue;

            if (!owner_member_matches_field_(req, alias))
                continue;

            const int sign = owner_alias_sign_(req, owner, alias);
            copy_owner_to_alias_local_(req, fid, owner, alias, sign);
        }
    }
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

TOPO::EquivDofKind Halo::owner_policy_to_equiv_kind_(OwnerSyncPolicy policy) const
{
    switch (policy)
    {
    case OwnerSyncPolicy::NodeOwner:
        return TOPO::EquivDofKind::Node;
    case OwnerSyncPolicy::EdgeOwner:
        return TOPO::EquivDofKind::Edge;
    case OwnerSyncPolicy::FaceOwner:
        return TOPO::EquivDofKind::Face;
    case OwnerSyncPolicy::None:
        break;
    }

    ERROR::Abort("[Halo] owner_policy_to_equiv_kind_: invalid OwnerSyncPolicy");
    return TOPO::EquivDofKind::Node;
}

void Halo::require_owner_equiv_available_(OwnerSyncPolicy policy,
                                          const std::string &field_name) const
{
    if (policy == OwnerSyncPolicy::None)
        return;

    if (!equiv_)
    {
        ERROR::Abort("[Halo] OwnerAliasSync requested but TopologyEquiv is not bound. field=" +
                     field_name);
    }

    if (policy == OwnerSyncPolicy::NodeOwner && !equiv_->has_node_equiv())
    {
        ERROR::Abort("[Halo] NodeOwner sync requested but node equivalence is not built. field=" +
                     field_name);
    }

    if (policy == OwnerSyncPolicy::EdgeOwner && !equiv_->has_edge_equiv())
    {
        ERROR::Abort("[Halo] EdgeOwner sync requested but edge equivalence is not built. field=" +
                     field_name);
    }

    if (policy == OwnerSyncPolicy::FaceOwner && !equiv_->has_face_equiv())
    {
        ERROR::Abort("[Halo] FaceOwner sync requested but face equivalence is not built. field=" +
                     field_name);
    }
}

int Halo::owner_alias_sign_(const HaloOwnerRequest &req,
                            const TOPO::EquivMember &owner,
                            const TOPO::EquivMember &alias) const
{
    if (!req.orientation_aware)
        return +1;

    if (req.policy == OwnerSyncPolicy::EdgeOwner &&
        req.value_kind == FieldValueKind::EdgeCovariant1Form)
    {
        // orient_sign is relative to the canonical orientation. For signs in
        // {+1, -1}, alias / owner is equivalent to alias * owner.
        return alias.orient_sign * owner.orient_sign;
    }

    if (req.policy == OwnerSyncPolicy::FaceOwner &&
        req.value_kind == FieldValueKind::FaceContravariant2Form)
    {
        // orient_sign is relative to the canonical orientation. For signs in
        // {+1, -1}, alias / owner is equivalent to alias * owner.
        return alias.orient_sign * owner.orient_sign;
    }

    return +1;
}

bool Halo::owner_member_matches_field_(const HaloOwnerRequest &req,
                                       const TOPO::EquivMember &member) const
{
    return member.location == req.location;
}

void Halo::copy_owner_to_alias_local_(const HaloOwnerRequest &req,
                                      int fid,
                                      const TOPO::EquivMember &owner,
                                      const TOPO::EquivMember &alias,
                                      int sign)
{
    FieldBlock &owner_block = fld_->field(fid, owner.block);
    FieldBlock &alias_block = fld_->field(fid, alias.block);

    if (!owner_block.is_allocated() || !alias_block.is_allocated())
        return;

    for (int m = 0; m < req.ncomp; ++m)
    {
        alias_block(alias.i, alias.j, alias.k, m) =
            static_cast<double>(sign) * owner_block(owner.i, owner.j, owner.k, m);
    }
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
