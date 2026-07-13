#include "4_halo/Halo.h"
#include "0_basic/Error.h"
#include "0_basic/LayoutTraits.h"

#include <algorithm>
#include <cstdlib>
#include <ostream>
#include <sstream>
#include <vector>

namespace
{
    const char *halo_level_name(HaloLevel level)
    {
        switch (level)
        {
        case HaloLevel::Corner1D:
            return "1DCorner";
        case HaloLevel::Corner2D:
            return "2DCorner";
        case HaloLevel::Corner3D:
            return "3DCorner";
        }

        return "Unknown";
    }

    int entity_axis_index(TOPO::EntityAxis axis)
    {
        switch (axis)
        {
        case TOPO::EntityAxis::Xi:
            return 0;
        case TOPO::EntityAxis::Eta:
            return 1;
        case TOPO::EntityAxis::Zeta:
            return 2;
        case TOPO::EntityAxis::None:
            break;
        }

        ERROR::Abort("[Halo] owner sync requires oriented edge/face axis.");
        return 0;
    }
}

const FieldHaloRequest &Halo::halo_request_(const std::string &field_name) const
{
    auto it = halo_registry_.find(field_name);
    if (it == halo_registry_.end())
        ERROR::Abort("[Halo] halo_request_: field is not registered: " + field_name);

    return it->second;
}

bool Halo::halo_level_includes_edge_(HaloLevel level) const
{
    return level == HaloLevel::Corner2D ||
           level == HaloLevel::Corner3D;
}

bool Halo::halo_level_includes_vertex_(HaloLevel level) const
{
    return level == HaloLevel::Corner3D;
}

void Halo::sync_component_copy_field_(const std::string &field_name)
{
    const FieldHaloRequest &req = halo_request_(field_name);

    sync_component_copy_field_stage_(field_name, HaloLevel::Corner1D);

    if (halo_level_includes_edge_(req.level))
        sync_component_copy_field_stage_(field_name, HaloLevel::Corner2D);

    if (halo_level_includes_vertex_(req.level))
        sync_component_copy_field_stage_(field_name, HaloLevel::Corner3D);
}

void Halo::sync_component_copy_field_stage_(const std::string &field_name, HaloLevel stage)
{
    const FieldHaloRequest &req = halo_request_(field_name);
    if (static_cast<int>(req.level) < static_cast<int>(stage))
        return;

    if (stage == HaloLevel::Corner1D)
    {
        exchange_inner(field_name);
        exchange_parallel(field_name);
    }
    else if (stage == HaloLevel::Corner2D)
    {
        exchange_inner_edge(field_name);
        exchange_parallel_edge(field_name);
    }
    else if (stage == HaloLevel::Corner3D)
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

void Halo::sync_component_copy_registered_stage_(HaloLevel stage)
{
    for (const std::string &field_name : component_copy_fields_)
        sync_component_copy_field_stage_(field_name, stage);
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

void Halo::sync_edge_1form_triplet_stage_(const HaloTripletRequest &tri, HaloLevel stage)
{
    if (static_cast<int>(tri.level) < static_cast<int>(stage))
        return;

    const std::vector<std::string> fields{tri.xi, tri.eta, tri.zeta};
    if (stage == HaloLevel::Corner1D)
    {
        exchange_inner_face_edge_1form_triplet_(fields);
        exchange_parallel_face_edge_1form_triplet_(fields);
    }
    else if (stage == HaloLevel::Corner2D)
    {
        exchange_inner_edge_edge_1form_triplet_(fields);
        exchange_parallel_edge_edge_1form_triplet_(fields);
    }
    else if (stage == HaloLevel::Corner3D)
    {
        exchange_inner_vertex_edge_1form_triplet_(fields);
        exchange_parallel_vertex_edge_1form_triplet_(fields);
    }
}

void Halo::sync_edge_1form_triplets_registered_()
{
    for (const auto &kv : edge_1form_triplets_)
        sync_edge_1form_triplet_(kv.second);
}

void Halo::sync_edge_1form_triplets_registered_stage_(HaloLevel stage)
{
    for (const auto &kv : edge_1form_triplets_)
        sync_edge_1form_triplet_stage_(kv.second, stage);
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

    if (halo_level_includes_edge_(tri.level))
        sync_face_2form_triplet_edge_level_(tri);

    if (halo_level_includes_vertex_(tri.level))
        sync_face_2form_triplet_vertex_level_(tri);
}

void Halo::sync_face_2form_triplet_stage_(const HaloTripletRequest &tri, HaloLevel stage)
{
    if (static_cast<int>(tri.level) < static_cast<int>(stage))
        return;

    if (stage == HaloLevel::Corner1D)
        sync_face_2form_triplet_face_level_(tri);
    else if (stage == HaloLevel::Corner2D)
        sync_face_2form_triplet_edge_level_(tri);
    else if (stage == HaloLevel::Corner3D)
        sync_face_2form_triplet_vertex_level_(tri);
}

void Halo::sync_face_2form_triplets_registered_()
{
    for (const auto &kv : face_2form_triplets_)
        sync_face_2form_triplet_(kv.second);
}

void Halo::sync_face_2form_triplets_registered_stage_(HaloLevel stage)
{
    for (const auto &kv : face_2form_triplets_)
        sync_face_2form_triplet_stage_(kv.second, stage);
}

void Halo::sync_owner_alias_request_(const HaloOwnerRequest &req)
{
    if (req.policy == OwnerSyncPolicy::None)
        return;

    auto it = owner_sync_patterns_.find(req.field_name);
    if (it == owner_sync_patterns_.end())
        ERROR::Abort("[Halo] owner sync pattern not built for field: " + req.field_name);

    OwnerSyncPattern &pat = it->second;

    execute_owner_sync_local_ops_(pat);
    execute_owner_sync_mpi_ops_(pat);
}

void Halo::sync_owner_alias_registered_()
{
    for (const auto &req : owner_sync_requests_)
        sync_owner_alias_request_(req);
}

bool Halo::sync_owner_alias_field_(const std::string &field_name)
{
    bool handled = false;
    for (const auto &own : owner_sync_requests_)
    {
        if (own.field_name == field_name)
        {
            sync_owner_alias_request_(own);
            handled = true;
        }
    }
    return handled;
}

bool Halo::sync_owner_alias_group_(const std::string &group_name)
{
    bool handled = false;
    for (const auto &own : owner_sync_requests_)
    {
        if (own.sync_group == group_name)
        {
            sync_owner_alias_request_(own);
            handled = true;
        }
    }
    return handled;
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

void Halo::dump_sync_registry(std::ostream &os) const
{
    os << "========== Halo Sync Registry ==========\n";

    os << "ComponentCopy fields (" << component_copy_fields_.size() << "):\n";
    for (const auto &name : component_copy_fields_)
    {
        const auto &req = halo_request_(name);
        const std::string group = req.sync_group.empty() ? req.field_name : req.sync_group;

        os << "  - " << name
           << " group=" << group
           << " loc=" << LAYOUT::location_name(req.location)
           << " kind=" << field_value_kind_name(req.value_kind)
           << " do_halo=" << (req.do_halo ? "true" : "false")
           << " level=" << halo_level_name(req.level)
           << " owner=" << owner_sync_policy_name(req.owner_sync)
           << " orientation=" << (req.orientation_aware ? "true" : "false")
           << "\n";
    }

    os << "Edge1FormTriplets (" << edge_1form_triplets_.size() << "):\n";
    for (const auto &kv : edge_1form_triplets_)
    {
        const auto &tri = kv.second;
        os << "  - " << kv.first
           << " xi=" << tri.xi
           << " eta=" << tri.eta
           << " zeta=" << tri.zeta
           << " level=" << halo_level_name(tri.level)
           << " nghost=" << tri.nghost
           << "\n";
    }

    os << "Face2FormTriplets (" << face_2form_triplets_.size() << "):\n";
    for (const auto &kv : face_2form_triplets_)
    {
        const auto &tri = kv.second;
        os << "  - " << kv.first
           << " xi=" << tri.xi
           << " eta=" << tri.eta
           << " zeta=" << tri.zeta
           << " level=" << halo_level_name(tri.level)
           << " nghost=" << tri.nghost
           << "\n";
    }

    os << "OwnerAlias requests (" << owner_sync_requests_.size() << "):\n";
    for (const auto &own : owner_sync_requests_)
    {
        os << "  - field=" << own.field_name
           << " group=" << own.sync_group
           << " policy=" << owner_sync_policy_name(own.policy)
           << " loc=" << LAYOUT::location_name(own.location)
           << " kind=" << field_value_kind_name(own.value_kind)
           << " orientation=" << (own.orientation_aware ? "true" : "false")
           << "\n";
    }

    os << "OwnerSyncPatterns (" << owner_sync_patterns_.size() << "):\n";
    for (const auto &kv : owner_sync_patterns_)
    {
        const auto &pat = kv.second;
        os << "  - field=" << pat.field_name
           << " group=" << pat.sync_group
           << " local_ops=" << pat.local_ops.size()
           << " send_ops=" << pat.send_ops.size()
           << " recv_ops=" << pat.recv_ops.size()
           << "\n";
    }

    os << "========================================\n";
}

TOPO::EntityDim Halo::owner_policy_to_entity_dim_(OwnerSyncPolicy policy) const
{
    switch (policy)
    {
    case OwnerSyncPolicy::NodeOwner:
        return TOPO::EntityDim::Node;
    case OwnerSyncPolicy::EdgeOwner:
        return TOPO::EntityDim::Edge;
    case OwnerSyncPolicy::FaceOwner:
        return TOPO::EntityDim::Face;
    case OwnerSyncPolicy::None:
        break;
    }

    ERROR::Abort("[Halo] owner_policy_to_entity_dim_: invalid OwnerSyncPolicy");
    return TOPO::EntityDim::Node;
}

void Halo::require_owner_equiv_available_(OwnerSyncPolicy policy,
                                          const std::string &field_name) const
{
    if (policy == OwnerSyncPolicy::None)
        return;

    if (policy == OwnerSyncPolicy::FaceOwner &&
        (!equiv_ || !equiv_->has_face_equiv()))
    {
        ERROR::Abort("[Halo] FaceOwner sync requested but face equivalence is not built. field=" +
                     field_name);
    }

    if (!equiv_)
    {
        ERROR::Abort("[Halo] OwnerAliasSync requested but Topology is not bound. field=" +
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
        // orient_sign is relative to the canonical face orientation.
        // alias_from_owner_sign = alias_sign / owner_sign. Since signs are
        // in {+1, -1}, division is equivalent to alias_sign * owner_sign.
        return alias.orient_sign * owner.orient_sign;
    }

    return +1;
}

bool Halo::owner_member_matches_field_(const HaloOwnerRequest &req,
                                       const TOPO::EquivMember &member) const
{
    if (req.policy == OwnerSyncPolicy::FaceOwner)
    {
        switch (req.location)
        {
        case StaggerLocation::FaceXi:
            return member.entity.axis == TOPO::EntityAxis::Xi;
        case StaggerLocation::FaceEt:
            return member.entity.axis == TOPO::EntityAxis::Eta;
        case StaggerLocation::FaceZe:
            return member.entity.axis == TOPO::EntityAxis::Zeta;
        default:
            return false;
        }
    }

    return TOPO::stagger_location(member.entity) == req.location;
}

void Halo::copy_owner_to_alias_local_(const HaloOwnerRequest &req,
                                      int fid,
                                      const TOPO::EquivMember &owner,
                                      const TOPO::EquivMember &alias,
                                      int sign)
{
    FieldBlock &owner_block = fld_->field(fid, owner.entity.block);
    FieldBlock &alias_block = fld_->field(fid, alias.entity.block);

    if (!owner_block.is_allocated() || !alias_block.is_allocated())
        return;

    for (int m = 0; m < req.ncomp; ++m)
    {
        alias_block(alias.entity.i, alias.entity.j, alias.entity.k, m) =
            static_cast<double>(sign) * owner_block(owner.entity.i, owner.entity.j, owner.entity.k, m);
    }
}

int Halo::owner_sync_member_field_id_(const HaloOwnerRequest &req,
                                      const TOPO::EquivMember &member) const
{
    if (req.value_kind == FieldValueKind::EdgeCovariant1Form ||
        req.value_kind == FieldValueKind::FaceContravariant2Form)
    {
        const int axis = entity_axis_index(member.entity.axis);
        return fld_->field_id(
            find_triplet_field_name_(req.field_name, req.value_kind, axis));
    }

    return fld_->field_id(req.field_name);
}

void Halo::build_owner_sync_patterns_()
{
    owner_sync_patterns_.clear();

    for (const auto &req : owner_sync_requests_)
    {
        OwnerSyncPattern pat = build_owner_sync_pattern_for_request_(req);
        owner_sync_patterns_[req.field_name] = pat;
    }
}

Halo::OwnerSyncPattern Halo::build_owner_sync_pattern_for_request_(const HaloOwnerRequest &req) const
{
    require_owner_equiv_available_(req.policy, req.field_name);

    OwnerSyncPattern pat;
    pat.field_name = req.field_name;
    pat.sync_group = req.sync_group;
    pat.policy = req.policy;
    pat.value_kind = req.value_kind;
    pat.location = req.location;
    pat.orientation_aware = req.orientation_aware;

    const int alias_fid = fld_->field_id(req.field_name);
    const TOPO::EntityDim dim = owner_policy_to_entity_dim_(req.policy);
    const auto &classes = equiv_->classes(dim);

    int my_rank = 0;
    PARALLEL::mpi_rank(&my_rank);

    for (const auto &cls : classes)
    {
        const TOPO::EquivMember &owner = cls.owner;
        const int owner_fid = owner_sync_member_field_id_(req, owner);

        for (const auto &alias : cls.members)
        {
            if (alias.is_owner)
                continue;

            if (!owner_member_matches_field_(req, alias))
                continue;

            const int sign = owner_alias_sign_(req, owner, alias);

            if (owner.entity.rank == my_rank && alias.entity.rank == my_rank)
            {
                OwnerSyncLocalOp op;
                op.owner_fid = owner_fid;
                op.alias_fid = alias_fid;
                op.owner_block = owner.entity.block;
                op.owner_i = owner.entity.i;
                op.owner_j = owner.entity.j;
                op.owner_k = owner.entity.k;
                op.alias_block = alias.entity.block;
                op.alias_i = alias.entity.i;
                op.alias_j = alias.entity.j;
                op.alias_k = alias.entity.k;
                op.ncomp = req.ncomp;
                op.sign = sign;
                pat.local_ops.push_back(op);
            }
            else if (owner.entity.rank == my_rank && alias.entity.rank != my_rank)
            {
                if (cls.global_id < 0)
                {
                    if (req.policy == OwnerSyncPolicy::FaceOwner)
                        ERROR::Abort("[Halo] cross-rank FaceOwner sync requires valid face global_id");

                    ERROR::Abort("[Halo] cross-rank OwnerAliasSync requires EquivClass::global_id >= 0 for field: " +
                                 req.field_name);
                }

                OwnerSyncSendOp op;
                op.owner_fid = owner_fid;
                op.class_gid = cls.global_id;
                op.owner_block = owner.entity.block;
                op.owner_i = owner.entity.i;
                op.owner_j = owner.entity.j;
                op.owner_k = owner.entity.k;
                op.alias_rank = alias.entity.rank;
                op.alias_block = alias.entity.block;
                op.alias_i = alias.entity.i;
                op.alias_j = alias.entity.j;
                op.alias_k = alias.entity.k;
                op.ncomp = req.ncomp;
                op.sign_for_alias = sign;
                op.dst_rank = alias.entity.rank;
                op.tag = owner_sync_tag_(req);
                pat.send_ops.push_back(op);
            }
            else if (owner.entity.rank != my_rank && alias.entity.rank == my_rank)
            {
                if (cls.global_id < 0)
                {
                    if (req.policy == OwnerSyncPolicy::FaceOwner)
                        ERROR::Abort("[Halo] cross-rank FaceOwner sync requires valid face global_id");

                    ERROR::Abort("[Halo] cross-rank OwnerAliasSync requires EquivClass::global_id >= 0 for field: " +
                                 req.field_name);
                }

                OwnerSyncRecvOp op;
                op.alias_fid = alias_fid;
                op.class_gid = cls.global_id;
                op.alias_block = alias.entity.block;
                op.alias_i = alias.entity.i;
                op.alias_j = alias.entity.j;
                op.alias_k = alias.entity.k;
                op.ncomp = req.ncomp;
                op.src_rank = owner.entity.rank;
                op.tag = owner_sync_tag_(req);
                pat.recv_ops.push_back(op);
            }
        }
    }

    std::sort(pat.send_ops.begin(), pat.send_ops.end(),
              [](const OwnerSyncSendOp &a, const OwnerSyncSendOp &b)
              {
                  if (a.dst_rank != b.dst_rank)
                      return a.dst_rank < b.dst_rank;
                  if (a.class_gid != b.class_gid)
                      return a.class_gid < b.class_gid;
                  if (a.alias_rank != b.alias_rank)
                      return a.alias_rank < b.alias_rank;
                  if (a.alias_block != b.alias_block)
                      return a.alias_block < b.alias_block;
                  if (a.alias_i != b.alias_i)
                      return a.alias_i < b.alias_i;
                  if (a.alias_j != b.alias_j)
                      return a.alias_j < b.alias_j;
                  if (a.alias_k != b.alias_k)
                      return a.alias_k < b.alias_k;
                  return a.owner_fid < b.owner_fid;
              });

    std::sort(pat.recv_ops.begin(), pat.recv_ops.end(),
              [](const OwnerSyncRecvOp &a, const OwnerSyncRecvOp &b)
              {
                  if (a.src_rank != b.src_rank)
                      return a.src_rank < b.src_rank;
                  if (a.class_gid != b.class_gid)
                      return a.class_gid < b.class_gid;
                  if (a.alias_block != b.alias_block)
                      return a.alias_block < b.alias_block;
                  if (a.alias_i != b.alias_i)
                      return a.alias_i < b.alias_i;
                  if (a.alias_j != b.alias_j)
                      return a.alias_j < b.alias_j;
                  if (a.alias_k != b.alias_k)
                      return a.alias_k < b.alias_k;
                  return a.alias_fid < b.alias_fid;
              });

    resize_owner_sync_buffers_(pat);
    if (std::getenv("Z0_OWNER_DEBUG") != nullptr)
    {
        int my_rank = 0;
        PARALLEL::mpi_rank(&my_rank);
        std::cout << "[OwnerPattern] rank=" << my_rank
                  << " field=" << pat.field_name
                  << " policy=" << owner_sync_policy_name(pat.policy)
                  << " local=" << pat.local_ops.size()
                  << " send=" << pat.send_ops.size()
                  << " recv=" << pat.recv_ops.size()
                  << "\n";
    }
    return pat;
}

int Halo::owner_sync_tag_(const HaloOwnerRequest &req) const
{
    const int fid = fld_ ? fld_->field_id(req.field_name) : 0;
    if (req.policy == OwnerSyncPolicy::FaceOwner)
        return owner_sync_tag_base_ + 100 + fid;
    if (req.policy == OwnerSyncPolicy::EdgeOwner)
        return owner_sync_tag_base_ + 200 + fid;
    if (req.policy == OwnerSyncPolicy::NodeOwner)
        return owner_sync_tag_base_ + 300 + fid;

    ERROR::Abort("[Halo] owner_sync_tag_: invalid OwnerSyncPolicy for field: " + req.field_name);
    return owner_sync_tag_base_;
}

void Halo::resize_owner_sync_buffers_(OwnerSyncPattern &pat) const
{
    int send_size = 0;
    for (auto &op : pat.send_ops)
    {
        op.buffer_offset = send_size;
        send_size += op.ncomp;
    }

    int recv_size = 0;
    for (auto &op : pat.recv_ops)
    {
        op.buffer_offset = recv_size;
        recv_size += op.ncomp;
    }

    pat.send_buffer.resize(send_size);
    pat.recv_buffer.resize(recv_size);
}

void Halo::execute_owner_sync_local_ops_(const OwnerSyncPattern &pat)
{
    for (const auto &op : pat.local_ops)
    {
        FieldBlock &owner_block = fld_->field(op.owner_fid, op.owner_block);
        FieldBlock &alias_block = fld_->field(op.alias_fid, op.alias_block);

        if (!owner_block.is_allocated() || !alias_block.is_allocated())
            continue;

        for (int m = 0; m < op.ncomp; ++m)
        {
            alias_block(op.alias_i, op.alias_j, op.alias_k, m) =
                static_cast<double>(op.sign) *
                owner_block(op.owner_i, op.owner_j, op.owner_k, m);
        }
    }
}

void Halo::pack_owner_sync_send_buffer_(OwnerSyncPattern &pat)
{
    for (const auto &op : pat.send_ops)
    {
        FieldBlock &owner_block = fld_->field(op.owner_fid, op.owner_block);

        if (!owner_block.is_allocated())
            ERROR::Abort("[Halo] owner sync send op references inactive owner field: " + pat.field_name);

        for (int m = 0; m < op.ncomp; ++m)
        {
            pat.send_buffer[op.buffer_offset + m] =
                static_cast<double>(op.sign_for_alias) *
                owner_block(op.owner_i, op.owner_j, op.owner_k, m);
        }
    }
}

void Halo::unpack_owner_sync_recv_buffer_(OwnerSyncPattern &pat)
{
    for (const auto &op : pat.recv_ops)
    {
        FieldBlock &alias_block = fld_->field(op.alias_fid, op.alias_block);

        if (!alias_block.is_allocated())
            continue;

        for (int m = 0; m < op.ncomp; ++m)
        {
            alias_block(op.alias_i, op.alias_j, op.alias_k, m) =
                pat.recv_buffer[op.buffer_offset + m];
        }
    }
}

void Halo::execute_owner_sync_mpi_ops_(OwnerSyncPattern &pat)
{
    pack_owner_sync_send_buffer_(pat);

    if(!pat.mpi_lengths_validated)
    {
        int nrank = 1;
        PARALLEL::mpi_size(&nrank);
        std::vector<int> send_lengths_by_rank(nrank, 0);
        std::vector<int> recv_lengths_by_rank(nrank, 0);
        std::vector<int> expected_recv_lengths_by_rank(nrank, 0);

        for (const auto &op : pat.send_ops)
            send_lengths_by_rank[op.dst_rank] += op.ncomp;
        for (const auto &op : pat.recv_ops)
            recv_lengths_by_rank[op.src_rank] += op.ncomp;

        PARALLEL::mpi_alltoall(send_lengths_by_rank.data(), 1,
                               expected_recv_lengths_by_rank.data(), 1);

        for (int r = 0; r < nrank; ++r)
        {
            if (recv_lengths_by_rank[r] != expected_recv_lengths_by_rank[r])
            {
                int my_rank = 0;
                PARALLEL::mpi_rank(&my_rank);
                std::ostringstream oss;
                oss << "[Halo] owner sync length mismatch field=" << pat.field_name
                    << " policy=" << owner_sync_policy_name(pat.policy)
                    << " rank=" << my_rank
                    << " peer_rank=" << r
                    << " posted_recv_len=" << recv_lengths_by_rank[r]
                    << " expected_from_peer_send_len=" << expected_recv_lengths_by_rank[r];
                ERROR::Abort(oss.str());
            }
        }
        pat.mpi_lengths_validated=true;
    }

    if (pat.send_ops.empty() && pat.recv_ops.empty())
        return;

    auto &recv_requests=owner_recv_requests_;
    auto &send_requests=owner_send_requests_;
    auto &recv_statuses=owner_recv_statuses_;
    auto &send_statuses=owner_send_statuses_;
    recv_requests.clear();
    send_requests.clear();

    recv_requests.reserve(pat.recv_ops.size());
    send_requests.reserve(pat.send_ops.size());

    for (std::size_t first = 0; first < pat.recv_ops.size();)
    {
        std::size_t last = first + 1;
        while (last < pat.recv_ops.size() &&
               pat.recv_ops[last].src_rank == pat.recv_ops[first].src_rank)
        {
            ++last;
        }

        int length = 0;
        for (std::size_t i = first; i < last; ++i)
            length += pat.recv_ops[i].ncomp;

        const OwnerSyncRecvOp &op = pat.recv_ops[first];
        recv_requests.push_back(MPI_Request{});
        PARALLEL::mpi_data_recv(op.src_rank,
                                op.tag,
                                pat.recv_buffer.data() + op.buffer_offset,
                                length,
                                &recv_requests.back());

        first = last;
    }

    for (std::size_t first = 0; first < pat.send_ops.size();)
    {
        std::size_t last = first + 1;
        while (last < pat.send_ops.size() &&
               pat.send_ops[last].dst_rank == pat.send_ops[first].dst_rank)
        {
            ++last;
        }

        int length = 0;
        for (std::size_t i = first; i < last; ++i)
            length += pat.send_ops[i].ncomp;

        const OwnerSyncSendOp &op = pat.send_ops[first];
        send_requests.push_back(MPI_Request{});
        PARALLEL::mpi_data_send(op.dst_rank,
                                op.tag,
                                pat.send_buffer.data() + op.buffer_offset,
                                length,
                                &send_requests.back());

        first = last;
    }

    int nrecv = static_cast<int>(recv_requests.size());
    int nsend = static_cast<int>(send_requests.size());
    recv_statuses.resize(nrecv);
    send_statuses.resize(nsend);

    if (nrecv > 0)
        PARALLEL::mpi_wait(nrecv, recv_requests.data(), recv_statuses.data());
    if (nsend > 0)
        PARALLEL::mpi_wait(nsend, send_requests.data(), send_statuses.data());

    unpack_owner_sync_recv_buffer_(pat);
}

void Halo::sync_registered()
{
    sync_component_copy_registered_();

    sync_edge_1form_triplets_registered_();

    sync_face_2form_triplets_registered_();

    sync_owner_alias_registered_();
}

void Halo::sync_registered(HaloLevel stage)
{
    sync_component_copy_registered_stage_(stage);

    sync_edge_1form_triplets_registered_stage_(stage);

    sync_face_2form_triplets_registered_stage_(stage);
}

void Halo::sync_owner_alias()
{
    sync_owner_alias_registered_();
}

void Halo::sync_owner_alias_field(const std::string &field_name)
{
    if (!sync_owner_alias_field_(field_name))
        ERROR::Abort("[Halo] sync_owner_alias_field: field has no owner sync request: " + field_name);
}

void Halo::sync_owner_alias_group(const std::string &group_name)
{
    if (!sync_owner_alias_group_(group_name))
        ERROR::Abort("[Halo] sync_owner_alias_group: group has no owner sync requests: " + group_name);
}

void Halo::sync_field(const std::string &field_name)
{
    const FieldHaloRequest &req = halo_request_(field_name);

    const std::string group =
        req.sync_group.empty() ? req.field_name : req.sync_group;

    if (req.do_halo)
    {
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
    }

    sync_owner_alias_field_(field_name);
}

void Halo::sync_field(const std::string &field_name, HaloLevel stage)
{
    const FieldHaloRequest &req = halo_request_(field_name);

    const std::string group =
        req.sync_group.empty() ? req.field_name : req.sync_group;

    if (req.do_halo)
    {
        const HaloSyncSemantics sem = sync_semantics_(req);

        if (sem == HaloSyncSemantics::ComponentCopy)
        {
            sync_component_copy_field_stage_(field_name, stage);
        }
        else if (sem == HaloSyncSemantics::Edge1FormTriplet)
        {
            auto it = edge_1form_triplets_.find(group);
            if (it == edge_1form_triplets_.end())
                ERROR::Abort("[Halo] sync_field(stage): missing edge 1-form group for field: " + field_name);

            sync_edge_1form_triplet_stage_(it->second, stage);
        }
        else if (sem == HaloSyncSemantics::Face2FormTriplet)
        {
            auto it = face_2form_triplets_.find(group);
            if (it == face_2form_triplets_.end())
                ERROR::Abort("[Halo] sync_field(stage): missing face 2-form group for field: " + field_name);

            sync_face_2form_triplet_stage_(it->second, stage);
        }
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

    handled = sync_owner_alias_group_(group_name) || handled;

    if (!handled)
        ERROR::Abort("[Halo] sync_group: unknown or empty group: " + group_name);
}

void Halo::sync_group(const std::string &group_name, HaloLevel stage)
{
    bool handled = false;

    auto eit = edge_1form_triplets_.find(group_name);
    if (eit != edge_1form_triplets_.end())
    {
        sync_edge_1form_triplet_stage_(eit->second, stage);
        handled = true;
    }

    auto fit = face_2form_triplets_.find(group_name);
    if (fit != face_2form_triplets_.end())
    {
        sync_face_2form_triplet_stage_(fit->second, stage);
        handled = true;
    }

    for (const auto &field_name : component_copy_fields_)
    {
        auto git = field_to_group_.find(field_name);

        if (git != field_to_group_.end() && git->second == group_name)
        {
            sync_component_copy_field_stage_(field_name, stage);
            handled = true;
        }
    }

    if (!handled)
        ERROR::Abort("[Halo] sync_group(stage): unknown or empty group: " + group_name);
}
