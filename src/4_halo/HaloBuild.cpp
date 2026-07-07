#include "4_halo/Halo.h"
#include "2_topology/TopologyView.h"
#include "0_basic/Error.h"

#include <algorithm>
#include <set>

namespace
{
    int triplet_axis_from_location(StaggerLocation loc)
    {
        switch (loc)
        {
        case StaggerLocation::EdgeXi:
        case StaggerLocation::FaceXi:
            return 0;

        case StaggerLocation::EdgeEt:
        case StaggerLocation::FaceEt:
            return 1;

        case StaggerLocation::EdgeZe:
        case StaggerLocation::FaceZe:
            return 2;

        default:
            return -1;
        }
    }

    bool same_halo_request_metadata(const FieldHaloRequest &a, const FieldHaloRequest &b)
    {
        return a.sync_group == b.sync_group &&
               a.location == b.location &&
               a.value_kind == b.value_kind &&
               a.ncomp == b.ncomp &&
               a.nghost == b.nghost &&
               a.owner_sync == b.owner_sync &&
               a.orientation_aware == b.orientation_aware;
    }

    bool is_edge_location(StaggerLocation loc)
    {
        return loc == StaggerLocation::EdgeXi ||
               loc == StaggerLocation::EdgeEt ||
               loc == StaggerLocation::EdgeZe;
    }

    bool is_face_location(StaggerLocation loc)
    {
        return loc == StaggerLocation::FaceXi ||
               loc == StaggerLocation::FaceEt ||
               loc == StaggerLocation::FaceZe;
    }

    bool location_matches_triplet_axis(StaggerLocation loc,
                                       int axis,
                                       bool face)
    {
        if (face)
        {
            if (axis == 0)
                return loc == StaggerLocation::FaceXi;
            if (axis == 1)
                return loc == StaggerLocation::FaceEt;
            if (axis == 2)
                return loc == StaggerLocation::FaceZe;
        }
        else
        {
            if (axis == 0)
                return loc == StaggerLocation::EdgeXi;
            if (axis == 1)
                return loc == StaggerLocation::EdgeEt;
            if (axis == 2)
                return loc == StaggerLocation::EdgeZe;
        }

        return false;
    }
}

Halo::HaloSyncSemantics Halo::sync_semantics_(const FieldHaloRequest &req) const
{
    if (!req.orientation_aware)
        return HaloSyncSemantics::ComponentCopy;

    if (req.value_kind == FieldValueKind::EdgeCovariant1Form)
        return HaloSyncSemantics::Edge1FormTriplet;

    if (req.value_kind == FieldValueKind::FaceContravariant2Form)
        return HaloSyncSemantics::Face2FormTriplet;

    return HaloSyncSemantics::ComponentCopy;
}

void Halo::classify_registered_request_(const FieldHaloRequest &req)
{
    const std::string group =
        req.sync_group.empty() ? req.field_name : req.sync_group;

    field_to_group_[req.field_name] = group;

    const HaloSyncSemantics sem = sync_semantics_(req);

    if (sem == HaloSyncSemantics::ComponentCopy)
    {
        component_copy_fields_.push_back(req.field_name);
    }
    else if (sem == HaloSyncSemantics::Edge1FormTriplet)
    {
        HaloTripletRequest &tri = edge_1form_triplets_[group];

        if (tri.group_name.empty())
            tri.group_name = group;

        tri.value_kind = req.value_kind;
        tri.level = req.level;
        tri.nghost = req.nghost;
        tri.orientation_aware = req.orientation_aware;

        const int ax = triplet_axis_from_location(req.location);
        if (ax == 0)
            tri.xi = req.field_name;
        else if (ax == 1)
            tri.eta = req.field_name;
        else if (ax == 2)
            tri.zeta = req.field_name;
        else
            ERROR::Abort("[Halo] invalid triplet component location for field: " + req.field_name);
    }
    else if (sem == HaloSyncSemantics::Face2FormTriplet)
    {
        HaloTripletRequest &tri = face_2form_triplets_[group];

        if (tri.group_name.empty())
            tri.group_name = group;

        tri.value_kind = req.value_kind;
        tri.level = req.level;
        tri.nghost = req.nghost;
        tri.orientation_aware = req.orientation_aware;

        const int ax = triplet_axis_from_location(req.location);
        if (ax == 0)
            tri.xi = req.field_name;
        else if (ax == 1)
            tri.eta = req.field_name;
        else if (ax == 2)
            tri.zeta = req.field_name;
        else
            ERROR::Abort("[Halo] invalid triplet component location for field: " + req.field_name);
    }

    if (req.owner_sync != OwnerSyncPolicy::None)
    {
        HaloOwnerRequest own;

        own.field_name = req.field_name;
        own.sync_group = group;
        own.policy = req.owner_sync;
        own.value_kind = req.value_kind;
        own.location = req.location;
        own.ncomp = req.ncomp;
        own.orientation_aware = req.orientation_aware;

        owner_sync_requests_.push_back(own);
    }
}

void Halo::rebuild_sync_registry_()
{
    component_copy_fields_.clear();
    edge_1form_triplets_.clear();
    face_2form_triplets_.clear();
    owner_sync_requests_.clear();
    field_to_group_.clear();
    owner_sync_patterns_.clear();

    for (const auto &kv : halo_registry_)
        classify_registered_request_(kv.second);

    std::sort(owner_sync_requests_.begin(), owner_sync_requests_.end(),
              [](const HaloOwnerRequest &a, const HaloOwnerRequest &b)
              {
                  return a.field_name < b.field_name;
              });
    std::sort(component_copy_fields_.begin(), component_copy_fields_.end());

    validate_triplet_registry_();
    validate_sync_registry_consistency_();
}

void Halo::validate_triplet_registry_() const
{
    for (const auto &kv : edge_1form_triplets_)
    {
        const HaloTripletRequest &tri = kv.second;

        if (tri.xi.empty() || tri.eta.empty() || tri.zeta.empty())
            ERROR::Abort("[Halo] incomplete edge 1-form triplet group: " + kv.first);

        if (tri.value_kind != FieldValueKind::EdgeCovariant1Form)
            ERROR::Abort("[Halo] edge 1-form triplet group is not EdgeCovariant1Form: " + kv.first);

        if (!tri.orientation_aware)
            ERROR::Abort("[Halo] edge 1-form triplet group is not orientation-aware: " + kv.first);
    }

    for (const auto &kv : face_2form_triplets_)
    {
        const HaloTripletRequest &tri = kv.second;

        if (tri.xi.empty() || tri.eta.empty() || tri.zeta.empty())
            ERROR::Abort("[Halo] incomplete face 2-form triplet group: " + kv.first);
    }
}

void Halo::validate_sync_registry_consistency_() const
{
    for (const auto &field_name : component_copy_fields_)
    {
        const FieldHaloRequest &req = halo_request_(field_name);
        if (sync_semantics_(req) != HaloSyncSemantics::ComponentCopy)
        {
            ERROR::Abort("[Halo] sync registry inconsistency: non-component field in component_copy_fields_: " +
                         field_name);
        }
    }

    for (const auto &kv : edge_1form_triplets_)
    {
        const HaloTripletRequest &tri = kv.second;
        const std::string fields[3] = {tri.xi, tri.eta, tri.zeta};

        for (int axis = 0; axis < 3; ++axis)
        {
            const std::string &field_name = fields[axis];
            auto it = halo_registry_.find(field_name);
            if (it == halo_registry_.end())
                ERROR::Abort("[Halo] edge 1-form triplet field is not registered: " + field_name);

            const FieldHaloRequest &req = it->second;
            const std::string group = req.sync_group.empty() ? req.field_name : req.sync_group;

            if (req.value_kind != FieldValueKind::EdgeCovariant1Form)
                ERROR::Abort("[Halo] edge 1-form triplet field has wrong value_kind: " + field_name);

            if (!req.orientation_aware)
                ERROR::Abort("[Halo] edge 1-form triplet field is not orientation-aware: " + field_name);

            if (group != kv.first)
                ERROR::Abort("[Halo] edge 1-form triplet field has wrong sync_group: " + field_name);

            if (req.nghost != tri.nghost)
                ERROR::Abort("[Halo] edge 1-form triplet field has inconsistent nghost: " + field_name);

            if (req.level != tri.level)
                ERROR::Abort("[Halo] edge 1-form triplet field has inconsistent halo level: " + field_name);

            if (!location_matches_triplet_axis(req.location, axis, false))
                ERROR::Abort("[Halo] edge 1-form triplet field has wrong location: " + field_name);
        }
    }

    for (const auto &kv : face_2form_triplets_)
    {
        const HaloTripletRequest &tri = kv.second;
        const std::string fields[3] = {tri.xi, tri.eta, tri.zeta};

        for (int axis = 0; axis < 3; ++axis)
        {
            const std::string &field_name = fields[axis];
            auto it = halo_registry_.find(field_name);
            if (it == halo_registry_.end())
                ERROR::Abort("[Halo] face 2-form triplet field is not registered: " + field_name);

            const FieldHaloRequest &req = it->second;
            const std::string group = req.sync_group.empty() ? req.field_name : req.sync_group;

            if (req.value_kind != FieldValueKind::FaceContravariant2Form)
                ERROR::Abort("[Halo] face 2-form triplet field has wrong value_kind: " + field_name);

            if (!req.orientation_aware)
                ERROR::Abort("[Halo] face 2-form triplet field is not orientation-aware: " + field_name);

            if (group != kv.first)
                ERROR::Abort("[Halo] face 2-form triplet field has wrong sync_group: " + field_name);

            if (req.nghost != tri.nghost)
                ERROR::Abort("[Halo] face 2-form triplet field has inconsistent nghost: " + field_name);

            if (req.level != tri.level)
                ERROR::Abort("[Halo] face 2-form triplet field has inconsistent halo level: " + field_name);

            if (!location_matches_triplet_axis(req.location, axis, true))
                ERROR::Abort("[Halo] face 2-form triplet field has wrong location: " + field_name);
        }
    }

    for (const auto &own : owner_sync_requests_)
    {
        if (halo_registry_.find(own.field_name) == halo_registry_.end())
            ERROR::Abort("[Halo] owner sync field is not registered: " + own.field_name);

        if (own.policy == OwnerSyncPolicy::None)
            ERROR::Abort("[Halo] owner sync request has OwnerSyncPolicy::None: " + own.field_name);

        if (own.policy == OwnerSyncPolicy::NodeOwner &&
            own.location != StaggerLocation::Node)
        {
            ERROR::Abort("[Halo] NodeOwner sync request must use Node location: " + own.field_name);
        }

        if (own.policy == OwnerSyncPolicy::EdgeOwner &&
            !is_edge_location(own.location))
        {
            ERROR::Abort("[Halo] EdgeOwner sync request must use edge location: " + own.field_name);
        }

        if (own.policy == OwnerSyncPolicy::FaceOwner &&
            !is_face_location(own.location))
        {
            ERROR::Abort("[Halo] FaceOwner sync request must use face location: " + own.field_name);
        }
    }
}

void Halo::upsert_halo_request_(const FieldHaloRequest &request)
{
    FieldHaloRequest normalized = request;
    if (normalized.sync_group.empty())
        normalized.sync_group = normalized.field_name;

    // 防止拼错名字导致 silent 插入
    if (!fld_->has_field(normalized.field_name))
        ERROR::Abort("[Halo] register_halo_field: field_name not registered in Field: " + normalized.field_name);

    // 升级策略：同名多次注册取更高等级，其他 metadata 必须一致。
    auto it = halo_registry_.find(normalized.field_name);
    if (it == halo_registry_.end())
    {
        halo_registry_.emplace(normalized.field_name, normalized);
        return;
    }

    FieldHaloRequest &old = it->second;
    if (!same_halo_request_metadata(old, normalized))
        ERROR::Abort("[Halo] register_halo_field: duplicate field has inconsistent halo request metadata: " + normalized.field_name);

    const int old_lv = static_cast<int>(old.level);
    const int new_lv = static_cast<int>(normalized.level);
    old.level = static_cast<HaloLevel>(std::max(old_lv, new_lv));
}

void Halo::register_halo_field(const FieldHaloRequest &request)
{
    upsert_halo_request_(request);
    rebuild_sync_registry_();
}

void Halo::register_halo_fields(const std::vector<FieldHaloRequest> &requests)
{
    for (const auto &request : requests)
        upsert_halo_request_(request);

    rebuild_sync_registry_();
}

void Halo::register_halo_field(const std::string &field_name, HaloLevel level)
{
    const FieldDescriptor &desc = fld_->descriptor(field_name);

    FieldHaloRequest request;
    request.field_name = field_name;
    request.sync_group = desc.sync.group.empty() ? field_name : desc.sync.group;
    request.location = desc.location;
    request.value_kind = desc.value_kind;
    request.ncomp = desc.ncomp;
    request.nghost = desc.nghost;
    request.level = level;
    request.owner_sync = desc.sync.owner_sync;
    request.orientation_aware = desc.sync.orientation_aware;

    register_halo_field(request);
}

void Halo::build_registered_patterns()
{
    rebuild_sync_registry_();

    // 清理旧 pattern（可重复 build）
    inner_patterns_.clear();
    parallel_patterns_.clear();
    inner_edge_patterns_.clear();
    parallel_edge_patterns_send.clear();
    parallel_edge_patterns_recv.clear();
    inner_vertex_patterns_.clear();
    parallel_vertex_patterns_send.clear();
    parallel_vertex_patterns_recv.clear();

    // coupling parallel corner patterns
    coupling_parallel_edge_patterns_send.clear();
    coupling_parallel_edge_patterns_recv.clear();
    coupling_parallel_vertex_patterns_send.clear();
    coupling_parallel_vertex_patterns_recv.clear();
    owner_sync_patterns_.clear();

    // 维度
    const int dim = fld_->grd->dimension;

    using PatternKey = std::pair<StaggerLocation, int>;
    std::set<PatternKey> face_keys, edge_keys, vertex_keys;

    // 1) 收集 keys
    for (const auto &kv : halo_registry_)
    {
        const FieldHaloRequest &req = kv.second;

        PatternKey key = {req.location, req.nghost};
        face_keys.insert(key);

        if (dim >= 2 && static_cast<int>(req.level) >= static_cast<int>(HaloLevel::Edge))
            edge_keys.insert(key);

        if (dim >= 3 && static_cast<int>(req.level) >= static_cast<int>(HaloLevel::Vertex))
            vertex_keys.insert(key);
    }

    // 2) Face patterns
    for (const auto &k : face_keys)
    {
        build_inner_1DCorner_pattern(k.first, k.second);
        build_parallel_1DCorner_pattern(k.first, k.second);
    }

    // 3) Edge patterns
    if (!edge_keys.empty())
    {
        for (const auto &k : edge_keys)
        {
            build_inner_2DCorner_pattern(k.first, k.second);
            build_parallel_2DCorner_pattern(k.first, k.second);
        }
    }

    // 4) Vertex patterns
    if (!vertex_keys.empty())
    {
        for (const auto &k : vertex_keys)
        {
            build_inner_3DCorner_pattern(k.first, k.second);
            build_parallel_3DCorner_pattern(k.first, k.second);
        }
    }

    // 5) Coupling parallel corner patterns (directed src -> dst)
    //    Build once here to avoid lazy rebuild during frequent coupling exchanges.
    const auto &cpairs = fld_->coupling_pairs();
    if (!cpairs.empty())
    {
        // small helpers: check whether any matching coupling patches exist on this rank
        auto has_parallel_coupling_edge = [&](const std::string &src, const std::string &dst) -> bool
        {
            const auto &parallel_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Parallel);
            for (const auto &ep : parallel_edges)
                if (ep.is_coupling && ep.nb_block_name == src && ep.this_block_name == dst)
                    return true;
            return false;
        };
        auto has_parallel_coupling_vertex = [&](const std::string &src, const std::string &dst) -> bool
        {
            const auto &parallel_vertices = TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Parallel);
            for (const auto &vp : parallel_vertices)
                if (vp.is_coupling && vp.nb_block_name == src && vp.this_block_name == dst)
                    return true;
            return false;
        };

        for (const auto &kv : cpairs)
        {
            const CouplingPairDesc &pd = kv.second;
            const std::string &src = pd.pair.src;
            const std::string &dst = pd.pair.dst;

            std::set<PatternKey> ckeys;
            for (const auto &ch : pd.channels)
            {
                ckeys.insert(PatternKey{ch.location, ch.nghost});
                if (ch.orientation_aware &&
                    (ch.value_kind == FieldValueKind::EdgeCovariant1Form ||
                     ch.value_kind == FieldValueKind::FaceContravariant2Form))
                {
                    const bool is_face = ch.value_kind == FieldValueKind::FaceContravariant2Form;
                    ckeys.insert(PatternKey{is_face ? StaggerLocation::FaceXi : StaggerLocation::EdgeXi, ch.nghost});
                    ckeys.insert(PatternKey{is_face ? StaggerLocation::FaceEt : StaggerLocation::EdgeEt, ch.nghost});
                    ckeys.insert(PatternKey{is_face ? StaggerLocation::FaceZe : StaggerLocation::EdgeZe, ch.nghost});
                }
            }

            if (dim >= 2 && has_parallel_coupling_edge(src, dst))
            {
                for (const auto &k : ckeys)
                    build_coupling_parallel_2DCorner_pattern(src, dst, k.first, k.second);
            }

            if (dim >= 3 && has_parallel_coupling_vertex(src, dst))
            {
                for (const auto &k : ckeys)
                    build_coupling_parallel_3DCorner_pattern(src, dst, k.first, k.second);
            }
        }
    }

    build_owner_sync_patterns_();
}
