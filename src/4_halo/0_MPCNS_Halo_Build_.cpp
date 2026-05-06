#include "4_halo/1_MPCNS_Halo.h"
#include "0_basic/Error.h"

#include <algorithm>
#include <set>

namespace
{
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
}

void Halo::register_halo_field(const FieldHaloRequest &request)
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
    }
    else
    {
        FieldHaloRequest &old = it->second;
        if (!same_halo_request_metadata(old, normalized))
            ERROR::Abort("[Halo] register_halo_field: duplicate field has inconsistent halo request metadata: " + normalized.field_name);

        const int old_lv = static_cast<int>(old.level);
        const int new_lv = static_cast<int>(normalized.level);
        old.level = static_cast<HaloLevel>(std::max(old_lv, new_lv));
    }
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
            for (const auto &ep : topo_->parallel_edge_patches)
                if (ep.is_coupling && ep.nb_block_name == src && ep.this_block_name == dst)
                    return true;
            return false;
        };
        auto has_parallel_coupling_vertex = [&](const std::string &src, const std::string &dst) -> bool
        {
            for (const auto &vp : topo_->parallel_vertex_patches)
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
                ckeys.insert(PatternKey{ch.location, ch.nghost});

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
}
