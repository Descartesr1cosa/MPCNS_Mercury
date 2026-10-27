#include "4_halo/Halo.h"

Halo::Halo(Field *field, TOPO::Topology *topo)
    : fld_(field), topo_(topo)
{
}

void Halo::set_topology_equiv(const TOPO::Topology *equiv)
{
    equiv_ = equiv;
}

const TOPO::Topology *Halo::topology_equiv() const
{
    return equiv_;
}

std::vector<HaloRegion> Halo::debug_halo_regions(StaggerLocation location,
                                                 int nghost,
                                                 HaloLevel stage) const
{
    const PatternKey key = {location, nghost};
    std::vector<HaloRegion> regions;

    auto append = [&regions](const std::map<PatternKey, HaloPattern> &patterns,
                             const PatternKey &pattern_key)
    {
        auto it = patterns.find(pattern_key);
        if (it == patterns.end())
            return;
        regions.insert(regions.end(), it->second.regions.begin(), it->second.regions.end());
    };

    if (stage == HaloLevel::Corner1D)
    {
        append(inner_patterns_, key);
        append(parallel_patterns_, key);
    }
    else if (stage == HaloLevel::Corner2D)
    {
        append(inner_edge_patterns_, key);
        append(parallel_edge_patterns_recv, key);
    }
    else if (stage == HaloLevel::Corner3D)
    {
        append(inner_vertex_patterns_, key);
        append(parallel_vertex_patterns_recv, key);
    }

    return regions;
}

std::vector<HaloRegion> Halo::debug_halo_send_regions(StaggerLocation location,
                                                       int nghost,
                                                       HaloLevel stage) const
{
    const PatternKey key = {location, nghost};
    std::vector<HaloRegion> regions;

    auto append = [&regions](const std::map<PatternKey, HaloPattern> &patterns,
                             const PatternKey &pattern_key)
    {
        auto it = patterns.find(pattern_key);
        if (it == patterns.end())
            return;
        regions.insert(regions.end(), it->second.regions.begin(), it->second.regions.end());
    };

    if (stage == HaloLevel::Corner1D)
    {
        append(parallel_patterns_, key);
    }
    else if (stage == HaloLevel::Corner2D)
    {
        append(parallel_edge_patterns_send, key);
    }
    else if (stage == HaloLevel::Corner3D)
    {
        append(parallel_vertex_patterns_send, key);
    }

    return regions;
}

std::vector<LegacyHaloPatternSnapshot>
Halo::debug_component_copy_patterns(StaggerLocation location, int nghost) const
{
    const PatternKey key{location, nghost};
    std::vector<LegacyHaloPatternSnapshot> result;
    const auto append = [&](const auto &patterns, LegacyHaloPatternKind kind, HaloLevel stage)
    {
        const auto it = patterns.find(key);
        if (it != patterns.end())
            result.push_back({kind, stage, it->second.location, it->second.nghost,
                              it->second.regions});
    };
    append(inner_patterns_, LegacyHaloPatternKind::InnerCopy, HaloLevel::Corner1D);
    append(parallel_patterns_, LegacyHaloPatternKind::ParallelFace, HaloLevel::Corner1D);
    append(inner_edge_patterns_, LegacyHaloPatternKind::InnerCopy, HaloLevel::Corner2D);
    append(parallel_edge_patterns_send, LegacyHaloPatternKind::ParallelEdgeSend, HaloLevel::Corner2D);
    append(parallel_edge_patterns_recv, LegacyHaloPatternKind::ParallelEdgeRecv, HaloLevel::Corner2D);
    append(inner_vertex_patterns_, LegacyHaloPatternKind::InnerCopy, HaloLevel::Corner3D);
    append(parallel_vertex_patterns_send, LegacyHaloPatternKind::ParallelVertexSend, HaloLevel::Corner3D);
    append(parallel_vertex_patterns_recv, LegacyHaloPatternKind::ParallelVertexRecv, HaloLevel::Corner3D);
    return result;
}

std::vector<LegacyOrientedTripletSnapshot>
Halo::debug_oriented_form_triplets() const
{
    std::vector<LegacyOrientedTripletSnapshot> result;
    const auto append = [&](const auto &groups, LegacyOrientedFormKind kind)
    {
        for (const auto &entry : groups)
        {
            const HaloTripletRequest &triplet = entry.second;
            LegacyOrientedTripletSnapshot snapshot;
            snapshot.group_name = triplet.group_name;
            snapshot.form_kind = kind;
            snapshot.field_names = {triplet.xi, triplet.eta, triplet.zeta};
            snapshot.level = triplet.level;
            snapshot.orientation_aware = triplet.orientation_aware;
            for (int axis = 0; axis < 3; ++axis)
            {
                const FieldDescriptor &descriptor =
                    fld_->descriptor(snapshot.field_names[axis]);
                snapshot.locations[axis] = descriptor.location;
                snapshot.value_kinds[axis] = descriptor.value_kind;
                snapshot.ncomp[axis] = descriptor.ncomp;
                snapshot.nghost[axis] = descriptor.nghost;
            }
            result.push_back(std::move(snapshot));
        }
    };

    append(edge_1form_triplets_, LegacyOrientedFormKind::Edge1Form);
    append(face_2form_triplets_, LegacyOrientedFormKind::Face2Form);
    return result;
}

void Halo::enable_edge_1form_corner_plan_capture(bool enabled)
{
    if (enabled)
        edge_1form_corner_plan_snapshots_.clear();
    capture_edge_1form_corner_plans_ = enabled;
}

std::optional<LegacyOrientedCornerPlanSnapshot>
Halo::debug_edge_1form_corner_plan(
    const std::string &group_name) const
{
    const auto it =
        edge_1form_corner_plan_snapshots_.find(group_name);
    if (it == edge_1form_corner_plan_snapshots_.end())
        return std::nullopt;
    return it->second;
}

void Halo::enable_vertex_edge_1form_plan_capture(bool enabled)
{
    if (enabled)
        vertex_edge_1form_plan_snapshots_.clear();
    capture_vertex_edge_1form_plans_ = enabled;
}

std::optional<LegacyOrientedVertexPlanSnapshot>
Halo::debug_vertex_edge_1form_plan(
    const std::string &group_name) const
{
    const auto it =
        vertex_edge_1form_plan_snapshots_.find(group_name);
    if (it == vertex_edge_1form_plan_snapshots_.end())
        return std::nullopt;
    return it->second;
}

void Halo::enable_edge_face_2form_plan_capture(bool enabled)
{
    if (enabled)
        edge_face_2form_plan_snapshots_.clear();
    capture_edge_face_2form_plans_ = enabled;
}

std::optional<LegacyFace2FormEdgePlanSnapshot>
Halo::debug_edge_face_2form_plan(
    const std::string &group_name) const
{
    const auto it =
        edge_face_2form_plan_snapshots_.find(group_name);
    if (it == edge_face_2form_plan_snapshots_.end())
        return std::nullopt;
    return it->second;
}

void Halo::enable_vertex_face_2form_plan_capture(bool enabled)
{
    if (enabled)
        vertex_face_2form_plan_snapshots_.clear();
    capture_vertex_face_2form_plans_ = enabled;
}

std::optional<LegacyFace2FormVertexPlanSnapshot>
Halo::debug_vertex_face_2form_plan(
    const std::string &group_name) const
{
    const auto it =
        vertex_face_2form_plan_snapshots_.find(group_name);
    if (it == vertex_face_2form_plan_snapshots_.end())
        return std::nullopt;
    return it->second;
}

void Halo::data_trans_1DCorner(std::string &field_name)
{
    exchange_inner(field_name);
    exchange_parallel(field_name);
}

void Halo::data_trans_2DCorner(std::string &field_name)
{
    exchange_inner_edge(field_name);
    exchange_parallel_edge(field_name);
}

void Halo::data_trans_3DCorner(std::string &field_name)
{
    exchange_inner_vertex(field_name);
    exchange_parallel_vertex(field_name);
}

void Halo::coupling_trans_1DCorner(std::string &src, std::string &dst)
{
    coupling_inner_face(src, dst);
    coupling_parallel_face(src, dst);
}

void Halo::coupling_trans_1DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_face(src, dst, field_cids);
    coupling_parallel_face(src, dst, field_cids);
}

void Halo::coupling_trans_2DCorner(std::string &src, std::string &dst)
{
    coupling_inner_edge(src, dst);
    coupling_parallel_edge(src, dst);
}

void Halo::coupling_trans_2DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_edge(src, dst, field_cids);
    coupling_parallel_edge(src, dst, field_cids);
}

void Halo::coupling_trans_3DCorner(std::string &src, std::string &dst)
{
    coupling_inner_vertex(src, dst);
    coupling_parallel_vertex(src, dst);
}

void Halo::coupling_trans_3DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_vertex(src, dst, field_cids);
    coupling_parallel_vertex(src, dst, field_cids);
}
