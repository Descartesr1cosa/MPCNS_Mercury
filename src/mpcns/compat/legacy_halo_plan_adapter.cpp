#include "mpcns/compat/legacy_halo_plan_adapter.hpp"

#include "3_field/FieldDescriptor.h"
#include "4_halo/Halo.h"

#include <stdexcept>
#include <utility>

namespace
{
using namespace mpcns;

parallel::FieldLocation convert_location(StaggerLocation location)
{
    switch (location)
    {
    case StaggerLocation::Cell: return parallel::FieldLocation::Cell;
    case StaggerLocation::Node: return parallel::FieldLocation::Node;
    case StaggerLocation::FaceXi: return parallel::FieldLocation::FaceXi;
    case StaggerLocation::FaceEt: return parallel::FieldLocation::FaceEta;
    case StaggerLocation::FaceZe: return parallel::FieldLocation::FaceZeta;
    case StaggerLocation::EdgeXi: return parallel::FieldLocation::EdgeXi;
    case StaggerLocation::EdgeEt: return parallel::FieldLocation::EdgeEta;
    case StaggerLocation::EdgeZe: return parallel::FieldLocation::EdgeZeta;
    }
    throw std::invalid_argument("unknown legacy field location");
}

core::LogicalBox convert_box(const Box3 &box)
{
    return {{box.lo.i, box.lo.j, box.lo.k},
            {box.hi.i, box.hi.j, box.hi.k}};
}

mesh::IndexTransform convert_transform(const TOPO::IndexTransform &transform)
{
    mesh::IndexTransform result;
    for (int axis = 0; axis < 3; ++axis)
    {
        result.permutation[axis] = transform.perm[axis];
        result.sign[axis] = transform.sign[axis];
    }
    result.offset = {transform.offset.i, transform.offset.j, transform.offset.k};
    return result;
}

parallel::HaloStage convert_stage(HaloLevel stage)
{
    if (stage == HaloLevel::Corner1D)
        return parallel::HaloStage::Face;
    if (stage == HaloLevel::Corner2D)
        return parallel::HaloStage::Edge;
    return parallel::HaloStage::Vertex;
}

parallel::ExchangeKind convert_kind(LegacyHaloPatternKind kind)
{
    switch (kind)
    {
    case LegacyHaloPatternKind::InnerCopy:
        return parallel::ExchangeKind::InnerCopy;
    case LegacyHaloPatternKind::ParallelFace:
        return parallel::ExchangeKind::ParallelFace;
    case LegacyHaloPatternKind::ParallelEdgeSend:
        return parallel::ExchangeKind::ParallelEdgeSend;
    case LegacyHaloPatternKind::ParallelEdgeRecv:
        return parallel::ExchangeKind::ParallelEdgeRecv;
    case LegacyHaloPatternKind::ParallelVertexSend:
        return parallel::ExchangeKind::ParallelVertexSend;
    case LegacyHaloPatternKind::ParallelVertexRecv:
        return parallel::ExchangeKind::ParallelVertexRecv;
    }
    throw std::invalid_argument("unknown legacy Halo pattern kind");
}

bool has_send_role(LegacyHaloPatternKind kind)
{
    return kind == LegacyHaloPatternKind::InnerCopy ||
           kind == LegacyHaloPatternKind::ParallelFace ||
           kind == LegacyHaloPatternKind::ParallelEdgeSend ||
           kind == LegacyHaloPatternKind::ParallelVertexSend;
}

bool has_recv_role(LegacyHaloPatternKind kind)
{
    return kind == LegacyHaloPatternKind::InnerCopy ||
           kind == LegacyHaloPatternKind::ParallelFace ||
           kind == LegacyHaloPatternKind::ParallelEdgeRecv ||
           kind == LegacyHaloPatternKind::ParallelVertexRecv;
}
} // namespace

namespace mpcns::compat
{
parallel::ExchangePlan LegacyHaloPlanAdapter::build(
    const mesh::MeshPartition &partition, const Halo &halo,
    const FieldHaloRequest &request)
{
    if (request.orientation_aware)
        throw std::invalid_argument("Unit 2A supports component-copy requests only");

    const auto patterns =
        halo.debug_component_copy_patterns(request.location, request.nghost);
    if (patterns.empty())
        throw std::invalid_argument("legacy Halo pattern key is not built");

    std::vector<parallel::ExchangeRegion> regions;
    for (const LegacyHaloPatternSnapshot &pattern : patterns)
    {
        for (std::size_t index = 0; index < pattern.regions.size(); ++index)
        {
            const HaloRegion &source = pattern.regions[index];
            if (source.this_rank != partition.rank() || source.this_block < 0 ||
                static_cast<std::size_t>(source.this_block) >=
                    partition.local_blocks().size())
                throw std::invalid_argument(
                    "legacy Halo region has invalid local block identity");
            if (source.neighbor_rank < 0 ||
                source.neighbor_rank >= partition.size() ||
                source.neighbor_block < 0)
                throw std::invalid_argument(
                    "legacy Halo region has invalid peer block identity");

            parallel::ExchangeRegion region;
            region.location = convert_location(request.location);
            region.ghost_layers = request.nghost;
            region.stage = convert_stage(pattern.stage);
            region.kind = convert_kind(pattern.kind);
            region.local_block = {source.this_rank, source.this_block};
            region.peer_block = {source.neighbor_rank, source.neighbor_block};
            if (has_send_role(pattern.kind))
                region.send_box = convert_box(source.send_box);
            if (has_recv_role(pattern.kind))
                region.recv_box = convert_box(source.recv_box);
            region.transform = convert_transform(source.trans);
            if (pattern.kind != LegacyHaloPatternKind::InnerCopy &&
                has_send_role(pattern.kind))
                region.send_tag = source.send_flag;
            if (pattern.kind != LegacyHaloPatternKind::InnerCopy &&
                has_recv_role(pattern.kind))
                region.recv_tag = source.recv_flag;
            region.pattern_order = index;
            regions.push_back(std::move(region));
        }
    }

    return parallel::ExchangePlan(convert_location(request.location),
                                  request.nghost, std::move(regions));
}
} // namespace mpcns::compat
