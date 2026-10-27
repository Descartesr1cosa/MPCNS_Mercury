#include "mpcns/compat/legacy_oriented_form_plan_adapter.hpp"

#include "4_halo/Halo.h"

#include <stdexcept>
#include <utility>

namespace
{
using namespace mpcns;

parallel::OrientedFormKind convert_kind(LegacyOrientedFormKind kind)
{
    return kind == LegacyOrientedFormKind::Edge1Form
               ? parallel::OrientedFormKind::Edge1Form
               : parallel::OrientedFormKind::Face2Form;
}

parallel::FieldLocation convert_location(StaggerLocation location)
{
    switch (location)
    {
    case StaggerLocation::EdgeXi: return parallel::FieldLocation::EdgeXi;
    case StaggerLocation::EdgeEt: return parallel::FieldLocation::EdgeEta;
    case StaggerLocation::EdgeZe: return parallel::FieldLocation::EdgeZeta;
    case StaggerLocation::FaceXi: return parallel::FieldLocation::FaceXi;
    case StaggerLocation::FaceEt: return parallel::FieldLocation::FaceEta;
    case StaggerLocation::FaceZe: return parallel::FieldLocation::FaceZeta;
    default:
        throw std::invalid_argument(
            "oriented-form triplet has a non-form location");
    }
}

mesh::IndexTransform convert_transform(const TOPO::IndexTransform &transform)
{
    mesh::IndexTransform result;
    for (int axis = 0; axis < 3; ++axis)
    {
        result.permutation[axis] = transform.perm[axis];
        result.sign[axis] = transform.sign[axis];
    }
    result.offset = {transform.offset.i, transform.offset.j,
                     transform.offset.k};
    return result;
}

core::LogicalBox convert_box(const Box3 &box)
{
    return {{box.lo.i, box.lo.j, box.lo.k},
            {box.hi.i, box.hi.j, box.hi.k}};
}

LegacyHaloPatternSnapshot face_pattern(
    const Halo &halo,
    StaggerLocation location,
    int ghost_layers,
    LegacyHaloPatternKind wanted)
{
    const auto patterns =
        halo.debug_component_copy_patterns(location, ghost_layers);
    for (const LegacyHaloPatternSnapshot &pattern : patterns)
        if (pattern.stage == HaloLevel::Corner1D &&
            pattern.kind == wanted)
            return pattern;

    throw std::invalid_argument(
        "required legacy Face-level Halo pattern is missing");
}

FieldValueKind expected_value_kind(LegacyOrientedFormKind kind)
{
    return kind == LegacyOrientedFormKind::Edge1Form
               ? FieldValueKind::EdgeCovariant1Form
               : FieldValueKind::FaceContravariant2Form;
}
} // namespace

namespace mpcns::compat
{
parallel::OrientedFormTransferPlan LegacyOrientedFormPlanAdapter::build(
    const mesh::MeshPartition &partition,
    const Halo &halo,
    const FieldHaloRequest &request)
{
    if (!request.do_halo || !request.orientation_aware)
        throw std::invalid_argument(
            "request is not an orientation-aware Halo transfer");
    if (request.value_kind != FieldValueKind::EdgeCovariant1Form &&
        request.value_kind != FieldValueKind::FaceContravariant2Form)
        throw std::invalid_argument(
            "orientation-aware request has unsupported value kind");

    const std::string group =
        request.sync_group.empty() ? request.field_name
                                   : request.sync_group;
    return build_group(partition, halo, group);
}

parallel::OrientedFormTransferPlan LegacyOrientedFormPlanAdapter::build_group(
    const mesh::MeshPartition &partition,
    const Halo &halo,
    const std::string &sync_group)
{
    const auto groups = halo.debug_oriented_form_triplets();
    const LegacyOrientedTripletSnapshot *legacy = nullptr;
    for (const LegacyOrientedTripletSnapshot &candidate : groups)
        if (candidate.group_name == sync_group)
        {
            legacy = &candidate;
            break;
        }

    if (legacy == nullptr)
        throw std::invalid_argument(
            "oriented-form triplet group is not registered");
    if (!legacy->orientation_aware)
        throw std::invalid_argument(
            "triplet group is not orientation-aware");

    const FieldValueKind value_kind =
        expected_value_kind(legacy->form_kind);
    const int ncomp = legacy->ncomp[0];
    const int ghost_layers = legacy->nghost[0];
    std::array<parallel::FormQuantityIdentity, 3> quantities;
    for (int axis = 0; axis < 3; ++axis)
    {
        if (legacy->field_names[axis].empty())
            throw std::invalid_argument("triplet member is missing");
        if (legacy->value_kinds[axis] != value_kind)
            throw std::invalid_argument("triplet mixes field value kinds");
        if (legacy->ncomp[axis] != ncomp ||
            legacy->nghost[axis] != ghost_layers)
            throw std::invalid_argument(
                "triplet members have inconsistent shape metadata");

        quantities[axis] = {
            legacy->field_names[axis],
            convert_kind(legacy->form_kind),
            convert_location(legacy->locations[axis]),
            axis};
    }

    std::vector<parallel::OrientedFormTransfer> transfers;
    const parallel::OrientedFormKind form_kind =
        convert_kind(legacy->form_kind);
    std::size_t transfer_order = 0;

    for (LegacyHaloPatternKind legacy_kind :
         {LegacyHaloPatternKind::InnerCopy,
          LegacyHaloPatternKind::ParallelFace})
    {
        for (int destination_axis = 0; destination_axis < 3;
             ++destination_axis)
        {
            const LegacyHaloPatternSnapshot destination_pattern =
                face_pattern(halo, legacy->locations[destination_axis],
                             ghost_layers, legacy_kind);

            for (std::size_t region_index = 0;
                 region_index < destination_pattern.regions.size();
                 ++region_index)
            {
                const HaloRegion &destination_region =
                    destination_pattern.regions[region_index];
                const mesh::IndexTransform transform =
                    convert_transform(destination_region.trans);
                const int source_axis =
                    parallel::inverse_permuted_axis(
                        transform, destination_axis);
                const LegacyHaloPatternSnapshot source_pattern =
                    face_pattern(halo, legacy->locations[source_axis],
                                 ghost_layers, legacy_kind);
                if (region_index >= source_pattern.regions.size())
                    throw std::invalid_argument(
                        "source and destination patterns have different order/size");

                const HaloRegion &source_region =
                    source_pattern.regions[region_index];
                if (source_region.this_rank !=
                        destination_region.this_rank ||
                    source_region.this_block !=
                        destination_region.this_block ||
                    source_region.neighbor_rank !=
                        destination_region.neighbor_rank ||
                    source_region.neighbor_block !=
                        destination_region.neighbor_block)
                    throw std::invalid_argument(
                        "source/destination pattern identities do not match");

                if (source_region.this_rank != partition.rank() ||
                    source_region.this_block < 0 ||
                    static_cast<std::size_t>(source_region.this_block) >=
                        partition.local_blocks().size())
                    throw std::invalid_argument(
                        "oriented transfer has invalid local source block");
                if (destination_region.neighbor_rank < 0 ||
                    destination_region.neighbor_rank >= partition.size() ||
                    destination_region.neighbor_block < 0)
                    throw std::invalid_argument(
                        "oriented transfer has invalid destination block");

                parallel::OrientedFormTransfer transfer;
                transfer.kind =
                    legacy_kind == LegacyHaloPatternKind::InnerCopy
                        ? parallel::OrientedTransferKind::Inner
                        : parallel::OrientedTransferKind::Parallel;
                transfer.source_block = {
                    source_region.this_rank, source_region.this_block};
                transfer.destination_block = {
                    destination_region.neighbor_rank,
                    destination_region.neighbor_block};
                transfer.source_axis = source_axis;
                transfer.destination_axis = destination_axis;
                transfer.source_quantity =
                    legacy->field_names[source_axis];
                transfer.destination_quantity =
                    legacy->field_names[destination_axis];
                transfer.ncomp = ncomp;
                transfer.ghost_layers = ghost_layers;
                transfer.source_region_order = region_index;
                transfer.destination_region_order = region_index;
                transfer.transfer_order = transfer_order++;
                transfer.send_box =
                    convert_box(source_region.send_box);
                transfer.recv_box =
                    convert_box(destination_region.recv_box);
                transfer.transform = transform;
                transfer.orientation_sign =
                    parallel::oriented_form_sign(
                        form_kind, transform, source_axis);
                if (transfer.kind ==
                    parallel::OrientedTransferKind::Parallel)
                {
                    transfer.send_tag =
                        destination_region.send_flag;
                    transfer.recv_tag =
                        destination_region.recv_flag;
                }
                transfers.push_back(std::move(transfer));
            }
        }
    }

    return parallel::OrientedFormTransferPlan(
        legacy->group_name, form_kind, std::move(quantities),
        ncomp, ghost_layers, std::move(transfers));
}
} // namespace mpcns::compat
