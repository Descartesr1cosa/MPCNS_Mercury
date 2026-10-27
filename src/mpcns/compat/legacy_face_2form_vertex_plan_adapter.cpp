#include "mpcns/compat/legacy_face_2form_vertex_plan_adapter.hpp"

#include "4_halo/Halo.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace
{
using namespace mpcns;

core::LogicalBox convert_box(const Box3 &box)
{
    return {{box.lo.i, box.lo.j, box.lo.k},
            {box.hi.i, box.hi.j, box.hi.k}};
}

mesh::IndexTransform convert_transform(
    const TOPO::IndexTransform &transform)
{
    mesh::IndexTransform result;
    for (int axis = 0; axis < 3; ++axis)
    {
        result.permutation[axis] = transform.perm[axis];
        result.sign[axis] = transform.sign[axis];
    }
    result.offset = {
        transform.offset.i,
        transform.offset.j,
        transform.offset.k};
    return result;
}

parallel::Face2FormVertexTransfer convert_transfer(
    const LegacyFace2FormVertexTransferSnapshot &legacy,
    const std::array<std::string, 3> &field_names)
{
    parallel::Face2FormVertexTransfer transfer;
    transfer.role =
        legacy.role == LegacyOrientedCornerRole::Receive
            ? parallel::OrientedCornerRole::Receive
            : parallel::OrientedCornerRole::Send;
    transfer.source_block = {
        legacy.send_rank, legacy.send_block};
    transfer.destination_block = {
        legacy.recv_rank, legacy.recv_block};
    transfer.send_axis = legacy.send_axis;
    transfer.recv_axis = legacy.recv_axis;
    if (legacy.send_axis >= 0 && legacy.send_axis < 3)
        transfer.source_quantity = field_names[legacy.send_axis];
    if (legacy.recv_axis >= 0 && legacy.recv_axis < 3)
        transfer.destination_quantity =
            field_names[legacy.recv_axis];
    transfer.orientation_sign = legacy.sign;
    transfer.destination_recv_box =
        convert_box(legacy.recv_box);
    if (legacy.send_box)
        transfer.source_send_box =
            convert_box(*legacy.send_box);
    transfer.source_node_box =
        convert_box(legacy.node_box_on_send);
    transfer.source_directions = {
        legacy.dir1_send, legacy.dir2_send,
        legacy.dir3_send};
    transfer.recv_to_send_transform =
        convert_transform(legacy.trans_recv_to_send);
    transfer.peer_rank = legacy.peer_rank;
    transfer.message_tag = legacy.tag;
    transfer.ncomp = legacy.ncomp;
    transfer.ghost_layers = legacy.nghost;
    transfer.message_length = legacy.message_length;
    transfer.local_plan_order = legacy.local_plan_order;
    transfer.peer_message_order =
        legacy.peer_message_order;
    return transfer;
}

void validate_local_identity(
    const mesh::MeshPartition &partition,
    const LegacyFace2FormVertexTransferSnapshot &legacy)
{
    const bool receive =
        legacy.role == LegacyOrientedCornerRole::Receive;
    const int local_rank =
        receive ? legacy.recv_rank : legacy.send_rank;
    const int local_block =
        receive ? legacy.recv_block : legacy.send_block;
    const int remote_rank =
        receive ? legacy.send_rank : legacy.recv_rank;
    const int remote_block =
        receive ? legacy.send_block : legacy.recv_block;

    if (local_rank != partition.rank() ||
        local_block < 0 ||
        static_cast<std::size_t>(local_block) >=
            partition.local_blocks().size())
        throw std::invalid_argument(
            "Face2Form Vertex snapshot has invalid local block identity");
    if (remote_rank < 0 ||
        remote_rank >= partition.size() ||
        remote_rank == partition.rank() ||
        remote_block < 0)
        throw std::invalid_argument(
            "Face2Form Vertex snapshot has invalid remote block identity");
}
} // namespace

namespace mpcns::compat
{
parallel::Face2FormVertexTransferPlan
LegacyFace2FormVertexPlanAdapter::build(
    const mesh::MeshPartition &partition,
    const Halo &halo,
    const FieldHaloRequest &request)
{
    if (!request.do_halo ||
        !request.orientation_aware ||
        request.value_kind !=
            FieldValueKind::FaceContravariant2Form)
        throw std::invalid_argument(
            "request is not an orientation-aware Face2Form transfer");
    if (static_cast<int>(request.level) <
        static_cast<int>(HaloLevel::Corner3D))
        throw std::invalid_argument(
            "request does not include the Vertex stage");

    const std::string group =
        request.sync_group.empty()
            ? request.field_name
            : request.sync_group;
    const auto captured =
        halo.debug_vertex_face_2form_plan(group);
    if (!captured)
        throw std::invalid_argument(
            "Vertex-stage Face2Form metadata has not been captured");
    if (captured->group_name != group ||
        captured->ncomp != request.ncomp ||
        captured->nghost != request.nghost)
        throw std::invalid_argument(
            "captured Face2Form Vertex metadata does not match the request");

    if (std::find(
            captured->field_names.begin(),
            captured->field_names.end(),
            request.field_name) ==
        captured->field_names.end())
        throw std::invalid_argument(
            "request field is not a captured triplet member");

    const parallel::FieldLocation locations[3] = {
        parallel::FieldLocation::FaceXi,
        parallel::FieldLocation::FaceEta,
        parallel::FieldLocation::FaceZeta};
    std::array<parallel::FormQuantityIdentity, 3>
        quantities;
    for (int axis = 0; axis < 3; ++axis)
    {
        if (captured->field_names[axis].empty())
            throw std::invalid_argument(
                "captured Face2Form triplet member is missing");
        quantities[axis] = {
            captured->field_names[axis],
            parallel::OrientedFormKind::Face2Form,
            locations[axis],
            axis};
    }

    std::vector<parallel::Face2FormVertexTransfer>
        receive_transfers;
    receive_transfers.reserve(
        captured->recv_plan.size());
    for (const auto &legacy : captured->recv_plan)
    {
        if (legacy.role !=
            LegacyOrientedCornerRole::Receive)
            throw std::invalid_argument(
                "receive snapshot has wrong role");
        validate_local_identity(partition, legacy);
        receive_transfers.push_back(
            convert_transfer(
                legacy, captured->field_names));
    }

    std::vector<parallel::Face2FormVertexTransfer>
        send_transfers;
    send_transfers.reserve(
        captured->send_plan.size());
    for (const auto &legacy : captured->send_plan)
    {
        if (legacy.role !=
            LegacyOrientedCornerRole::Send)
            throw std::invalid_argument(
                "send snapshot has wrong role");
        validate_local_identity(partition, legacy);
        send_transfers.push_back(
            convert_transfer(
                legacy, captured->field_names));
    }

    return parallel::Face2FormVertexTransferPlan(
        partition.rank(), partition.size(), group,
        std::move(quantities), captured->ncomp,
        captured->nghost,
        std::move(receive_transfers),
        std::move(send_transfers));
}
} // namespace mpcns::compat
