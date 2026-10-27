#include "mpcns/parallel/face_2form_edge_plan.hpp"

#include <cstdlib>
#include <map>
#include <stdexcept>
#include <utility>

namespace
{
using mpcns::core::LogicalBox;
using mpcns::mesh::IndexTransform;
using mpcns::parallel::FieldLocation;
using mpcns::parallel::OrientedCornerRole;
using mpcns::parallel::OrientedFormKind;
using mpcns::parallel::Face2FormEdgeTransfer;
using mpcns::parallel::inverse_permuted_axis;
using mpcns::parallel::oriented_form_sign;

bool valid_box(const LogicalBox &box)
{
    return box.hi.i >= box.lo.i &&
           box.hi.j >= box.lo.j &&
           box.hi.k >= box.lo.k;
}

long long box_volume(const LogicalBox &box)
{
    return static_cast<long long>(box.hi.i - box.lo.i) *
           static_cast<long long>(box.hi.j - box.lo.j) *
           static_cast<long long>(box.hi.k - box.lo.k);
}

bool valid_transform(const IndexTransform &transform)
{
    bool seen[3] = {false, false, false};
    for (int axis = 0; axis < 3; ++axis)
    {
        const int mapped = transform.permutation[axis];
        if (mapped < 0 || mapped > 2 || seen[mapped])
            return false;
        seen[mapped] = true;
        if (transform.sign[axis] != -1 &&
            transform.sign[axis] != 1)
            return false;
    }
    return true;
}

FieldLocation face_location(int axis)
{
    const FieldLocation locations[3] = {
        FieldLocation::FaceXi,
        FieldLocation::FaceEta,
        FieldLocation::FaceZeta};
    return locations[axis];
}

bool valid_direction(int direction)
{
    const int axis = std::abs(direction);
    return axis >= 1 && axis <= 3;
}

void validate_transfers(
    const std::vector<Face2FormEdgeTransfer> &transfers,
    OrientedCornerRole expected_role,
    int rank,
    int communicator_size,
    const std::array<mpcns::parallel::FormQuantityIdentity, 3> &quantities,
    int ncomp,
    int ghost_layers)
{
    std::map<int, std::size_t> next_peer_order;
    for (std::size_t order = 0; order < transfers.size(); ++order)
    {
        const Face2FormEdgeTransfer &transfer = transfers[order];
        if (transfer.stage != mpcns::parallel::HaloStage::Edge ||
            transfer.role != expected_role)
            throw std::invalid_argument(
                "Face2Form Edge transfer has wrong stage or role");
        if (!transfer.source_block.valid() ||
            !transfer.destination_block.valid())
            throw std::invalid_argument(
                "Face2Form Edge transfer has invalid block address");
        if (transfer.source_block.rank >= communicator_size ||
            transfer.destination_block.rank >= communicator_size ||
            transfer.source_block.rank ==
                transfer.destination_block.rank)
            throw std::invalid_argument(
                "Face2Form Edge transfer has invalid rank relation");

        const bool receive = expected_role == OrientedCornerRole::Receive;
        const int local_rank = receive
                                   ? transfer.destination_block.rank
                                   : transfer.source_block.rank;
        const int expected_peer = receive
                                      ? transfer.source_block.rank
                                      : transfer.destination_block.rank;
        if (local_rank != rank || transfer.peer_rank != expected_peer)
            throw std::invalid_argument(
                "Face2Form Edge transfer local/peer identity mismatch");

        if (transfer.send_axis < 0 || transfer.send_axis > 2 ||
            transfer.recv_axis < 0 || transfer.recv_axis > 2)
            throw std::invalid_argument(
                "Face2Form Edge transfer has invalid axis");
        if (!valid_transform(transfer.recv_to_send_transform) ||
            transfer.send_axis != inverse_permuted_axis(
                transfer.recv_to_send_transform,
                transfer.recv_axis))
            throw std::invalid_argument(
                "Face2Form Edge-stage send axis must use inverse permutation");
        if (transfer.orientation_sign != oriented_form_sign(
                OrientedFormKind::Face2Form,
                transfer.recv_to_send_transform,
                transfer.send_axis) ||
            (transfer.orientation_sign != -1 &&
             transfer.orientation_sign != 1))
            throw std::invalid_argument(
                "Face2Form Edge-stage orientation sign mismatch");

        if (transfer.source_quantity !=
                quantities[transfer.send_axis].name ||
            transfer.destination_quantity !=
                quantities[transfer.recv_axis].name)
            throw std::invalid_argument(
                "Face2Form Edge quantity/axis mismatch");
        if (transfer.ncomp != ncomp ||
            transfer.ghost_layers != ghost_layers)
            throw std::invalid_argument(
                "Face2Form Edge transfer shape mismatch");

        if (!valid_box(transfer.destination_recv_box))
            throw std::invalid_argument(
                "Face2Form Edge transfer has invalid receive box");
        if (!valid_box(transfer.source_node_box))
            throw std::invalid_argument(
                "Face2Form Edge transfer has invalid source node box");
        if (receive && transfer.source_send_box)
            throw std::invalid_argument(
                "receive Face2Form Edge transfer cannot expose a local send box");
        if (!receive &&
            (!transfer.source_send_box ||
             !valid_box(*transfer.source_send_box)))
            throw std::invalid_argument(
                "send Face2Form Edge transfer requires a valid source send box");

        if (!valid_direction(transfer.source_directions[0]) ||
            !valid_direction(transfer.source_directions[1]) ||
            std::abs(transfer.source_directions[0]) ==
                std::abs(transfer.source_directions[1]))
            throw std::invalid_argument(
                "Face2Form Edge transfer has invalid source directions");

        if (transfer.message_tag < 0 ||
            transfer.message_length < 0 ||
            transfer.message_length !=
                box_volume(transfer.destination_recv_box) * ncomp)
            throw std::invalid_argument(
                "Face2Form Edge transfer has invalid message contract");
        if (transfer.local_plan_order != order)
            throw std::invalid_argument(
                "Face2Form Edge local plan order is not contiguous");

        std::size_t &peer_order =
            next_peer_order[transfer.peer_rank];
        if (transfer.peer_message_order != peer_order ||
            transfer.message_tag !=
                static_cast<int>(peer_order))
            throw std::invalid_argument(
                "Face2Form Edge peer message order/tag mismatch");
        ++peer_order;
    }
}
} // namespace

namespace mpcns::parallel
{
Face2FormEdgeTransferPlan::Face2FormEdgeTransferPlan(
    int rank,
    int communicator_size,
    std::string sync_group,
    std::array<FormQuantityIdentity, 3> quantities,
    int ncomp,
    int ghost_layers,
    std::vector<Face2FormEdgeTransfer> receive_transfers,
    std::vector<Face2FormEdgeTransfer> send_transfers)
    : rank_(rank),
      communicator_size_(communicator_size),
      sync_group_(std::move(sync_group)),
      quantities_(std::move(quantities)),
      ncomp_(ncomp),
      ghost_layers_(ghost_layers),
      receive_transfers_(std::move(receive_transfers)),
      send_transfers_(std::move(send_transfers))
{
    if (communicator_size_ <= 0 ||
        rank_ < 0 || rank_ >= communicator_size_)
        throw std::invalid_argument(
            "Face2Form Edge plan has invalid communicator identity");
    if (sync_group_.empty())
        throw std::invalid_argument(
            "Face2Form Edge sync group is empty");
    if (ncomp_ <= 0 || ghost_layers_ < 0)
        throw std::invalid_argument(
            "Face2Form Edge plan has invalid shape");

    for (int axis = 0; axis < 3; ++axis)
    {
        const FormQuantityIdentity &quantity = quantities_[axis];
        if (quantity.name.empty() ||
            quantity.form_kind != OrientedFormKind::Face2Form ||
            quantity.location != face_location(axis) ||
            quantity.axis != axis)
            throw std::invalid_argument(
                "Face2Form Edge quantity identity is inconsistent");
    }

    validate_transfers(
        receive_transfers_, OrientedCornerRole::Receive,
        rank_, communicator_size_, quantities_, ncomp_, ghost_layers_);
    validate_transfers(
        send_transfers_, OrientedCornerRole::Send,
        rank_, communicator_size_, quantities_, ncomp_, ghost_layers_);
}
} // namespace mpcns::parallel
