#pragma once

#include "mpcns/parallel/oriented_corner_plan.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace mpcns::parallel
{
struct Face2FormEdgeTransfer
{
    HaloStage stage = HaloStage::Edge;
    OrientedCornerRole role = OrientedCornerRole::Receive;

    mesh::BlockAddress source_block;
    mesh::BlockAddress destination_block;
    int send_axis = -1;
    int recv_axis = -1;
    std::string source_quantity;
    std::string destination_quantity;
    int orientation_sign = 0;

    // The destination receive box controls pack and unpack traversal. Legacy
    // metadata retains valid empty half-open boxes for zero-length messages.
    core::LogicalBox destination_recv_box;
    // Valid only for Send records. It records the legacy box computed on the
    // source rank, although packing maps destination_recv_box coordinates.
    std::optional<core::LogicalBox> source_send_box;
    core::LogicalBox source_node_box;
    std::array<int, 2> source_directions{};
    mesh::IndexTransform recv_to_send_transform;

    int peer_rank = -1;
    int message_tag = -1;
    int ncomp = 0;
    int ghost_layers = 0;
    int message_length = 0;

    std::size_t local_plan_order = 0;
    std::size_t peer_message_order = 0;
    OrientedCornerPackingOrder packing_order =
        OrientedCornerPackingOrder::DestinationRecvIJKComponent;
};

class Face2FormEdgeTransferPlan
{
public:
    Face2FormEdgeTransferPlan(
        int rank,
        int communicator_size,
        std::string sync_group,
        std::array<FormQuantityIdentity, 3> quantities,
        int ncomp,
        int ghost_layers,
        std::vector<Face2FormEdgeTransfer> receive_transfers,
        std::vector<Face2FormEdgeTransfer> send_transfers);

    int rank() const noexcept { return rank_; }
    int communicator_size() const noexcept { return communicator_size_; }
    const std::string &sync_group() const noexcept { return sync_group_; }
    OrientedFormKind form_kind() const noexcept
    {
        return OrientedFormKind::Face2Form;
    }
    HaloStage stage() const noexcept { return HaloStage::Edge; }
    const std::array<FormQuantityIdentity, 3> &quantities() const noexcept
    {
        return quantities_;
    }
    int ncomp() const noexcept { return ncomp_; }
    int ghost_layers() const noexcept { return ghost_layers_; }
    const std::vector<Face2FormEdgeTransfer> &
    receive_transfers() const noexcept
    {
        return receive_transfers_;
    }
    const std::vector<Face2FormEdgeTransfer> &
    send_transfers() const noexcept
    {
        return send_transfers_;
    }

private:
    int rank_;
    int communicator_size_;
    std::string sync_group_;
    std::array<FormQuantityIdentity, 3> quantities_;
    int ncomp_;
    int ghost_layers_;
    std::vector<Face2FormEdgeTransfer> receive_transfers_;
    std::vector<Face2FormEdgeTransfer> send_transfers_;
};
} // namespace mpcns::parallel
