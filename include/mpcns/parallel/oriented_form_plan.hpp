#pragma once

#include "mpcns/parallel/exchange_plan.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace mpcns::parallel
{
enum class OrientedFormKind { Edge1Form, Face2Form };
enum class OrientedTransferKind { Inner, Parallel };
enum class PackingOrder { TransformedDestinationIJKComponent };

struct FormQuantityIdentity
{
    std::string name;
    OrientedFormKind form_kind = OrientedFormKind::Edge1Form;
    FieldLocation location = FieldLocation::Cell;
    int axis = -1;
};

struct OrientedFormTransfer
{
    HaloStage stage = HaloStage::Face;
    OrientedTransferKind kind = OrientedTransferKind::Inner;

    mesh::BlockAddress source_block;
    mesh::BlockAddress destination_block;

    int source_axis = -1;
    int destination_axis = -1;
    std::string source_quantity;
    std::string destination_quantity;

    int ncomp = 0;
    int ghost_layers = 0;

    std::size_t source_region_order = 0;
    std::size_t destination_region_order = 0;
    std::size_t transfer_order = 0;

    // For Inner, source/destination are the directed local copy endpoints.
    // For Parallel, they are the outbound local source and peer destination.
    // send_box is in source_block coordinates; recv_box is the local inbound
    // box from the same bidirectional legacy Face pattern.
    core::LogicalBox send_box;
    core::LogicalBox recv_box;
    mesh::IndexTransform transform;

    int orientation_sign = 0;
    std::optional<int> send_tag;
    std::optional<int> recv_tag;

    PackingOrder packing_order =
        PackingOrder::TransformedDestinationIJKComponent;
};

class OrientedFormTransferPlan
{
public:
    OrientedFormTransferPlan(
        std::string sync_group,
        OrientedFormKind form_kind,
        std::array<FormQuantityIdentity, 3> quantities,
        int ncomp,
        int ghost_layers,
        std::vector<OrientedFormTransfer> transfers);

    const std::string &sync_group() const noexcept { return sync_group_; }
    OrientedFormKind form_kind() const noexcept { return form_kind_; }
    const std::array<FormQuantityIdentity, 3> &quantities() const noexcept
    {
        return quantities_;
    }
    int ncomp() const noexcept { return ncomp_; }
    int ghost_layers() const noexcept { return ghost_layers_; }
    const std::vector<OrientedFormTransfer> &transfers() const noexcept
    {
        return transfers_;
    }

private:
    std::string sync_group_;
    OrientedFormKind form_kind_;
    std::array<FormQuantityIdentity, 3> quantities_;
    int ncomp_;
    int ghost_layers_;
    std::vector<OrientedFormTransfer> transfers_;
};

int inverse_permuted_axis(const mesh::IndexTransform &transform,
                          int destination_axis);
int oriented_form_sign(OrientedFormKind kind,
                       const mesh::IndexTransform &transform,
                       int source_axis);
} // namespace mpcns::parallel
