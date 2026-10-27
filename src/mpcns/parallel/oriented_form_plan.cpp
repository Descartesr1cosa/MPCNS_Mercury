#include "mpcns/parallel/oriented_form_plan.hpp"

#include <stdexcept>
#include <utility>

namespace
{
using mpcns::parallel::FieldLocation;
using mpcns::parallel::OrientedFormKind;

FieldLocation expected_location(OrientedFormKind kind, int axis)
{
    if (kind == OrientedFormKind::Edge1Form)
    {
        const FieldLocation locations[3] = {
            FieldLocation::EdgeXi, FieldLocation::EdgeEta,
            FieldLocation::EdgeZeta};
        return locations[axis];
    }

    const FieldLocation locations[3] = {
        FieldLocation::FaceXi, FieldLocation::FaceEta,
        FieldLocation::FaceZeta};
    return locations[axis];
}

bool valid_transform(const mpcns::mesh::IndexTransform &transform)
{
    bool seen[3] = {false, false, false};
    for (int axis = 0; axis < 3; ++axis)
    {
        const int mapped = transform.permutation[axis];
        if (mapped < 0 || mapped > 2 || seen[mapped])
            return false;
        seen[mapped] = true;
        if (transform.sign[axis] != -1 && transform.sign[axis] != 1)
            return false;
    }
    return true;
}
} // namespace

namespace mpcns::parallel
{
int inverse_permuted_axis(const mesh::IndexTransform &transform,
                          int destination_axis)
{
    if (destination_axis < 0 || destination_axis > 2 ||
        !valid_transform(transform))
        throw std::invalid_argument("invalid oriented-form axis transform");

    for (int source_axis = 0; source_axis < 3; ++source_axis)
        if (transform.permutation[source_axis] == destination_axis)
            return source_axis;

    throw std::invalid_argument("destination axis is absent from permutation");
}

int oriented_form_sign(OrientedFormKind kind,
                       const mesh::IndexTransform &transform,
                       int source_axis)
{
    if (source_axis < 0 || source_axis > 2 || !valid_transform(transform))
        throw std::invalid_argument("invalid oriented-form sign transform");

    if (kind == OrientedFormKind::Edge1Form)
        return transform.sign[source_axis];

    const int determinant_sign =
        transform.sign[0] * transform.sign[1] * transform.sign[2];
    return determinant_sign * transform.sign[source_axis] >= 0 ? 1 : -1;
}

OrientedFormTransferPlan::OrientedFormTransferPlan(
    std::string sync_group,
    OrientedFormKind form_kind,
    std::array<FormQuantityIdentity, 3> quantities,
    int ncomp,
    int ghost_layers,
    std::vector<OrientedFormTransfer> transfers)
    : sync_group_(std::move(sync_group)),
      form_kind_(form_kind),
      quantities_(std::move(quantities)),
      ncomp_(ncomp),
      ghost_layers_(ghost_layers),
      transfers_(std::move(transfers))
{
    if (sync_group_.empty())
        throw std::invalid_argument("oriented-form sync group is empty");
    if (ncomp_ <= 0)
        throw std::invalid_argument("oriented-form ncomp must be positive");
    if (ghost_layers_ < 0)
        throw std::invalid_argument("oriented-form ghost depth is negative");

    for (int axis = 0; axis < 3; ++axis)
    {
        const FormQuantityIdentity &quantity = quantities_[axis];
        if (quantity.name.empty() || quantity.form_kind != form_kind_ ||
            quantity.axis != axis ||
            quantity.location != expected_location(form_kind_, axis))
            throw std::invalid_argument(
                "oriented-form quantity identity is incomplete or inconsistent");
    }

    for (std::size_t order = 0; order < transfers_.size(); ++order)
    {
        const OrientedFormTransfer &transfer = transfers_[order];
        if (transfer.stage != HaloStage::Face)
            throw std::invalid_argument("Unit 2B-1 supports Face stage only");
        if (!transfer.source_block.valid() ||
            !transfer.destination_block.valid())
            throw std::invalid_argument(
                "oriented-form transfer requires valid block addresses");
        if (transfer.source_axis < 0 || transfer.source_axis > 2 ||
            transfer.destination_axis < 0 || transfer.destination_axis > 2)
            throw std::invalid_argument("oriented-form transfer has invalid axis");
        if (transfer.source_quantity !=
                quantities_[transfer.source_axis].name ||
            transfer.destination_quantity !=
                quantities_[transfer.destination_axis].name)
            throw std::invalid_argument(
                "oriented-form transfer quantity/axis mismatch");
        if (transfer.ncomp != ncomp_ ||
            transfer.ghost_layers != ghost_layers_)
            throw std::invalid_argument(
                "oriented-form transfer shape contract mismatch");
        if (!valid_transform(transfer.transform) ||
            inverse_permuted_axis(transfer.transform,
                                  transfer.destination_axis) !=
                transfer.source_axis)
            throw std::invalid_argument(
                "oriented-form transfer axis permutation mismatch");
        if (transfer.orientation_sign !=
            oriented_form_sign(form_kind_, transfer.transform,
                               transfer.source_axis))
            throw std::invalid_argument(
                "oriented-form transfer orientation sign mismatch");
        if (transfer.transfer_order != order)
            throw std::invalid_argument(
                "oriented-form transfer order is not contiguous");

        const bool parallel =
            transfer.kind == OrientedTransferKind::Parallel;
        if (parallel && (!transfer.send_tag || !transfer.recv_tag))
            throw std::invalid_argument(
                "parallel oriented-form transfer requires both tags");
        if (!parallel && (transfer.send_tag || transfer.recv_tag))
            throw std::invalid_argument(
                "inner oriented-form transfer cannot carry MPI tags");
        if (parallel &&
            transfer.source_block.rank == transfer.destination_block.rank)
            throw std::invalid_argument(
                "parallel oriented-form transfer must cross ranks");
        if (!parallel &&
            transfer.source_block.rank != transfer.destination_block.rank)
            throw std::invalid_argument(
                "inner oriented-form transfer must stay on one rank");
    }
}
} // namespace mpcns::parallel
