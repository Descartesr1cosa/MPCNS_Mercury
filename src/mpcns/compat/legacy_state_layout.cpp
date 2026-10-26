#include "mpcns/compat/legacy_state_layout.hpp"

#include "3_field/Field.h"

#include <stdexcept>
#include <utility>

namespace mpcns::compat
{
namespace
{

state::QuantityKind quantity_kind(FieldValueKind kind)
{
    switch (kind)
    {
    case FieldValueKind::Scalar: return state::QuantityKind::Scalar;
    case FieldValueKind::CartesianVector: return state::QuantityKind::CartesianVector;
    case FieldValueKind::ConservedVector: return state::QuantityKind::ConservativeState;
    case FieldValueKind::EdgeCovariant1Form: return state::QuantityKind::PrimalEdgeForm;
    case FieldValueKind::FaceContravariant2Form: return state::QuantityKind::PrimalFaceForm;
    case FieldValueKind::GeometryMetric: return state::QuantityKind::GeometryMetric;
    // These legacy tags did not encode a mathematical representation.  Keep
    // their historical scalar fallback while their storage purpose is mapped
    // independently below.
    case FieldValueKind::CouplingBufferOnly:
    case FieldValueKind::Temporary: return state::QuantityKind::Scalar;
    }
    throw std::logic_error("unrecognised legacy field value kind");
}

state::StorageRole storage_role(const FieldDescriptor &legacy)
{
    switch (legacy.storage_role)
    {
    case FieldStorageRole::State: break;
    case FieldStorageRole::Derived: return state::StorageRole::Derived;
    case FieldStorageRole::Workspace: return state::StorageRole::Workspace;
    case FieldStorageRole::CouplingBuffer: return state::StorageRole::CouplingBuffer;
    }
    // Preserve meaningful roles for unannotated historical descriptors.
    if (legacy.value_kind == FieldValueKind::Temporary)
        return state::StorageRole::Workspace;
    if (legacy.value_kind == FieldValueKind::CouplingBufferOnly)
        return state::StorageRole::CouplingBuffer;
    if (legacy.value_kind == FieldValueKind::GeometryMetric)
        return state::StorageRole::Derived;
    return state::StorageRole::State;
}

state::OwnerPolicy owner_policy(OwnerSyncPolicy policy)
{
    switch (policy)
    {
    case OwnerSyncPolicy::None: return state::OwnerPolicy::None;
    case OwnerSyncPolicy::NodeOwner: return state::OwnerPolicy::Node;
    case OwnerSyncPolicy::EdgeOwner: return state::OwnerPolicy::Edge;
    case OwnerSyncPolicy::FaceOwner: return state::OwnerPolicy::Face;
    }
    throw std::logic_error("unrecognised legacy owner-sync policy");
}

state::QuantityDescriptor quantity_descriptor(const FieldDescriptor &legacy)
{
    return state::QuantityDescriptor{
        legacy.name, legacy.name, quantity_kind(legacy.value_kind), storage_role(legacy), legacy.location,
        legacy.ncomp, legacy.nghost, legacy.physics,
        state::SynchronizationContract{legacy.sync.group, legacy.sync.do_coupling,
                                       legacy.sync.do_physical, legacy.sync.do_halo,
                                       legacy.sync.halo_level,
                                       owner_policy(legacy.sync.owner_sync),
                                       legacy.sync.orientation_aware}};
}

}

int LegacyStateLayoutBinding::legacy_field_id(state::QuantityId id) const
{
    if (!id.is_valid() || static_cast<std::size_t>(id.value()) >= legacy_field_ids.size())
        throw std::out_of_range("state quantity id has no legacy Field binding");
    return legacy_field_ids[static_cast<std::size_t>(id.value())];
}

LegacyStateLayoutBinding borrow_legacy_state_layout(const Field &fields)
{
    std::vector<state::QuantityDescriptor> quantities;
    std::vector<int> ids;
    quantities.reserve(static_cast<std::size_t>(fields.num_fields()));
    ids.reserve(static_cast<std::size_t>(fields.num_fields()));
    for (int field_id = 0; field_id < fields.num_fields(); ++field_id)
    {
        quantities.push_back(quantity_descriptor(fields.descriptor(field_id)));
        ids.push_back(field_id);
    }
    return LegacyStateLayoutBinding{state::StateLayout{std::move(quantities)},
                                    std::move(ids)};
}

} // namespace mpcns::compat
