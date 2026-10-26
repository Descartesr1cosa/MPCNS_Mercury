#pragma once

#include "0_basic/StaggerLocation.h"

#include <cstddef>
#include <string>
#include <vector>

namespace mpcns::state
{

// Mathematical meaning is deliberately independent from storage and runtime
// exchange implementation.  Legacy FieldValueKind is mapped at compat/.
enum class QuantityKind
{
    Scalar,
    CartesianVector,
    ConservativeState,
    PrimalZeroForm,
    PrimalEdgeForm,
    PrimalFaceForm,
    GeometryMetric,
    Workspace
};

enum class OwnerPolicy
{
    None,
    Node,
    Edge,
    Face
};

struct SynchronizationContract
{
    std::string group;
    bool needs_coupling = false;
    bool needs_physical_boundary = false;
    bool needs_halo = false;
    HaloLevel halo_level = HaloLevel::Corner3D;
    OwnerPolicy owner_policy = OwnerPolicy::None;
    bool orientation_aware = false;
};

// This describes the stable mathematical/layout contract.  It owns neither
// data nor MPI plans; SimulationState and parallel adapters provide those.
struct QuantityDescriptor
{
    std::string canonical_name;
    std::string external_name;
    QuantityKind kind = QuantityKind::Scalar;
    StaggerLocation location = StaggerLocation::Cell;
    int component_count = 1;
    int ghost_layers = 0;
    std::string allocation_domain;
    SynchronizationContract synchronization;
};

class QuantityId
{
public:
    constexpr QuantityId() = default;
    explicit constexpr QuantityId(int value) : value_(value) {}

    constexpr int value() const { return value_; }
    constexpr bool is_valid() const { return value_ >= 0; }

private:
    int value_ = -1;
};

// Immutable description of a state backend.  It has no dependency on Field,
// Grid, MPI, or an array implementation, and can therefore outlive the
// legacy storage migration.
class StateLayout
{
public:
    StateLayout() = default;
    explicit StateLayout(std::vector<QuantityDescriptor> quantities);

    std::size_t size() const { return quantities_.size(); }
    const QuantityDescriptor &descriptor(QuantityId id) const;
    QuantityId find(const std::string &canonical_name) const;
    const std::vector<QuantityDescriptor> &quantities() const { return quantities_; }

private:
    std::vector<QuantityDescriptor> quantities_;
};

} // namespace mpcns::state
