#pragma once

#include "mpcns/state/state_layout.hpp"
#include "mpcns/state/state_view.hpp"

#include <vector>

namespace mpcns::state
{

class SimulationState;
class Workspace;

// A checked, typed reference to a layout quantity.  It owns no storage and
// remains valid only while the StateLayout used to create it remains alive.
class StateHandle
{
public:
    constexpr StateHandle() = default;

    constexpr QuantityId quantity() const { return quantity_; }
    constexpr bool is_valid() const { return owner_ != nullptr && quantity_.is_valid(); }

private:
    friend class SimulationState;
    friend class Workspace;
    constexpr StateHandle(const SimulationState *owner, QuantityId quantity)
        : owner_(owner), quantity_(quantity) {}

    const SimulationState *owner_ = nullptr;
    QuantityId quantity_{};
};

// Lightweight state identity/lifetime anchor.  Production arrays remain in
// the selected backend (currently legacy FieldStorage); this object merely
// provides checked handles into its immutable StateLayout.
class SimulationState
{
public:
    explicit SimulationState(const StateLayout &layout);

    const StateLayout &layout() const;
    StateHandle handle(QuantityId quantity) const;

private:
    const StateLayout *layout_ = nullptr;
};

// Explicitly managed, non-owning scratch bindings.  Workspace never allocates
// or copies an array: callers bind a StateView borrowed from the backend, then
// release it before the backend can be reconfigured or destroyed.
class Workspace
{
public:
    explicit Workspace(const SimulationState &state);

    void bind(StateHandle handle, StateView view);
    bool is_bound(StateHandle handle) const;
    StateView view(StateHandle handle) const;
    void release(StateHandle handle);
    void clear();

private:
    struct Binding
    {
        QuantityId quantity;
        StateView view;
    };

    const SimulationState *state_ = nullptr;
    std::vector<Binding> bindings_;

    void validate_workspace_handle(StateHandle handle) const;
};

} // namespace mpcns::state
