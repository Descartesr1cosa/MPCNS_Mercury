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
    SimulationState(const SimulationState &) = delete;
    SimulationState &operator=(const SimulationState &) = delete;
    SimulationState(SimulationState &&) = delete;
    SimulationState &operator=(SimulationState &&) = delete;

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
    Workspace(const Workspace &) = delete;
    Workspace &operator=(const Workspace &) = delete;
    Workspace(Workspace &&) = delete;
    Workspace &operator=(Workspace &&) = delete;

    // Binding the same (handle, block) replaces only that borrowed View after
    // validation. These setup-time operations may grow bindings_; numerical
    // kernels only retrieve an already-bound view.
    void bind(StateHandle handle, int local_block_id, StateView view);
    bool is_bound(StateHandle handle, int local_block_id) const;
    StateView view(StateHandle handle, int local_block_id) const;
    void release(StateHandle handle, int local_block_id);
    void clear();

private:
    struct Binding
    {
        QuantityId quantity;
        int local_block_id = -1;
        StateView view;
    };

    const SimulationState *state_ = nullptr;
    std::vector<Binding> bindings_;

    void validate_workspace_handle(StateHandle handle) const;
    void validate_view(StateHandle handle, int local_block_id, const StateView &view) const;
};

} // namespace mpcns::state
