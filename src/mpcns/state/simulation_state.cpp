#include "mpcns/state/simulation_state.hpp"

#include <algorithm>
#include <stdexcept>

namespace mpcns::state
{

SimulationState::SimulationState(const StateLayout &layout)
    : layout_(&layout)
{
}

const StateLayout &SimulationState::layout() const
{
    if (!layout_)
        throw std::logic_error("SimulationState has no StateLayout");
    return *layout_;
}

StateHandle SimulationState::handle(QuantityId quantity) const
{
    (void)layout().descriptor(quantity);
    return StateHandle{this, quantity};
}

Workspace::Workspace(const SimulationState &state)
    : state_(&state)
{
}

void Workspace::validate_workspace_handle(StateHandle handle) const
{
    if (!state_ || !handle.is_valid() || handle.owner_ != state_)
        throw std::out_of_range("invalid Workspace StateHandle");
    if (state_->layout().descriptor(handle.quantity()).kind != QuantityKind::Workspace)
        throw std::invalid_argument("Workspace accepts only Workspace quantities");
}

void Workspace::bind(StateHandle handle, StateView view)
{
    validate_workspace_handle(handle);
    if (!view.is_bound())
        throw std::invalid_argument("Workspace requires a bound StateView");

    const auto existing = std::find_if(bindings_.begin(), bindings_.end(),
        [handle](const Binding &binding)
        { return binding.quantity.value() == handle.quantity().value(); });
    if (existing != bindings_.end())
        existing->view = view;
    else
        bindings_.push_back({handle.quantity(), view});
}

bool Workspace::is_bound(StateHandle handle) const
{
    validate_workspace_handle(handle);
    return std::any_of(bindings_.begin(), bindings_.end(),
        [handle](const Binding &binding)
        { return binding.quantity.value() == handle.quantity().value(); });
}

StateView Workspace::view(StateHandle handle) const
{
    validate_workspace_handle(handle);
    const auto existing = std::find_if(bindings_.begin(), bindings_.end(),
        [handle](const Binding &binding)
        { return binding.quantity.value() == handle.quantity().value(); });
    if (existing == bindings_.end())
        throw std::logic_error("Workspace quantity is not bound");
    return existing->view;
}

void Workspace::release(StateHandle handle)
{
    validate_workspace_handle(handle);
    bindings_.erase(std::remove_if(bindings_.begin(), bindings_.end(),
        [handle](const Binding &binding)
        { return binding.quantity.value() == handle.quantity().value(); }),
        bindings_.end());
}

void Workspace::clear()
{
    bindings_.clear();
}

} // namespace mpcns::state
