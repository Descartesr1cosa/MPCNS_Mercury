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
    if (state_->layout().descriptor(handle.quantity()).storage_role != StorageRole::Workspace)
        throw std::invalid_argument("Workspace accepts only Workspace-role quantities");
}

void Workspace::validate_view(StateHandle handle, int local_block_id,
                              const StateView &view) const
{
    if (local_block_id < 0 || !view.is_bound())
        throw std::invalid_argument("Workspace requires a bound StateView");
    const QuantityDescriptor &quantity = state_->layout().descriptor(handle.quantity());
    if (view.quantity().value() != handle.quantity().value() ||
        view.location() != quantity.location ||
        view.component_count() != quantity.component_count)
        throw std::invalid_argument("Workspace View does not match its StateHandle layout");
}

void Workspace::bind(StateHandle handle, int local_block_id, StateView view)
{
    validate_workspace_handle(handle);
    validate_view(handle, local_block_id, view);

    const auto existing = std::find_if(bindings_.begin(), bindings_.end(),
        [handle, local_block_id](const Binding &binding)
        { return binding.quantity.value() == handle.quantity().value() &&
                 binding.local_block_id == local_block_id; });
    if (existing != bindings_.end())
        existing->view = view;
    else
        bindings_.push_back({handle.quantity(), local_block_id, view});
}

bool Workspace::is_bound(StateHandle handle, int local_block_id) const
{
    validate_workspace_handle(handle);
    if (local_block_id < 0)
        throw std::out_of_range("invalid Workspace local block id");
    return std::any_of(bindings_.begin(), bindings_.end(),
        [handle, local_block_id](const Binding &binding)
        { return binding.quantity.value() == handle.quantity().value() &&
                 binding.local_block_id == local_block_id; });
}

StateView Workspace::view(StateHandle handle, int local_block_id) const
{
    validate_workspace_handle(handle);
    if (local_block_id < 0)
        throw std::out_of_range("invalid Workspace local block id");
    const auto existing = std::find_if(bindings_.begin(), bindings_.end(),
        [handle, local_block_id](const Binding &binding)
        { return binding.quantity.value() == handle.quantity().value() &&
                 binding.local_block_id == local_block_id; });
    if (existing == bindings_.end())
        throw std::logic_error("Workspace quantity is not bound");
    return existing->view;
}

void Workspace::release(StateHandle handle, int local_block_id)
{
    validate_workspace_handle(handle);
    if (local_block_id < 0)
        throw std::out_of_range("invalid Workspace local block id");
    bindings_.erase(std::remove_if(bindings_.begin(), bindings_.end(),
        [handle, local_block_id](const Binding &binding)
        { return binding.quantity.value() == handle.quantity().value() &&
                 binding.local_block_id == local_block_id; }),
        bindings_.end());
}

void Workspace::clear()
{
    bindings_.clear();
}

} // namespace mpcns::state
