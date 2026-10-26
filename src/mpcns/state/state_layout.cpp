#include "mpcns/state/state_layout.hpp"

#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace mpcns::state
{
namespace
{

bool is_edge_location(StaggerLocation location)
{
    return location == StaggerLocation::EdgeXi ||
           location == StaggerLocation::EdgeEt ||
           location == StaggerLocation::EdgeZe;
}

bool is_face_location(StaggerLocation location)
{
    return location == StaggerLocation::FaceXi ||
           location == StaggerLocation::FaceEt ||
           location == StaggerLocation::FaceZe;
}

void validate_quantity(const QuantityDescriptor &quantity)
{
    if (quantity.canonical_name.empty() || quantity.external_name.empty())
        throw std::invalid_argument("state quantity names must not be empty");
    if (quantity.component_count <= 0 || quantity.ghost_layers < 0)
        throw std::invalid_argument("state quantity has invalid component or ghost count: " +
                                    quantity.canonical_name);
    if (quantity.kind == QuantityKind::PrimalZeroForm &&
        quantity.location != StaggerLocation::Node)
        throw std::invalid_argument("primal zero-form must be node-centred: " +
                                    quantity.canonical_name);
    if (quantity.kind == QuantityKind::PrimalEdgeForm && !is_edge_location(quantity.location))
        throw std::invalid_argument("primal edge-form must be edge-centred: " +
                                    quantity.canonical_name);
    if (quantity.kind == QuantityKind::PrimalFaceForm && !is_face_location(quantity.location))
        throw std::invalid_argument("primal face-form must be face-centred: " +
                                    quantity.canonical_name);
}

}

StateLayout::StateLayout(std::vector<QuantityDescriptor> quantities)
    : quantities_(std::move(quantities))
{
    std::unordered_set<std::string> canonical_names;
    std::unordered_set<std::string> external_names;
    for (const QuantityDescriptor &quantity : quantities_)
    {
        validate_quantity(quantity);
        if (!canonical_names.insert(quantity.canonical_name).second)
            throw std::invalid_argument("duplicate state canonical name: " +
                                        quantity.canonical_name);
        if (!external_names.insert(quantity.external_name).second)
            throw std::invalid_argument("duplicate state external name: " +
                                        quantity.external_name);
    }
}

const QuantityDescriptor &StateLayout::descriptor(QuantityId id) const
{
    if (!id.is_valid() || static_cast<std::size_t>(id.value()) >= quantities_.size())
        throw std::out_of_range("state quantity id lies outside the layout");
    return quantities_[static_cast<std::size_t>(id.value())];
}

QuantityId StateLayout::find(const std::string &canonical_name) const
{
    for (std::size_t index = 0; index < quantities_.size(); ++index)
        if (quantities_[index].canonical_name == canonical_name)
            return QuantityId{static_cast<int>(index)};
    return QuantityId{};
}

} // namespace mpcns::state
