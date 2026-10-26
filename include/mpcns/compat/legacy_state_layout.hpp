#pragma once

#include "mpcns/state/state_layout.hpp"

#include <vector>

class Field;

namespace mpcns::compat
{

// Transitional binding only.  StateLayout is storage-independent; this maps
// its stable IDs back to the current FieldCatalog while legacy storage remains
// the production backend.
struct LegacyStateLayoutBinding
{
    state::StateLayout layout;
    std::vector<int> legacy_field_ids;
    std::vector<state::QuantityId> state_quantity_ids;

    int legacy_field_id(state::QuantityId id) const;
    state::QuantityId quantity_id_for_legacy_field_id(int legacy_field_id) const;
};

LegacyStateLayoutBinding borrow_legacy_state_layout(const Field &fields);

} // namespace mpcns::compat
