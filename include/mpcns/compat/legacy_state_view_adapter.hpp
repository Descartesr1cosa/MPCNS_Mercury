#pragma once

#include "mpcns/compat/legacy_state_layout.hpp"
#include "mpcns/state/state_view.hpp"

class Field;

namespace mpcns::compat
{

// Zero-copy adapter from the legacy FieldStorage backend. Views are invalid if
// the Field is destroyed, re-bound, or its storage is reallocated.
class LegacyStateViewAdapter
{
public:
    LegacyStateViewAdapter(Field &fields, const LegacyStateLayoutBinding &binding);

    state::StateView view(state::QuantityId quantity, int block) const;
    state::ConstStateView const_view(state::QuantityId quantity, int block) const;

private:
    Field *fields_ = nullptr;
    const LegacyStateLayoutBinding *binding_ = nullptr;
};

} // namespace mpcns::compat
