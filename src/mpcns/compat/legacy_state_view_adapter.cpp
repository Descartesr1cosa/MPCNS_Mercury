#include "mpcns/compat/legacy_state_view_adapter.hpp"

#include "3_field/Field.h"

#include <stdexcept>

namespace mpcns::compat
{
namespace
{

core::LogicalBox logical_box(const Int3 &lo, const Int3 &hi)
{
    return core::LogicalBox{{lo.i, lo.j, lo.k}, {hi.i, hi.j, hi.k}};
}

void validate_binding(const FieldBlock &field, const state::QuantityDescriptor &quantity)
{
    if (!field.is_allocated())
        throw std::logic_error("cannot borrow StateView from inactive field: " +
                               quantity.canonical_name);
    const FieldDescriptor &legacy = field.descriptor();
    if (legacy.location != quantity.location || legacy.ncomp != quantity.component_count ||
        legacy.nghost != quantity.ghost_layers)
        throw std::logic_error("legacy FieldBlock no longer matches StateLayout: " +
                               quantity.canonical_name);
}

template <class Value>
state::BasicStateView<Value> make_view(const FieldBlock &field,
                                       const state::QuantityDescriptor &quantity,
                                       state::QuantityId quantity_id)
{
    validate_binding(field, quantity);
    const Vector &array = field.raw_data();
    const Int3 lo = field.get_lo();
    const Int3 hi = field.get_hi();
    const Int3 inner_lo = field.inner_lo();
    const Int3 inner_hi = field.inner_hi();
    const Value *origin = array.LogicalOriginData();
    return state::BasicStateView<Value>{origin, logical_box(lo, hi),
                                        logical_box(inner_lo, inner_hi),
                                        quantity.component_count,
                                        array.StrideI(), array.StrideJ(),
                                        array.StrideK(), array.StrideComponent(),
                                        quantity.location, quantity_id};
}

template <>
state::StateView make_view<double>(const FieldBlock &field,
                                   const state::QuantityDescriptor &quantity,
                                   state::QuantityId quantity_id)
{
    validate_binding(field, quantity);
    FieldBlock &mutable_field = const_cast<FieldBlock &>(field);
    Vector &array = mutable_field.raw_data();
    const Int3 lo = field.get_lo();
    const Int3 hi = field.get_hi();
    const Int3 inner_lo = field.inner_lo();
    const Int3 inner_hi = field.inner_hi();
    return state::StateView{array.LogicalOriginData(), logical_box(lo, hi),
                            logical_box(inner_lo, inner_hi), quantity.component_count,
                            array.StrideI(), array.StrideJ(), array.StrideK(),
                            array.StrideComponent(), quantity.location, quantity_id};
}

}

LegacyStateViewAdapter::LegacyStateViewAdapter(
    Field &fields, const LegacyStateLayoutBinding &binding)
    : fields_(&fields), binding_(&binding)
{
}

state::StateView LegacyStateViewAdapter::view(state::QuantityId quantity, int block) const
{
    if (!fields_ || !binding_ || block < 0 || block >= fields_->num_blocks())
        throw std::out_of_range("invalid StateView quantity or block request");
    return make_view<double>(fields_->field(binding_->legacy_field_id(quantity), block),
                             binding_->layout.descriptor(quantity), quantity);
}

state::ConstStateView LegacyStateViewAdapter::const_view(state::QuantityId quantity, int block) const
{
    if (!fields_ || !binding_ || block < 0 || block >= fields_->num_blocks())
        throw std::out_of_range("invalid const StateView quantity or block request");
    return make_view<const double>(fields_->field(binding_->legacy_field_id(quantity), block),
                                   binding_->layout.descriptor(quantity), quantity);
}

} // namespace mpcns::compat
