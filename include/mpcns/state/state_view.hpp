#pragma once

#include "mpcns/core/logical_index.hpp"

#include <cstddef>
#include <stdexcept>
#include <type_traits>

namespace mpcns::state
{

// A non-owning strided view over one structured quantity block.  It never
// allocates or copies. The caller keeps the backing storage alive and stable
// for the entire lifetime of this view.
template <class Value>
class BasicStateView
{
public:
    static_assert(std::is_same<Value, double>::value ||
                      std::is_same<Value, const double>::value,
                  "StateView currently exposes double-valued production storage");

    using value_type = Value;

    constexpr BasicStateView() = default;

    BasicStateView(Value *logical_origin, core::LogicalBox extent,
                   core::LogicalBox inner_extent, int component_count,
                   std::ptrdiff_t stride_i, std::ptrdiff_t stride_j,
                   std::ptrdiff_t stride_k, std::ptrdiff_t stride_component)
        : logical_origin_(logical_origin), extent_(extent), inner_extent_(inner_extent),
          component_count_(component_count), stride_i_(stride_i), stride_j_(stride_j),
          stride_k_(stride_k), stride_component_(stride_component)
    {
        if (!logical_origin_ || component_count_ <= 0 || stride_component_ <= 0)
            throw std::invalid_argument("StateView requires valid contiguous component storage");
    }

    bool is_bound() const { return logical_origin_ != nullptr; }
    const core::LogicalBox &extent() const { return extent_; }
    const core::LogicalBox &inner_extent() const { return inner_extent_; }
    int component_count() const { return component_count_; }

    // Unchecked, inline hot-path access. Callers must use indices inside
    // extent() and a component in [0, component_count()).
    Value &operator()(int i, int j, int k, int component = 0) const
    {
        return logical_origin_[i * stride_i_ + j * stride_j_ +
                               k * stride_k_ + component * stride_component_];
    }

    // Checked access for setup, adapters, and tests.
    Value &at(int i, int j, int k, int component = 0) const
    {
        const core::LogicalIndex index{i, j, k};
        if (!is_bound() || !extent_.contains(index) ||
            component < 0 || component >= component_count_)
            throw std::out_of_range("StateView index lies outside its logical extent");
        return (*this)(i, j, k, component);
    }

private:
    Value *logical_origin_ = nullptr;
    core::LogicalBox extent_{};
    core::LogicalBox inner_extent_{};
    int component_count_ = 0;
    std::ptrdiff_t stride_i_ = 0;
    std::ptrdiff_t stride_j_ = 0;
    std::ptrdiff_t stride_k_ = 0;
    std::ptrdiff_t stride_component_ = 1;
};

using StateView = BasicStateView<double>;
using ConstStateView = BasicStateView<const double>;

} // namespace mpcns::state
