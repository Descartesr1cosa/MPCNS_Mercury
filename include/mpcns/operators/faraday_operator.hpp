#pragma once
#include "mpcns/state/cochain_view.hpp"
namespace mpcns::operators
{
class FaradayOperator
{
public:
    static void compute_rhs(const state::PrimalEdgeFormView &emf, state::PrimalFaceFormView rhs);
    static void compute_rhs(const state::PrimalEdgeFormStateView &emf,
                            state::PrimalFaceFormStateView rhs);
    static void advance(const state::PrimalEdgeFormView &emf, state::PrimalFaceFormView flux, double dt);
};
}
