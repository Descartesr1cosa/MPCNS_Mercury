#include "mpcns/operators/faraday_operator.hpp"
#include "mpcns/discrete/incidence_operator.hpp"
namespace mpcns::operators
{
void FaradayOperator::compute_rhs(const state::PrimalEdgeFormView &emf,
                                  state::PrimalFaceFormView rhs)
{
    discrete::apply_d1(emf, rhs, -1.0);
}

void FaradayOperator::compute_rhs(const state::PrimalEdgeFormStateView &emf,
                                  state::PrimalFaceFormStateView rhs)
{
    discrete::apply_d1(emf, rhs, -1.0);
}
void FaradayOperator::advance(const state::PrimalEdgeFormView &emf,
                              state::PrimalFaceFormView flux, double dt)
{
    discrete::subtract_d1(emf, flux, dt);
}
}
