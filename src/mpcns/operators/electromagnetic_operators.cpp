#include "mpcns/operators/electromagnetic_operators.hpp"

#include "mpcns/discrete/incidence_operator.hpp"

namespace mpcns::operators
{

void ElectromagneticOperators::compute_current(
    const state::PrimalFaceFormView &magnetic_flux,
    MagneticCurrentWorkspace workspace,
    state::PrimalEdgeFormView current) const
{
    hodge_.apply_star_2(magnetic_flux, workspace.star_2_flux);
    compute_current_from_star_2(workspace, current);
}

void ElectromagneticOperators::compute_current_from_star_2(
    MagneticCurrentWorkspace workspace,
    state::PrimalEdgeFormView current) const
{
    discrete::apply_d1_transpose(
        workspace.star_2_flux, workspace.curl_star_2_flux);
    hodge_.apply_inverse_star_1(workspace.curl_star_2_flux, current);
}

} // namespace mpcns::operators
