#pragma once

#include "mpcns/discrete/lumped_hodge.hpp"

namespace mpcns::operators
{

// Caller-owned scratch keeps this operator allocation-free and makes the
// synchronization boundary explicit. Input flux/EMF must have valid halos.
struct MagneticCurrentWorkspace
{
    state::DualEdgeFormView star_2_flux;
    state::DualFaceFormView curl_star_2_flux;
};

// Reusable local DEC composition for induction and current. It contains no
// Field catalog, Grid, topology construction, MPI, or solver-time decision.
class ElectromagneticOperators
{
public:
    explicit ElectromagneticOperators(const discrete::HodgeOperator &hodge)
        : hodge_(hodge)
    {
    }

    void advance_magnetic_flux(const state::PrimalEdgeFormView &emf,
                               state::PrimalFaceFormView magnetic_flux,
                               double dt) const;

    void compute_current(const state::PrimalFaceFormView &magnetic_flux,
                         MagneticCurrentWorkspace workspace,
                         state::PrimalEdgeFormView current) const;

private:
    const discrete::HodgeOperator &hodge_;
};

} // namespace mpcns::operators
