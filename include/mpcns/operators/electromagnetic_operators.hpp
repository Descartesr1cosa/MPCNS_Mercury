#pragma once

#include "mpcns/discrete/lumped_hodge.hpp"
#include "mpcns/compat/legacy_cochain_view.hpp"

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

    void compute_current(const state::PrimalFaceFormView &magnetic_flux,
                         MagneticCurrentWorkspace workspace,
                         state::PrimalEdgeFormView current) const;

    // Continue from a caller-computed star_2 result. This keeps d1 transpose
    // and inverse star-1 in one legacy-compatible implementation while a
    // producer such as SnapshotLumpedHodge may write star_2 through StateView.
    void compute_current_from_star_2(MagneticCurrentWorkspace workspace,
                                     state::PrimalEdgeFormView current) const;

    // Fully StateView local continuation. The dual-face and primal-edge views
    // may intentionally alias the same borrowed Jedge allocation.
    void compute_current_from_star_2(
        const state::DualEdgeFormStateView &star_2_flux,
        state::DualFaceFormStateView curl_star_2_flux,
        state::PrimalEdgeFormStateView current) const;

private:
    void finish_current_from_curl_(state::DualFaceFormView curl_star_2_flux,
                                   state::PrimalEdgeFormView current) const;
    void finish_current_from_curl_(state::DualFaceFormStateView curl_star_2_flux,
                                   state::PrimalEdgeFormStateView current) const;

    const discrete::HodgeOperator &hodge_;
};

} // namespace mpcns::operators
