#include "MarsChemistry.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
    MARS::ChemistryInput in;
    in.neutral_cm3.fill(1.0e6);
    in.ion_cm3.fill(1.0e2);
    in.photo_cm3_s.fill(1.0);
    in.ion_temperature_K.fill(1000.0);
    in.electron_temperature_K = 1500.0;
    const auto rates = MARS::EvaluateChemistry(in);
    for (std::size_t s = 0; s < MARS::NumSpecies; ++s)
    {
        assert(std::isfinite(rates.production_cm3_s[s]));
        assert(std::isfinite(rates.loss_cm3_s[s]));
        assert(rates.production_cm3_s[s] >= 0.0);
        assert(rates.loss_cm3_s[s] >= 0.0);
    }
    assert(rates.electron_production_cm3_s == 14.0);

    MARS::ChemistryInput empty;
    empty.ion_temperature_K.fill(50.0);
    const auto zero = MARS::EvaluateChemistry(empty);
    for (std::size_t s = 0; s < MARS::NumSpecies; ++s)
    {
        assert(zero.production_cm3_s[s] == 0.0);
        assert(zero.loss_cm3_s[s] == 0.0);
    }

    MARS::ChemistryInput stiff=in;
    stiff.neutral_cm3.fill(1.0e12);
    stiff.photo_cm3_s.fill(1.0e3);
    const auto initial=stiff.ion_cm3;
    std::array<double,MARS::NumSpecies> floors;
    floors.fill(1.0e-20);
    constexpr double stiff_dt=1.0e-2;
    const auto implicit=MARS::SolvePointImplicitChemistry(stiff,stiff_dt,floors);
    assert(implicit.converged);
    for(std::size_t s=0;s<MARS::NumSpecies;++s)
    {
        assert(std::isfinite(implicit.ion_cm3[s]));
        assert(implicit.ion_cm3[s]>=floors[s]);
        const double residual=implicit.ion_cm3[s]-initial[s]-stiff_dt*(
            implicit.rates.production_cm3_s[s]-implicit.rates.loss_cm3_s[s]);
        const double scale=1.0+std::abs(initial[s])+std::abs(implicit.ion_cm3[s]);
        assert(std::abs(residual)/scale<1.0e-8);
    }

    // A species already at the numerical floor may have a net destruction
    // rate.  The constrained backward-Euler solution must keep that species
    // on the active lower bound instead of reporting a false Newton failure.
    MARS::ChemistryInput bounded{};
    bounded.neutral_cm3.fill(1.0e12);
    bounded.ion_cm3.fill(1.0e-11);
    bounded.ion_temperature_K.fill(1500.0);
    bounded.electron_temperature_K=1500.0;
    std::array<double,MARS::NumSpecies> active_floors;
    active_floors.fill(1.0e-11);
    const auto active=MARS::SolvePointImplicitChemistry(bounded,1.0e-3,active_floors);
    assert(active.converged);
    bool active_bound_seen=false;
    for(std::size_t s=0;s<MARS::NumSpecies;++s)
    {
        assert(active.ion_cm3[s]>=active_floors[s]);
        if(active.ion_cm3[s]<=active_floors[s]*(1.0+1.0e-8))
            active_bound_seen=true;
    }
    assert(active_bound_seen);
    std::cout << "Mars chemistry tests passed\n";
}
