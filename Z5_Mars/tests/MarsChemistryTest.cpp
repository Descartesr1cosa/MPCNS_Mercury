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
    std::cout << "Mars chemistry tests passed\n";
}
