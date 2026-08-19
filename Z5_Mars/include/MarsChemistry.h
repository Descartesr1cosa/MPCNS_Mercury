#pragma once

#include "MarsSpecies.h"
#include <array>

namespace MARS
{
struct ChemistryInput
{
    std::array<double, NumNeutral> neutral_cm3{};
    std::array<double, NumSpecies> ion_cm3{};
    std::array<double, NumPhoto> photo_cm3_s{};
    std::array<double, NumSpecies> ion_temperature_K{};
    double electron_temperature_K = 50.0;
};

struct ChemistryRates
{
    std::array<double, NumSpecies> production_cm3_s{};
    std::array<double, NumSpecies> loss_cm3_s{};
    double electron_production_cm3_s = 0.0;
    double electron_loss_cm3_s = 0.0;
};

struct PointImplicitChemistryResult
{
    std::array<double, NumSpecies> ion_cm3{};
    ChemistryRates rates;
    int iterations = 0;
    bool converged = false;
};

ChemistryRates EvaluateChemistry(const ChemistryInput &in);
PointImplicitChemistryResult SolvePointImplicitChemistry(
    const ChemistryInput &input, double dt_seconds,
    const std::array<double, NumSpecies> &density_floor_cm3,
    int max_iterations = 10, double relative_tolerance = 1.0e-10);
}
