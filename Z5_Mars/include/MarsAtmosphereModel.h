#pragma once

#include "MarsChemistry.h"
#include "MarsSpecies.h"

#include <array>

namespace MARS
{
struct AtmosphereConfig
{
    double mars_radius_m = 3.396e6;
    double lower_boundary_altitude_km = 100.0;
    double heliocentric_distance_au = 1.524;
    std::array<double, 3> sun_direction = {1.0, 0.0, 0.0};
    double ion_temperature_K = 1500.0;
    double electron_temperature_K = 1500.0;
};

struct AtmosphereState
{
    std::array<double, NumNeutral> neutral_cm3{};
    std::array<double, NumPhoto> photo_cm3_s{};
};

struct PceSolution
{
    std::array<double, NumSpecies> ion_cm3{};
    double relative_residual = 0.0;
    int iterations = 0;
    bool converged = false;
};

// Grid-independent neutral atmosphere, grey-EUV attenuation and lower-boundary
// photochemical-equilibrium model. Positions are Mars-centred Cartesian metres.
class AtmosphereModel
{
public:
    explicit AtmosphereModel(AtmosphereConfig config = {});

    const AtmosphereConfig &config() const { return config_; }
    double NeutralTemperature(double altitude_km) const;
    std::array<double, NumNeutral> NeutralDensity(double altitude_km) const;
    AtmosphereState Evaluate(const std::array<double, 3> &position_m) const;
    PceSolution WallIonDensity(const std::array<double, 3> &position_m) const;

    bool IsInOpticalShadow(const std::array<double, 3> &position_m) const;
    double OpticalDepth(const std::array<double, 3> &position_m) const;

private:
    AtmosphereConfig config_;
};

PceSolution SolvePhotochemicalEquilibrium(const AtmosphereState &state,
                                           double ion_temperature_K,
                                           double electron_temperature_K);
}
