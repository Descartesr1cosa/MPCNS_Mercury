#include "MarsAtmosphereModel.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace MARS
{
namespace
{
constexpr double Boltzmann = 1.380649e-23;
constexpr double MarsGM = 4.282837e13;

// Reference densities at 100 km [cm^-3]. The molecular minor species retain
// the original model's lower-atmosphere composition while their vertical
// continuation is now computed by diffusive equilibrium.
constexpr std::array<double, NumSpecies> BaseDensity = {
    1.80e10, 4.49e11, 1.93e10, 3.49e2, 4.67e5,
    2.32e10, 4.50e6, 1.04e6, 3.98e9, 1.03e9};

// Optically-thin effective EUVAC channel frequencies [s^-1] at the reference
// Mars orbit. These preserve the fourteen-channel legacy EUVAC preprocessing.
constexpr std::array<double, NumPhoto> ChannelFrequency = {
    2.5826e-7, 4.5064e-7, 2.9182e-7, 4.0303e-7, 5.1862e-8,
    2.8442e-7, 1.8133e-7, 4.9796e-7, 2.1321e-7, 4.2292e-7,
    2.4412e-8, 4.1207e-8, 3.2316e-8, 5.1480e-8};
constexpr std::array<std::size_t, NumPhoto> PhotoParent = {
    Ar, CO2, CO, C, H, N2, N, NO, O, O2, CO2, N2, CO2, CO2};

// Band-integrated photoabsorption cross sections [cm^2]. Only abundant
// absorbers are retained in the grey reduction of the 37-bin EUVAC model.
constexpr std::array<double, NumSpecies> AbsorptionCrossSection = {
    1.75e-18, 3.50e-18, 2.625e-18, 0.0, 0.0,
    2.625e-18, 0.0, 0.0, 1.75e-18, 3.50e-18};

double clamp_exp(double x)
{
    return std::exp(std::clamp(x, -700.0, 700.0));
}

double norm3(const std::array<double, 3> &x)
{
    return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

double dot3(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
}

AtmosphereModel::AtmosphereModel(AtmosphereConfig config) : config_(config)
{
    const double sn = norm3(config_.sun_direction);
    if (!(sn > 0.0))
        config_.sun_direction = {1.0, 0.0, 0.0};
    else
        for (double &v : config_.sun_direction)
            v /= sn;
    config_.heliocentric_distance_au =
        std::max(config_.heliocentric_distance_au, 0.1);
}

double AtmosphereModel::NeutralTemperature(double altitude_km) const
{
    const double u = (altitude_km - 115.7) / 20.14;
    return std::max(80.0, 196.95 - 64.56 * std::exp(-0.5 * u * u));
}

std::array<double, NumNeutral> AtmosphereModel::NeutralDensity(double altitude_km) const
{
    // Do not extrapolate the thermospheric fits into the solid planet.
    const double z = std::max(altitude_km, config_.lower_boundary_altitude_km);
    std::array<double, NumNeutral> n{};

    // CO2 and cold O: two-scale-height analytic fits to a Mars thermosphere
    // profile. z is in km and densities are cm^-3.
    n[CO2] = 6.04e18 * clamp_exp(-z / 6.98) +
             1.67e15 * clamp_exp(-z / 11.49);
    const double oxygen_cold = 5.85e13 * clamp_exp(-z / 10.56) +
                               7.02e9 * clamp_exp(-z / 33.97);

    // Chamberlain-like two-component cold H exosphere (z and radii in km).
    const double r_km = z + 3393.5;
    n[H] = 1.50e5 * clamp_exp(25965.0 * (1.0 / r_km - 1.0 / 3593.5)) +
           1.90e4 * clamp_exp(10365.0 * (1.0 / r_km - 1.0 / 3593.5));

    // The gravity/temperature integral is common to all minor species and is
    // evaluated once per point.
    const double z0 = config_.lower_boundary_altitude_km;
    constexpr int panels = 32;
    const double dz = std::abs(z-z0)/panels;
    double integral = 0.0;
    for(int p=0;p<=panels;++p)
    {
        const double zp=std::min(z,z0)+p*dz;
        const double rp=config_.mars_radius_m+zp*1.0e3;
        integral+=((p==0||p==panels)?0.5:1.0)/(NeutralTemperature(zp)*rp*rp);
    }
    integral*=dz*1.0e3*(z>=z0?1.0:-1.0);
    const double temperature_ratio=NeutralTemperature(z0)/NeutralTemperature(z);
    for (std::size_t s : {std::size_t(Ar), std::size_t(CO), std::size_t(C),
                          std::size_t(N2), std::size_t(N), std::size_t(NO),
                          std::size_t(O2)})
        n[s] = BaseDensity[s]*temperature_ratio*
               clamp_exp(-ParticleMass[s]*MarsGM/Boltzmann*integral);

    // Three-component hot oxygen corona (Dieval et al., converted m^-3 to
    // cm^-3). It is kept separate in column 12 and included in total O.
    const double hot_o =
        (5.23e9 * clamp_exp(-(z - 0.0) / 626.2) +
         9.76e8 * clamp_exp(-(z - 0.0) / 2790.0) +
         3.71e10 * clamp_exp(-(z - 0.0) / 88.47)) * 1.0e-6;
    n[O] = oxygen_cold + hot_o;
    n[10] = oxygen_cold;
    n[11] = hot_o;
    return n;
}

bool AtmosphereModel::IsInOpticalShadow(const std::array<double, 3> &x) const
{
    const double a = dot3(x, config_.sun_direction);
    if (a >= 0.0)
        return false;
    const double r2 = dot3(x, x);
    const double impact2 = std::max(0.0, r2 - a * a);
    return impact2 < config_.mars_radius_m * config_.mars_radius_m;
}

double AtmosphereModel::OpticalDepth(const std::array<double, 3> &x) const
{
    if (IsInOpticalShadow(x))
        return std::numeric_limits<double>::infinity();

    const double r = norm3(x);
    const double top = config_.mars_radius_m + 1500.0e3;
    if (r >= top)
        return 0.0;
    const double a = dot3(x, config_.sun_direction);
    const double disc = a * a + top * top - r * r;
    if (disc <= 0.0)
        return 0.0;
    const double s_end = -a + std::sqrt(disc);

    // Growing ray steps resolve the dense lower thermosphere without paying a
    // uniform one-kilometre cost throughout the extended corona.
    double tau = 0.0;
    double s = 0.0;
    double ds = 1.0e3;
    while (s < s_end)
    {
        const double h = std::min(ds, s_end - s);
        const double sm = s + 0.5 * h;
        std::array<double, 3> xm{};
        for (int d = 0; d < 3; ++d)
            xm[d] = x[d] + sm * config_.sun_direction[d];
        const double z = (norm3(xm) - config_.mars_radius_m) * 1.0e-3;
        const auto n = NeutralDensity(z);
        double extinction_cm = 0.0;
        for (std::size_t sp = 0; sp < NumSpecies; ++sp)
            extinction_cm += AbsorptionCrossSection[sp] * n[sp];
        tau += extinction_cm * h * 100.0; // h [m] -> [cm]
        if (tau > 80.0)
            return tau;
        s += h;
        ds = std::min(ds * 1.18, 100.0e3);
    }
    return tau;
}

AtmosphereState AtmosphereModel::Evaluate(const std::array<double, 3> &position_m) const
{
    AtmosphereState state;
    const double z = (norm3(position_m) - config_.mars_radius_m) * 1.0e-3;
    state.neutral_cm3 = NeutralDensity(z);
    const double tau = OpticalDepth(position_m);
    const double attenuation = std::isfinite(tau) ? std::exp(-std::min(tau, 700.0)) : 0.0;
    const double distance_scale = (1.524 * 1.524) /
                                  (config_.heliocentric_distance_au *
                                   config_.heliocentric_distance_au);
    for (std::size_t c = 0; c < NumPhoto; ++c)
        state.photo_cm3_s[c] = state.neutral_cm3[PhotoParent[c]] *
                               ChannelFrequency[c] * distance_scale * attenuation;
    return state;
}

PceSolution SolvePhotochemicalEquilibrium(const AtmosphereState &state,
                                           double ion_temperature_K,
                                           double electron_temperature_K)
{
    PceSolution out;
    double direct = 0.0;
    for (double value : state.photo_cm3_s)
        direct += value;
    if (direct <= 1.0e-20)
    {
        out.ion_cm3.fill(1.0e-12);
        out.converged = true;
        return out;
    }
    const double ne0 = std::sqrt(direct / 2.0e-7);
    std::array<double, NumSpecies> seed{};
    for (std::size_t s = 0; s < NumSpecies; ++s)
        seed[s] = state.photo_cm3_s[s];
    seed[C] += state.photo_cm3_s[10];
    seed[N] += state.photo_cm3_s[11];
    seed[CO] += state.photo_cm3_s[12];
    seed[O] += state.photo_cm3_s[13];
    double qsum = 0.0;
    for (double v : seed) qsum += v;
    for (std::size_t s = 0; s < NumSpecies; ++s)
        out.ion_cm3[s] = std::max(1.0e-12, ne0 * seed[s] / std::max(qsum, 1.0e-30));

    ChemistryInput in;
    in.neutral_cm3 = state.neutral_cm3;
    in.photo_cm3_s = state.photo_cm3_s;
    in.electron_temperature_K = electron_temperature_K;
    in.ion_temperature_K.fill(ion_temperature_K);

    for (int iter = 0; iter < 600; ++iter)
    {
        in.ion_cm3 = out.ion_cm3;
        const auto rates = EvaluateChemistry(in);
        double worst = 0.0;
        std::array<double, NumSpecies> next = out.ion_cm3;
        for (std::size_t s = 0; s < NumSpecies; ++s)
        {
            const double y = std::max(out.ion_cm3[s], 1.0e-30);
            const double loss_coefficient = rates.loss_cm3_s[s] / y;
            const double target = std::clamp(
                rates.production_cm3_s[s] / std::max(loss_coefficient, 1.0e-30),
                1.0e-12, 1.0e12);
            const double log_step = std::clamp(std::log(target / y), -2.0, 2.0);
            next[s] = y * std::exp(0.35 * log_step);
            const double scale = std::max({rates.production_cm3_s[s],
                                           rates.loss_cm3_s[s], 1.0e-30});
            worst = std::max(worst, std::abs(rates.production_cm3_s[s] -
                                              rates.loss_cm3_s[s]) / scale);
        }
        out.ion_cm3 = next;
        out.relative_residual = worst;
        out.iterations = iter + 1;
        if (iter > 4 && worst < 1.0e-7)
        {
            out.converged = true;
            break;
        }
    }
    return out;
}

PceSolution AtmosphereModel::WallIonDensity(const std::array<double, 3> &position_m) const
{
    const double r = norm3(position_m);
    std::array<double, 3> wall = position_m;
    const double rw = config_.mars_radius_m +
                      config_.lower_boundary_altitude_km * 1.0e3;
    const double scale = rw / std::max(r, 1.0);
    for (double &v : wall)
        v *= scale;
    return SolvePhotochemicalEquilibrium(Evaluate(wall), config_.ion_temperature_K,
                                         config_.electron_temperature_K);
}
}
