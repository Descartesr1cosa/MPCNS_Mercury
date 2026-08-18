#pragma once

#include <array>
#include <cstddef>
#include <string_view>

namespace MARS
{
constexpr std::size_t NumSpecies = 10;
// The Fortran atmosphere file carries twelve columns. Only the first ten
// enter source_species; the final two are retained for restart/post-processing.
constexpr std::size_t NumNeutral = 12;
constexpr std::size_t NumReactiveNeutral = 10;
constexpr std::size_t NumPhoto = 14;

enum Species : std::size_t
{
    Ar = 0,
    CO2,
    CO,
    C,
    H,
    N2,
    N,
    NO,
    O,
    O2
};

inline constexpr std::array<std::string_view, NumSpecies> SpeciesNames = {
    "Ar", "CO2", "CO", "C", "H", "N2", "N", "NO", "O", "O2"};

// Particle masses [kg], copied from fortran-mars/initiation.for.
inline constexpr std::array<double, NumSpecies> ParticleMass = {
    6.645e-26, 7.309e-26, 4.651e-26, 1.993e-26, 1.673e-27,
    4.651e-26, 2.326e-26, 4.983e-26, 2.658e-26, 5.316e-26};
}
