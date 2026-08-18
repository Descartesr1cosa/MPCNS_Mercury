#include "MarsChemistry.h"

#include <algorithm>
#include <cmath>

namespace MARS
{
ChemistryRates EvaluateChemistry(const ChemistryInput &in)
{
    ChemistryRates out;
    const auto &r = in.neutral_cm3;
    const auto &y = in.ion_cm3;
    const auto &q = in.photo_cm3_s;
    double ne = 0.0;
    for (double value : y)
        ne += std::max(value, 0.0);

    const double Te = std::max(in.electron_temperature_K, 50.0);
    const double TiO = std::max(in.ion_temperature_K[O], 50.0);
    const double TiC = std::max(in.ion_temperature_K[C], 50.0);
    const double TiN = std::max(in.ion_temperature_K[N], 50.0);

    // cm^3/s. Names and values follow source_species in time_advance2.for.
    const double k4 = 1.64e-10, k5 = 9.6e-11;
    const double k6 = 1.1e-9 * std::pow(800.0 / TiO, 0.39);
    const double k7 = 6.4e-10, k8 = 5.08e-10;
    const double k60 = 7.38e-8 * std::pow(1200.0 / Te, 0.56);
    const double k10 = 3.1e-7 * std::pow(300.0 / Te, 0.5);
    const double k26 = 9.9e-10, k15 = 1.1e-9, k24 = 4.0e-11;
    const double k25 = 15.0e-10;
    const double k58 = 4.82e-6 * std::pow(Te, -0.55);
    const double k67 = 5.5e-11, k17 = 1.23e-10, k18 = 5.53e-11;
    const double k19 = 8.0e-12, k20 = 1.1e-10, k13 = 8.0e-10;
    const double k14 = 4.8e-10, k16 = 9.18e-10, k72 = 1.1e-10;
    const double k21 = 4.0e-10, k22 = 1.0e-11, k23 = 4.2e-10;
    const double k27 = 7.5e-10;
    const double k68 = 7.39e-18 * std::pow(TiC, -0.15) * std::exp(-68.0 / TiC);
    const double k28 = 3.48e-10, k29 = 5.22e-10;
    const double k30 = 1.9e-9, k31 = 1.17e-10, k32 = 2.0e-13;
    const double k33 = 1.1e-10, k34 = 7.30e-11, k35 = 1.4e-11;
    const double k36 = 4.4e-10, k37 = 1.3e-10, k38 = 9.8e-12;
    const double k39 = 5.0e-11, k40 = 4.93e-10, k41 = 5.6e-12;
    const double k42 = 6.16e-11, k43 = 2.02e-10;
    const double k69 = 9.44e-19 * std::pow(TiN, 0.24) * std::exp(-26.1 / TiN);
    const double k44 = 8.33e-11, k45 = 4.72e-10, k46 = 2.23e-10;
    const double k47 = 3.07e-10, k48 = 4.46e-11, k49 = 5.2e-11;
    const double k50 = 5.2e-11, k11 = 1.0e-10, k12 = 4.6e-10;
    const double k51 = 8.0e-13, k52 = 2.1e-11, k53 = 4.4e-11;
    const double k54 = 4.25e-9 * std::pow(TiO, -1.37) * std::exp(-28.592 / TiO);
    const double k55 = 1.1e-11, k56 = 3.1e-10, k57 = 4.6e-11;
    const double k63 = 4.0e-12 * std::pow(Te, 0.6);
    const double k59 = 5.09e-7 * std::pow(Te, -0.39);
    const double k61 = 8.52e-7 * std::pow(Te, -0.37);
    const double k62 = 1.4e-10 * std::pow(Te, -0.66);
    const double k64 = 6.46e-14 * std::pow(Te, 0.7);
    const double k65 = 1.90e-10 * std::pow(Te, -0.7);
    const double k66 = 6.28e-10 * std::pow(Te, -0.59);
    const double k70 = 4.53e-9 * std::pow(Te, 0.75);
    constexpr double k71 = 0.0;

    auto &p = out.production_cm3_s;
    auto &l = out.loss_cm3_s;

    p[Ar] = q[0] + y[N2] * r[Ar] * k32;
    l[Ar] = y[Ar] * (r[CO2] * k14 + r[CO] * k53 + r[N2] * k55 +
                     r[NO] * k56 + r[O2] * k57 + ne * k63);

    p[CO2] = q[1] + y[C] * r[CO2] * k20 + y[N2] * r[CO2] * k13 +
             y[Ar] * r[CO2] * k14 + y[CO] * r[CO2] * k15 + y[N] * r[CO2] * k16;
    l[CO2] = y[CO2] * (r[O] * (k4 + k5) + r[O2] * k67 + r[NO] * k17 +
                       r[H] * k18 + r[N] * k19 + ne * k10);

    p[CO] = q[12] + q[2] + y[C] * r[CO2] * k26 + y[C] * r[O] * k68 +
            y[Ar] * r[CO] * k53 + y[N] * r[CO] * k40 + y[N] * r[CO2] * k43 +
            y[O2] * r[C] * k49 + y[C] * r[O2] * k28 + y[N2] * r[CO] * k34;
    l[CO] = y[CO] * (r[CO2] * k15 + r[C] * k72 + r[H] * k21 + r[N] * k22 +
                     r[NO] * k23 + r[O] * k24 + r[O2] * k25 + ne * k58);

    p[C] = q[10] + q[3] + y[CO] * r[C] * k72 + y[N2] * r[C] * k33 +
           y[N] * r[CO] * k41 + y[O2] * r[C] * k50;
    l[C] = y[C] * (r[CO2] * (k20 + k26) + r[NO] * k27 + r[O] * k68 +
                   r[O2] * (k28 + k29) + ne * k66);

    p[H] = q[4] + y[CO2] * r[H] * k18 + y[CO] * r[H] * k21 + y[O] * r[H] * k7;
    l[H] = y[H] * (r[NO] * k30 + r[O] * k8 + r[O2] * k31 + r[CO2] * k71 + ne * k64);

    p[N2] = q[5] + y[N] * r[N] * k69 + y[N] * r[NO] * k44 + y[Ar] * r[N2] * k55;
    l[N2] = y[N2] * (r[CO2] * k13 + r[Ar] * k32 + r[C] * k33 + r[CO] * k34 +
                     r[N] * k35 + r[NO] * k36 + r[O] * (k37 + k38) +
                     r[O2] * k39 + ne * k59);

    p[N] = q[11] + q[6] + y[CO2] * r[N] * k19 + y[N2] * r[N] * k35;
    l[N] = y[N] * (r[CO2] * (k16 + k43) + r[CO] * (k40 + k41 + k42) +
                   r[N] * k69 + r[NO] * (k44 + k45) +
                   r[O2] * (k46 + k47 + k48) + ne * k65);

    p[NO] = q[7] + y[CO2] * r[NO] * k17 + y[CO] * r[N] * k22 +
            y[CO] * r[NO] * k23 + y[C] * r[NO] * k27 + y[H] * r[NO] * k30 +
            y[N2] * r[NO] * k36 + y[N2] * r[O] * k37 + y[N] * r[CO] * k42 +
            y[N] * r[NO] * k45 + y[N] * r[O2] * k46 + y[O2] * r[N] * k11 +
            y[O2] * r[NO] * k12 + y[O] * r[NO] * k51 + y[O] * r[N2] * k54 +
            y[Ar] * r[NO] * k56;
    l[NO] = y[NO] * ne * (k61 + k70);

    p[O] = q[13] + q[8] + y[CO2] * r[O] * k5 + y[CO] * r[O] * k24 +
           y[C] * r[O2] * k29 + y[H] * r[O] * k8 + y[N2] * r[O] * k38 +
           y[N2] * r[O2] * k48;
    l[O] = y[O] * (r[CO2] * k6 + r[H] * k7 + r[NO] * k51 + r[O2] * k52 +
                   r[N2] * k54 + ne * k62);

    p[O2] = q[9] + y[CO2] * r[O] * k4 + y[CO2] * r[O2] * k67 +
            y[CO] * r[O2] * k25 + y[H] * r[O2] * k31 + y[N2] * r[O2] * k39 +
            y[N] * r[O2] * k47 + y[O] * r[CO2] * k6 + y[O] * r[O2] * k52 +
            y[Ar] * r[O2] * k57;
    l[O2] = y[O2] * (r[C] * (k49 + k50) + r[N] * k11 + r[NO] * k12 + ne * k60);

    for (double value : q)
        out.electron_production_cm3_s += value;
    out.electron_loss_cm3_s = ne * (y[Ar] * k63 + y[CO2] * k10 + y[CO] * k58 +
                                     y[C] * k66 + y[H] * k64 + y[N2] * k59 +
                                     y[N] * k65 + y[NO] * (k61 + k70) +
                                     y[O] * k62 + y[O2] * k60);
    return out;
}
}
