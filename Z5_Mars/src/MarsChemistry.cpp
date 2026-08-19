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

namespace
{
using SpeciesVector = std::array<double, NumSpecies>;
using SpeciesMatrix = std::array<std::array<double, NumSpecies>, NumSpecies>;

double ChemistryResidual(const ChemistryInput &base, const SpeciesVector &state,
                         const SpeciesVector &initial, double dt_seconds,
                         const SpeciesVector &density_floor,
                         SpeciesVector &residual, ChemistryRates *rates_out,
                         std::array<bool,NumSpecies> *active_out=nullptr)
{
    ChemistryInput input=base;
    input.ion_cm3=state;
    const ChemistryRates rates=EvaluateChemistry(input);
    double norm=0.0;
    for(std::size_t s=0;s<NumSpecies;++s)
    {
        const double backward_euler_residual=state[s]-initial[s]-dt_seconds*(
            rates.production_cm3_s[s]-rates.loss_cm3_s[s]);
        const bool active=state[s]<=density_floor[s]*(1.0+1.0e-8)+1.0e-30 &&
                          backward_euler_residual>0.0;
        residual[s]=active?state[s]-density_floor[s]:backward_euler_residual;
        if(active_out) (*active_out)[s]=active;
        const double scale=1.0+std::abs(initial[s])+std::abs(state[s]);
        norm=std::max(norm,std::abs(residual[s])/scale);
    }
    if(rates_out) *rates_out=rates;
    return norm;
}

bool SolveDenseSystem(SpeciesMatrix matrix, SpeciesVector rhs, SpeciesVector &solution)
{
    for(std::size_t column=0;column<NumSpecies;++column)
    {
        std::size_t pivot=column;
        for(std::size_t row=column+1;row<NumSpecies;++row)
            if(std::abs(matrix[row][column])>std::abs(matrix[pivot][column])) pivot=row;
        if(!std::isfinite(matrix[pivot][column]) || std::abs(matrix[pivot][column])<1.0e-14)
            return false;
        if(pivot!=column){std::swap(matrix[pivot],matrix[column]);std::swap(rhs[pivot],rhs[column]);}
        const double diagonal=matrix[column][column];
        for(std::size_t row=column+1;row<NumSpecies;++row)
        {
            const double factor=matrix[row][column]/diagonal;
            matrix[row][column]=0.0;
            for(std::size_t q=column+1;q<NumSpecies;++q) matrix[row][q]-=factor*matrix[column][q];
            rhs[row]-=factor*rhs[column];
        }
    }
    for(int row=static_cast<int>(NumSpecies)-1;row>=0;--row)
    {
        double value=rhs[static_cast<std::size_t>(row)];
        for(std::size_t q=static_cast<std::size_t>(row)+1;q<NumSpecies;++q)
            value-=matrix[static_cast<std::size_t>(row)][q]*solution[q];
        solution[static_cast<std::size_t>(row)]=value/matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    }
    return true;
}
}

PointImplicitChemistryResult SolvePointImplicitChemistry(
    const ChemistryInput &input, double dt_seconds,
    const std::array<double, NumSpecies> &density_floor_cm3,
    int max_iterations, double relative_tolerance)
{
    PointImplicitChemistryResult result;
    SpeciesVector initial=input.ion_cm3;
    for(std::size_t s=0;s<NumSpecies;++s)
        initial[s]=std::max(initial[s],density_floor_cm3[s]);
    result.ion_cm3=initial;
    result.rates=EvaluateChemistry(input);
    if(!(dt_seconds>0.0) || !std::isfinite(dt_seconds))
    {
        result.converged=dt_seconds==0.0;
        return result;
    }

    // Positive Patankar predictor: destruction is implicit from the outset.
    for(std::size_t s=0;s<NumSpecies;++s)
    {
        const double loss_coefficient=result.rates.loss_cm3_s[s]/std::max(initial[s],1.0e-300);
        result.ion_cm3[s]=std::max(density_floor_cm3[s],
            (initial[s]+dt_seconds*result.rates.production_cm3_s[s])/
            (1.0+dt_seconds*loss_coefficient));
    }

    SpeciesVector residual{},trial_residual{};
    for(int iteration=0;iteration<max_iterations;++iteration)
    {
        ChemistryRates rates;
        std::array<bool,NumSpecies> active{};
        const double norm=ChemistryResidual(input,result.ion_cm3,initial,dt_seconds,
                                             density_floor_cm3,residual,&rates,&active);
        result.iterations=iteration+1;
        result.rates=rates;
        if(std::isfinite(norm) && norm<=relative_tolerance)
        {
            result.converged=true;
            return result;
        }
        if(!std::isfinite(norm)) return result;

        SpeciesMatrix jacobian{};
        for(std::size_t column=0;column<NumSpecies;++column)
        {
            SpeciesVector perturbed=result.ion_cm3;
            const double step=1.0e-7*std::max(std::abs(perturbed[column]),1.0e-6);
            perturbed[column]+=step;
            SpeciesVector perturbed_residual{};
            ChemistryResidual(input,perturbed,initial,dt_seconds,density_floor_cm3,
                              perturbed_residual,nullptr);
            for(std::size_t row=0;row<NumSpecies;++row)
                jacobian[row][column]=active[row]?(row==column?1.0:0.0):
                    (perturbed_residual[row]-residual[row])/step;
        }

        SpeciesVector delta{};
        SpeciesVector rhs{};
        for(std::size_t s=0;s<NumSpecies;++s) rhs[s]=-residual[s];
        if(!SolveDenseSystem(jacobian,rhs,delta)) return result;

        double alpha=1.0;
        for(std::size_t s=0;s<NumSpecies;++s)
            if(delta[s]<0.0)
                alpha=std::min(alpha,(result.ion_cm3[s]-density_floor_cm3[s])/(-delta[s]));
        alpha=std::clamp(alpha,1.0e-6,1.0);

        bool accepted=false;
        for(int line_search=0;line_search<16;++line_search)
        {
            SpeciesVector trial=result.ion_cm3;
            for(std::size_t s=0;s<NumSpecies;++s)
                trial[s]=std::max(density_floor_cm3[s],trial[s]+alpha*delta[s]);
            const double trial_norm=ChemistryResidual(input,trial,initial,dt_seconds,
                                                       density_floor_cm3,trial_residual,nullptr);
            if(std::isfinite(trial_norm) && trial_norm<norm)
            {
                result.ion_cm3=trial;
                accepted=true;
                break;
            }
            alpha*=0.5;
        }
        if(!accepted) return result;
    }

    SpeciesVector final_residual{};
    const double final_norm=ChemistryResidual(input,result.ion_cm3,initial,dt_seconds,
                                               density_floor_cm3,final_residual,&result.rates);
    result.converged=std::isfinite(final_norm) && final_norm<=std::max(relative_tolerance,1.0e-8);
    return result;
}
}
