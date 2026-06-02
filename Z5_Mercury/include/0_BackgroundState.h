#pragma once

#include "0_basic/1_MPCNS_Parameter.h"
#include "00_Mercury_Const.h"

#include <cmath>

struct MercuryBackgroundState
{
    double q_pv_inf[5]{};
    double q_pv_infs[5]{};
    double qinf[5]{};
    double qinfs[5]{};
    double B_imf[3]{};
    double gamma{0.0};
};

inline MercuryBackgroundState BuildMercuryBackgroundState(Param *par)
{
    MercuryBackgroundState state;

    state.gamma = par->GetDou_List("constant").data["gamma"];
    const double NA = par->GetDou_List("constant").data["NA"];
    const double R_uni = par->GetDou_List("constant").data["R_uni"];
    const double k_Boltz = R_uni / NA;

    List<double> ini = par->GetDou_List("INITIAL");
    List<double> ref = par->GetDou_List("REF");

    const double B_ref = ref.data["B_ref"];
    state.B_imf[0] = ini.data["Bx"] / B_ref;
    state.B_imf[1] = ini.data["By"] / B_ref;
    state.B_imf[2] = ini.data["Bz"] / B_ref;

    const double U_ref = ref.data["U"];
    const double n_ref = ref.data["n"];
    const double T_ref = ref.data["T"];
    const double Molecular_mass_ref = ref.data["Molecular_mass"];
    const double rho_ref = Molecular_mass_ref / NA * n_ref;

    const double rho0 = (ini.data["n"] / n_ref) *
                        (ini.data["Molecular_mass"] / Molecular_mass_ref);

    const double c_y = ini.data["c_y"];
    const double c_z = ini.data["c_z"];
    const double c_x = -std::sqrt(1.0 - c_y * c_y - c_z * c_z);
    const double u0 = c_x * ini.data["U"] / U_ref;
    const double v0 = c_y * ini.data["U"] / U_ref;
    const double w0 = c_z * ini.data["U"] / U_ref;

    const double p_ini = ini.data["n"] * k_Boltz * ini.data["T"];
    const double p0 = p_ini / (rho_ref * U_ref * U_ref);
    const double T0 = ini.data["T"] / T_ref;

    state.q_pv_inf[0] = u0;
    state.q_pv_inf[1] = v0;
    state.q_pv_inf[2] = w0;
    state.q_pv_inf[3] = p0;
    state.q_pv_inf[4] = T0;

    state.qinf[0] = rho0;
    state.qinf[1] = rho0 * u0;
    state.qinf[2] = rho0 * v0;
    state.qinf[3] = rho0 * w0;
    state.qinf[4] = 0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0) +
                    p0 / (state.gamma - 1.0);

    state.q_pv_infs[0] = 0.0;
    state.q_pv_infs[1] = 0.0;
    state.q_pv_infs[2] = 0.0;
    state.q_pv_infs[3] = p0 * rho_small / 23.0;
    state.q_pv_infs[4] = T0;

    state.qinfs[0] = rho_small * rho0;
    state.qinfs[1] = 0.0;
    state.qinfs[2] = 0.0;
    state.qinfs[3] = 0.0;
    state.qinfs[4] = state.q_pv_infs[3] / (state.gamma - 1.0) +
                     0.5 * state.qinfs[0] *
                         (state.q_pv_infs[0] * state.q_pv_infs[0] +
                          state.q_pv_infs[1] * state.q_pv_infs[1] +
                          state.q_pv_infs[2] * state.q_pv_infs[2]);

    return state;
}
