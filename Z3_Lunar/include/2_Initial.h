#pragma once

#include "3_field/Field.h"

#include "00_Lunar_Const.h"
#include "0_SolverFields.h"
#include "operators/Dipole.h"

class Lunar_Initial
{
public:
    double qinf[5], q_pv_inf[5], B_imf[3]; // SW: conservative quanities, primitive quantities, IMF

    DipoleField dipB;

    // Remark：Initial here will not be responsible for restart reading；restart has been done by IOModule
    void Initialization(Field *fld, SolverFields &fid)
    {
        Param *par = fld->par;

        // 1) Background state（IMF / Inflow etc.）
        build_background_state(par);

        // 2) Initialization from the beginning: U will be calculated；restart: U will be read by IOModule
        const bool if_continue = par->GetBoo("continue_calc");
        if (!if_continue)
            Common_Initial(fld, fid);

        // 3) No matter from the beginning/restarting，rebuild not writed(saved) dependent fields
        Initial_Badd(fld, fid);

        if (par->GetInt("myid") == 0)
        {
            std::cout << "********************Finish the Initial Process! !*********************\n\n\n";
            std::cout << "====================================="
                      << "===============================================================\n";
            std::cout << "\t \t Multi-Physics Coupling Numerical Simulation Solver Begins NOW!\n";
            std::cout << "====================================="
                      << "===============================================================\n\n\n";
        }
    }

    void build_background_state(Param *par)
    {
        // ---- Constants ----
        double gamma = par->GetDou_List("constant").data["gamma"];
        double NA = par->GetDou_List("constant").data["NA"];
        double R_uni = par->GetDou_List("constant").data["R_uni"];
        double q_e = par->GetDou_List("constant").data["q_e"];
        double k_Boltz = R_uni / NA;
        double mu_mag = par->GetDou_List("constant").data["mu_mag"];

        List<double> ini = par->GetDou_List("INITIAL");
        List<double> ref = par->GetDou_List("REF");

        // ---- IMF（Tesla）----
        double Bx_phy = ini.data["Bx"];
        double By_phy = ini.data["By"];
        double Bz_phy = ini.data["Bz"];

        double B_ref = ref.data["B_ref"];

        // Nondimensional IMF
        double Bx = Bx_phy / B_ref;
        double By = By_phy / B_ref;
        double Bz = Bz_phy / B_ref;

        // Fluid reference physical quantities
        double L_ref = ref.data["L_ref"];
        double U_ref = ref.data["U"];
        double n_ref = ref.data["n"];
        double T_ref = ref.data["T"];
        double Molecular_mass_ref = ref.data["Molecular_mass"];
        double rho_ref = Molecular_mass_ref / NA * n_ref;

        // ---- nondimensional primitive vars ----
        double rho0 = (ini.data["n"] / n_ref) * (ini.data["Molecular_mass"] / Molecular_mass_ref);

        double c_y = ini.data["c_y"];
        double c_z = ini.data["c_z"];
        double c_x = -std::sqrt(1.0 - c_y * c_y - c_z * c_z);
        double u0 = c_x * ini.data["U"] / U_ref;
        double v0 = c_y * ini.data["U"] / U_ref;
        double w0 = c_z * ini.data["U"] / U_ref;

        double p_ini = ini.data["n"] * k_Boltz * ini.data["T"];
        double p0 = p_ini / (rho_ref * U_ref * U_ref);

        double T0 = ini.data["T"] / T_ref;

        // Solar Wind state
        q_pv_inf[0] = u0;
        q_pv_inf[1] = v0;
        q_pv_inf[2] = w0;
        q_pv_inf[3] = p0;
        q_pv_inf[4] = T0;

        qinf[0] = rho0;
        qinf[1] = rho0 * u0;
        qinf[2] = rho0 * v0;
        qinf[3] = rho0 * w0;
        qinf[4] = 0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0) // Kinetic energy
                  + p0 / (gamma - 1.0);                      // Inertial energy

        // IMF
        B_imf[0] = Bx;
        B_imf[1] = By;
        B_imf[2] = Bz;

        dipB.load_from_param(par);
    }

public:
    Lunar_Initial() {};
    ~Lunar_Initial() = default;

private:
    void Common_Initial(Field *fld, SolverFields &fid)
    {
        Param *par = fld->par;
        if (par->GetInt("myid") == 0)
            std::cout << "---->Starting the Initial Process...\n";

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            FieldBlock &UH = fld->field(fid.fid_U_H, iblock);
            if (!UH.is_allocated())
                continue;
            const Int3 &sub = UH.get_lo();
            const Int3 &sup = UH.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        for (int32_t ll = 0; ll < UH.descriptor().ncomp; ll++)
                            UH(i, j, k, ll) = qinf[ll];
                    }
        }

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            FieldBlock &Ub_xi = fld->field(fid.fid_B.xi, iblock);
            FieldBlock &Ub_et = fld->field(fid.fid_B.eta, iblock);
            FieldBlock &Ub_ze = fld->field(fid.fid_B.zeta, iblock);
            if (!Ub_xi.is_allocated() || !Ub_et.is_allocated() || !Ub_ze.is_allocated())
                continue;
            {
                const Int3 &sub = Ub_xi.get_lo();
                const Int3 &sup = Ub_xi.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            for (int32_t ll = 0; ll < Ub_xi.descriptor().ncomp; ll++)
                                Ub_xi(i, j, k, ll) = 0.0; // induced Magnetic Field
                        }
            }

            {
                const Int3 &sub = Ub_et.get_lo();
                const Int3 &sup = Ub_et.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            for (int32_t ll = 0; ll < Ub_et.descriptor().ncomp; ll++)
                                Ub_et(i, j, k, ll) = 0.0; // induced Magnetic Field
                        }
            }

            {
                const Int3 &sub = Ub_ze.get_lo();
                const Int3 &sup = Ub_ze.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            for (int32_t ll = 0; ll < Ub_ze.descriptor().ncomp; ll++)
                                Ub_ze(i, j, k, ll) = 0.0; // induced Magnetic Field
                        }
            }
        }
    }

    // Build B_add：IMF + local Dipoles
    void Initial_Badd(Field *fld, const SolverFields &fid)
    {
        // IMF
        // Dipoles
        dipB.Build_Badd_FaceFlux(fld->grd, fld, fld->par,
                                 fid.fid_Badd.xi, fid.fid_Badd.eta, fid.fid_Badd.zeta);
    }

};
