#include <limits>
#include <algorithm>
#include <cstdio>
#include "MercurySolver.h"

// 统计 rho/p (H, Na) 与 Bx/By/Bz 的范围，并打印（rank0）
void MercurySolver::PrintMinMaxDiagnostics_()
{
    // --- local init ---
    double rhoH_min_l = std::numeric_limits<double>::infinity();
    double rhoH_max_l = -std::numeric_limits<double>::infinity();
    double pH_min_l = std::numeric_limits<double>::infinity();
    double pH_max_l = -std::numeric_limits<double>::infinity();

    double rhoNa_min_l = std::numeric_limits<double>::infinity();
    double rhoNa_max_l = -std::numeric_limits<double>::infinity();
    double pNa_min_l = std::numeric_limits<double>::infinity();
    double pNa_max_l = -std::numeric_limits<double>::infinity();

    double bx_min_l = std::numeric_limits<double>::infinity();
    double bx_max_l = -std::numeric_limits<double>::infinity();
    double by_min_l = std::numeric_limits<double>::infinity();
    double by_max_l = -std::numeric_limits<double>::infinity();
    double bz_min_l = std::numeric_limits<double>::infinity();
    double bz_max_l = -std::numeric_limits<double>::infinity();
    double b2_min_l = std::numeric_limits<double>::infinity();
    double b2_max_l = -std::numeric_limits<double>::infinity();
    double divb_max_l = -std::numeric_limits<double>::infinity();

    const int nblock = fld_->num_blocks();
    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &UNa = fld_->field(fid_.fid_U_Na, ib);
        auto &PVH = fld_->field(fid_.fid_PV_H, ib);
        auto &PVNa = fld_->field(fid_.fid_PV_Na, ib);
        auto &Bcel = fld_->field(fid_.fid_Bcell, ib);
        auto &divB = fld_->field(fid_.fid_divB, ib);

        if (!UH.is_allocated() || !UNa.is_allocated() || !PVH.is_allocated() || !PVNa.is_allocated() || !Bcel.is_allocated())
            continue;

        const Int3 lo = UH.inner_lo();
        const Int3 hi = UH.inner_hi();

        for (int k = lo.k; k < hi.k; ++k)
            for (int j = lo.j; j < hi.j; ++j)
                for (int i = lo.i; i < hi.i; ++i)
                {
                    const double rhoH = UH(i, j, k, 0);
                    const double pH = PVH(i, j, k, 3);

                    const double rhoNa = UNa(i, j, k, 0);
                    const double pNa = PVNa(i, j, k, 3);

                    const double bx = Bcel(i, j, k, 0);
                    const double by = Bcel(i, j, k, 1);
                    const double bz = Bcel(i, j, k, 2);
                    const double divb = fabs(divB(i, j, k, 0));

                    rhoH_min_l = std::min(rhoH_min_l, rhoH);
                    rhoH_max_l = std::max(rhoH_max_l, rhoH);
                    pH_min_l = std::min(pH_min_l, pH);
                    pH_max_l = std::max(pH_max_l, pH);

                    rhoNa_min_l = std::min(rhoNa_min_l, rhoNa);
                    rhoNa_max_l = std::max(rhoNa_max_l, rhoNa);
                    pNa_min_l = std::min(pNa_min_l, pNa);
                    pNa_max_l = std::max(pNa_max_l, pNa);

                    bx_min_l = std::min(bx_min_l, bx);
                    bx_max_l = std::max(bx_max_l, bx);
                    by_min_l = std::min(by_min_l, by);
                    by_max_l = std::max(by_max_l, by);
                    bz_min_l = std::min(bz_min_l, bz);
                    bz_max_l = std::max(bz_max_l, bz);

                    b2_min_l = std::min(b2_min_l, 0.5 * (bx * bx + by * by + bz * bz));
                    b2_max_l = std::max(b2_max_l, 0.5 * (bx * bx + by * by + bz * bz));

                    divb_max_l = std::max(divb_max_l, divb);
                }
    }

    // --- MPI reduction ---
    double mins_l[8] = {rhoH_min_l, pH_min_l, rhoNa_min_l, pNa_min_l, bx_min_l, by_min_l, bz_min_l, b2_min_l};
    double mins_g[8];

    double maxs_l[9] = {rhoH_max_l, pH_max_l, rhoNa_max_l, pNa_max_l, bx_max_l, by_max_l, bz_max_l, b2_max_l, divb_max_l};
    double maxs_g[9];

    PARALLEL::mpi_min(mins_l, mins_g, 8);
    PARALLEL::mpi_max(maxs_l, maxs_g, 9);

    const int myid = par_->GetInt("myid");
    if (myid == 0)
    {
        // -----------------------------
        // NEW: compute |B|_max from Bmag_max (=0.5|B|^2)
        // -----------------------------
        const double Bmag_max = maxs_g[7]; // 0.5|B|^2
        const double Babs_max = std::sqrt(std::max(0.0, 2.0 * Bmag_max));

        // -----------------------------
        // NEW: Omega0_s and dt*Omega_max
        // Omega0_s = (q_e * B_ref * L_ref) / (m_s * U_ref)
        // (dimensionless w.r.t. t' = t * U_ref/L_ref)
        // -----------------------------
        const auto &C = par_->GetDou_List("constant").data;
        const auto &R = par_->GetDou_List("REF").data;

        const double NA = C.at("NA");
        const double q_e = C.at("q_e");
        const double Lref = R.at("L_ref");
        const double Uref = R.at("U");
        const double Bref = R.at("B_ref");

        const double mH = par_->GetDou("mole_mass1") / NA;  // kg/particle
        const double mNa = par_->GetDou("mole_mass2") / NA; // kg/particle

        const double Omega0_H = (q_e * Bref * Lref) / (mH * Uref);
        const double Omega0_Na = (q_e * Bref * Lref) / (mNa * Uref);

        const double OmegaH_max = Omega0_H * Babs_max;
        const double OmegaNa_max = Omega0_Na * Babs_max;
        const double Omega_max = std::max(OmegaH_max, OmegaNa_max);

        const double dtOmega = dt * Omega_max;

        // 与输出对齐："[  Diag  ] " 长度是 11，所以这里用 11 个空格缩进
        // std::printf("           rhoH=[%.3e, %.3e]  pH=[%.3e, %.3e]  rhoNa=[%.3e, %.3e]  pNa=[%.3e, %.3e]\n",
        //             mins_g[0], maxs_g[0], mins_g[1], maxs_g[1], mins_g[2], maxs_g[2], mins_g[3], maxs_g[3]);
        std::printf("\n           rhoH=[%.3e, %.3e]  pH=[%.3e, %.3e]\n",
                    mins_g[0], maxs_g[0], mins_g[1], maxs_g[1]);
        std::printf("           rhoN=[%.3e, %.3e]  pN=[%.3e, %.3e]\n",
                    mins_g[2], maxs_g[2], mins_g[3], maxs_g[3]);

        std::printf("           Bx  =[%.3e, %.3e]  By=[%.3e, %.3e]  Bz   =[%.3e, %.3e]  \n",
                    mins_g[4], maxs_g[4], mins_g[5], maxs_g[5], mins_g[6], maxs_g[6]);
        std::printf("           Pmag=[%.3e, %.3e]  |B|max=%.3e         divB_max =%.3e\n\n",
                    mins_g[7], maxs_g[7], Babs_max, maxs_g[8]);

        // // -----------------------------
        // // NEW: stiffness diagnostics
        // // -----------------------------
        // std::printf("           |B|max=%.3e  Omega0_H=%.3e  Omega0_Na=%.3e\n",
        //             Babs_max, Omega0_H, Omega0_Na);
        // std::printf("           OmegaH_max=%.3e  OmegaNa_max=%.3e  dt*Omega_max=%.3e\n",
        //             OmegaH_max, OmegaNa_max, dtOmega);

        std::fflush(stdout);
    }
}
