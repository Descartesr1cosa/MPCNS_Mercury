#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "LunarSolver.h"

// Print only the basic global state ranges used for routine run monitoring.
void LunarSolver::PrintMinMaxDiagnostics_()
{
    double rhoH_min_l = std::numeric_limits<double>::infinity();
    double rhoH_max_l = -std::numeric_limits<double>::infinity();
    double pH_min_l = std::numeric_limits<double>::infinity();
    double pH_max_l = -std::numeric_limits<double>::infinity();

    double bx_min_l = std::numeric_limits<double>::infinity();
    double bx_max_l = -std::numeric_limits<double>::infinity();
    double by_min_l = std::numeric_limits<double>::infinity();
    double by_max_l = -std::numeric_limits<double>::infinity();
    double bz_min_l = std::numeric_limits<double>::infinity();
    double bz_max_l = -std::numeric_limits<double>::infinity();
    double b2_min_l = std::numeric_limits<double>::infinity();
    double b2_max_l = -std::numeric_limits<double>::infinity();
    double divb_max_l = 0.0;

    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &PVH = fld_->field(fid_.fid_PV_H, ib);
        auto &Bcell = fld_->field(fid_.fid_Bcell, ib);
        auto &divB = fld_->field(fid_.fid_divB, ib);

        if (!UH.is_allocated() || !PVH.is_allocated() ||
            !Bcell.is_allocated() || !divB.is_allocated())
            continue;

        const Int3 lo = UH.inner_lo();
        const Int3 hi = UH.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double rhoH = UH(i, j, k, 0);
                    const double pH = PVH(i, j, k, 3);
                    const double bx = Bcell(i, j, k, 0);
                    const double by = Bcell(i, j, k, 1);
                    const double bz = Bcell(i, j, k, 2);
                    const double b2 = 0.5 * (bx * bx + by * by + bz * bz);

                    rhoH_min_l = std::min(rhoH_min_l, rhoH);
                    rhoH_max_l = std::max(rhoH_max_l, rhoH);
                    pH_min_l = std::min(pH_min_l, pH);
                    pH_max_l = std::max(pH_max_l, pH);
                    bx_min_l = std::min(bx_min_l, bx);
                    bx_max_l = std::max(bx_max_l, bx);
                    by_min_l = std::min(by_min_l, by);
                    by_max_l = std::max(by_max_l, by);
                    bz_min_l = std::min(bz_min_l, bz);
                    bz_max_l = std::max(bz_max_l, bz);
                    b2_min_l = std::min(b2_min_l, b2);
                    b2_max_l = std::max(b2_max_l, b2);
                    divb_max_l = std::max(divb_max_l, std::abs(divB(i, j, k, 0)));
                }
    }

    double mins_l[6] = {
        rhoH_min_l, pH_min_l,
        bx_min_l, by_min_l, bz_min_l, b2_min_l};
    double mins_g[6];
    double maxs_l[7] = {
        rhoH_max_l, pH_max_l,
        bx_max_l, by_max_l, bz_max_l, b2_max_l, divb_max_l};
    double maxs_g[7];
    PARALLEL::mpi_min(mins_l, mins_g, 6);
    PARALLEL::mpi_max(maxs_l, maxs_g, 7);

    if (par_->GetInt("myid") != 0)
        return;

    const double Babs_max = std::sqrt(std::max(0.0, 2.0 * maxs_g[5]));
    std::printf("\n           rhoH=[%.3e, %.3e]  pH=[%.3e, %.3e]\n",
                mins_g[0], maxs_g[0], mins_g[1], maxs_g[1]);
    std::printf("           Bx  =[%.3e, %.3e]  By=[%.3e, %.3e]  Bz=[%.3e, %.3e]\n",
                mins_g[2], maxs_g[2], mins_g[3], maxs_g[3], mins_g[4], maxs_g[4]);
    std::printf("           Pmag=[%.3e, %.3e]  |B|max=%.3e  divB_max=%.3e\n\n",
                mins_g[5], maxs_g[5], Babs_max, maxs_g[6]);
    std::fflush(stdout);
}
