#include "LunarSolver.h"

void LunarSolver::AddAmbipolarEdgeEMF_()
{
    if (!ambipolar_control.enabled)
        return;

    const double pressure_coef =
        (rho_ref * U_ref) / (q_e * L_ref * B_ref * n_ref);

    auto pe_cell = [](FieldBlock &PVH, int i, int j, int k) -> double
    {
        return PVH(i, j, k, 3);
    };

    auto pe_node = [&](FieldBlock &PVH, int i, int j, int k) -> double
    {
        double sum = 0.0;
        for (int di = -1; di <= 0; ++di)
            for (int dj = -1; dj <= 0; ++dj)
                for (int dk = -1; dk <= 0; ++dk)
                    sum += pe_cell(PVH, i + di, j + dj, k + dk);
        return 0.125 * sum;
    };

    auto edge_num_xi = [&](FieldBlock &UH, int i, int j, int k) -> NumInfo
    {
        double rhoH = 0.0;

        for (int dj = -1; dj <= 0; ++dj)
            for (int dk = -1; dk <= 0; ++dk)
            {
                rhoH += UH(i, j + dj, k + dk, 0);
            }

        return Hall_Num_Limiter(0.25 * rhoH);
    };

    auto edge_num_eta = [&](FieldBlock &UH, int i, int j, int k) -> NumInfo
    {
        double rhoH = 0.0;

        for (int di = -1; di <= 0; ++di)
            for (int dk = -1; dk <= 0; ++dk)
            {
                rhoH += UH(i + di, j, k + dk, 0);
            }

        return Hall_Num_Limiter(0.25 * rhoH);
    };

    auto edge_num_zeta = [&](FieldBlock &UH, int i, int j, int k) -> NumInfo
    {
        double rhoH = 0.0;

        for (int di = -1; di <= 0; ++di)
            for (int dj = -1; dj <= 0; ++dj)
            {
                rhoH += UH(i + di, j + dj, k, 0);
            }

        return Hall_Num_Limiter(0.25 * rhoH);
    };

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &PVH = fld_->field(fid_.fid_PV_H, ib);

        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Ezeta = fld_->field(fid_.fid_E.zeta, ib);

        if (!UH.is_allocated() || !PVH.is_allocated() ||
            !Exi.is_allocated() || !Eeta.is_allocated() || !Ezeta.is_allocated())
        {
            continue;
        }

        {
            const Int3 lo = Exi.inner_lo();
            const Int3 hi = Exi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double pe0 = pe_node(PVH, i, j, k);
                        const double pe1 = pe_node(PVH, i + 1, j, k);
                        const NumInfo num = edge_num_xi(UH, i, j, k);

                        Exi(i, j, k, 0) -= pressure_coef * (pe1 - pe0) / num.ne_eff;
                    }
        }

        {
            const Int3 lo = Eeta.inner_lo();
            const Int3 hi = Eeta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double pe0 = pe_node(PVH, i, j, k);
                        const double pe1 = pe_node(PVH, i, j + 1, k);
                        const NumInfo num = edge_num_eta(UH, i, j, k);

                        Eeta(i, j, k, 0) -= pressure_coef * (pe1 - pe0) / num.ne_eff;
                    }
        }

        {
            const Int3 lo = Ezeta.inner_lo();
            const Int3 hi = Ezeta.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double pe0 = pe_node(PVH, i, j, k);
                        const double pe1 = pe_node(PVH, i, j, k + 1);
                        const NumInfo num = edge_num_zeta(UH, i, j, k);

                        Ezeta(i, j, k, 0) -= pressure_coef * (pe1 - pe0) / num.ne_eff;
                    }
        }
    }
}
