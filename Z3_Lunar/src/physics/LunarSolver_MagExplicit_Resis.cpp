#include "LunarSolver.h"

#include <algorithm>
#include <cmath>

void LunarSolver::AddArtificialResistivityToEdgeEMF_()
{
    const double eta_max = arti_resist_control.eta_max;
    if (eta_max <= 0.0)
        return;

    constexpr double small = 1.0e-30;
    const double J_start = arti_resist_control.J_range_start;
    const double J_on = arti_resist_control.J_range_on;

    auto eta_art = [&](double Jmag) -> double
    {
        if (J_on <= J_start)
            return (Jmag > J_start) ? eta_max : 0.0;

        double q = (Jmag - J_start) / (J_on - J_start);
        q = std::max(0.0, std::min(1.0, q));
        const double S = q * q * (3.0 - 2.0 * q);
        return eta_max * S;
    };

    auto add_one_edge = [&](FieldBlock &E, FieldBlock &J, FieldBlock &dr)
    {
        if (!E.is_allocated() || !J.is_allocated() || !dr.is_allocated())
            return;

        const Int3 lo = E.inner_lo();
        const Int3 hi = E.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double J1form = J(i, j, k, 0);
                    const double dx = dr(i, j, k, 0);
                    const double dy = dr(i, j, k, 1);
                    const double dz = dr(i, j, k, 2);
                    const double dr_len = std::sqrt(dx * dx + dy * dy + dz * dz);
                    const double Jmag = std::abs(J1form) / std::max(dr_len, small);

                    E(i, j, k, 0) += eta_art(Jmag) * J1form;
                }
    };

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        add_one_edge(fld_->field(fid_.fid_E.xi, ib),
                     fld_->field(fid_.fid_J.xi, ib),
                     fld_->field(fid_.Edge_dr.xi, ib));

        add_one_edge(fld_->field(fid_.fid_E.eta, ib),
                     fld_->field(fid_.fid_J.eta, ib),
                     fld_->field(fid_.Edge_dr.eta, ib));

        add_one_edge(fld_->field(fid_.fid_E.zeta, ib),
                     fld_->field(fid_.fid_J.zeta, ib),
                     fld_->field(fid_.Edge_dr.zeta, ib));
    }
}
void LunarSolver::AddLocalArtificialResistivityToEdgeEMF_()
{
    if (!arti_resist_control.local_enabled || arti_resist_control.local_eta_max <= 0.0 ||
        arti_resist_control.local_r_decay <= 0.0)
        return;

    const int nb = fld_->num_blocks();

    auto setup_like_edge = [](Scalar &eta, FieldBlock &E)
    {
        if (!E.is_allocated())
            return;

        const Int3 lo = E.get_lo();
        const Int3 hi = E.get_hi();
        const int ghost = -lo.i;
        eta.SetSize(hi.i - lo.i, hi.j - lo.j, hi.k - lo.k, ghost);
        eta = 0.0;
    };

    auto edge_midpoint = [&](int ib, int axis, int i, int j, int k,
                             double &xm, double &ym, double &zm)
    {
        const int di = (axis == 0) ? 1 : 0;
        const int dj = (axis == 1) ? 1 : 0;
        const int dk = (axis == 2) ? 1 : 0;

        auto &x = grd_->grids(ib).x;
        auto &y = grd_->grids(ib).y;
        auto &z = grd_->grids(ib).z;

        xm = 0.5 * (x(i, j, k) + x(i + di, j + dj, k + dk));
        ym = 0.5 * (y(i, j, k) + y(i + di, j + dj, k + dk));
        zm = 0.5 * (z(i, j, k) + z(i + di, j + dj, k + dk));
    };

    auto compute_eta = [&](double xm, double ym, double zm) -> double
    {
        const double dx = xm - arti_resist_control.local_center[0];
        const double dy = ym - arti_resist_control.local_center[1];
        const double dz = zm - arti_resist_control.local_center[2];
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (arti_resist_control.local_r_cutoff > 0.0 &&
            r >= arti_resist_control.local_r_cutoff)
            return 0.0;

        return arti_resist_control.local_eta_max *
               std::exp(-r / arti_resist_control.local_r_decay);
    };

    if (!local_arti_eta_ready_)
    {
        local_arti_eta_xi_.resize(nb);
        local_arti_eta_eta_.resize(nb);
        local_arti_eta_ze_.resize(nb);

        for (int ib = 0; ib < nb; ++ib)
        {
            auto &Exi = fld_->field(fid_.fid_E.xi, ib);
            auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
            auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

            setup_like_edge(local_arti_eta_xi_[ib], Exi);
            setup_like_edge(local_arti_eta_eta_[ib], Eeta);
            setup_like_edge(local_arti_eta_ze_[ib], Eze);

            auto fill_one_edge = [&](FieldBlock &E, Scalar &eta, int axis)
            {
                if (!E.is_allocated())
                    return;

                const Int3 lo = E.inner_lo();
                const Int3 hi = E.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            double xm = 0.0;
                            double ym = 0.0;
                            double zm = 0.0;
                            edge_midpoint(ib, axis, i, j, k, xm, ym, zm);
                            eta(i, j, k) = compute_eta(xm, ym, zm);
                        }
            };

            fill_one_edge(Exi, local_arti_eta_xi_[ib], 0);
            fill_one_edge(Eeta, local_arti_eta_eta_[ib], 1);
            fill_one_edge(Eze, local_arti_eta_ze_[ib], 2);
        }

        local_arti_eta_ready_ = true;
    }

    auto add_one_edge = [](FieldBlock &E, FieldBlock &J, Scalar &eta)
    {
        if (!E.is_allocated() || !J.is_allocated() || eta.GetA3().empty())
            return;

        const Int3 lo = E.inner_lo();
        const Int3 hi = E.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    E(i, j, k, 0) += eta(i, j, k) * J(i, j, k, 0);
    };

    for (int ib = 0; ib < nb; ++ib)
    {
        add_one_edge(fld_->field(fid_.fid_E.xi, ib),
                     fld_->field(fid_.fid_J.xi, ib),
                     local_arti_eta_xi_[ib]);
        add_one_edge(fld_->field(fid_.fid_E.eta, ib),
                     fld_->field(fid_.fid_J.eta, ib),
                     local_arti_eta_eta_[ib]);
        add_one_edge(fld_->field(fid_.fid_E.zeta, ib),
                     fld_->field(fid_.fid_J.zeta, ib),
                     local_arti_eta_ze_[ib]);
    }
}
