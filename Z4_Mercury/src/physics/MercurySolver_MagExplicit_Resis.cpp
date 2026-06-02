#include "MercurySolver.h"

#include <algorithm>
#include <cmath>

void MercurySolver::AddArtificialResistivityToEdgeEMF_()
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

void MercurySolver::AddResistiveEdgeEMF_To_(const IdTriplet &fid_Etarget)
{
    if (!resist_control.is_Mercury_resistance)
        return;

    const int nb = fld_->num_blocks();

    const double r_cut_in = 0.84;
    const double r_cut_out = 1.04;
    const double r0 = 0.84;
    const double r1 = 1.04;
    const double w = 0.02;

    auto yita0_of_r = [&](double r) -> double
    {
        if (r <= r_cut_in || r >= r_cut_out)
            return 0.0;
        return 0.5 * (std::tanh((r - r0) / w) - std::tanh((r - r1) / w));
    };

    // ----  add resistive E on solid blocks: E += invRem8 * yita0_edge * J ----
    for (int ib = 0; ib < nb; ++ib)
    {
        Block &blk = fld_->grd->grids(ib);
        // if (blk.block_name != "Solid")
        //     continue;

        auto &Exi = fld_->field(fid_Etarget.xi, ib);
        auto &Eeta = fld_->field(fid_Etarget.eta, ib);
        auto &Eze = fld_->field(fid_Etarget.zeta, ib);

        auto &dr_xi = fld_->field(fid_.Edge_dr.xi, ib);
        auto &dr_eta = fld_->field(fid_.Edge_dr.eta, ib);
        auto &dr_zeta = fld_->field(fid_.Edge_dr.zeta, ib);

        auto &Jcell = fld_->field(fid_.fid_Jcell, ib);

        if (!Jcell.is_allocated())
            continue;

        auto &x = grd_->grids(ib).x;
        auto &y = grd_->grids(ib).y;
        auto &z = grd_->grids(ib).z;

        auto get_edge_cells = [](int edge_axis, int i, int j, int k)
        {
            std::array<Int3, 4> c;

            if (edge_axis == 0)
            {
                c[0] = {i, j, k};
                c[1] = {i, j - 1, k};
                c[2] = {i, j, k - 1};
                c[3] = {i, j - 1, k - 1};
            }
            else if (edge_axis == 1)
            {
                c[0] = {i, j, k};
                c[1] = {i - 1, j, k};
                c[2] = {i, j, k - 1};
                c[3] = {i - 1, j, k - 1};
            }
            else
            {
                c[0] = {i, j, k};
                c[1] = {i - 1, j, k};
                c[2] = {i, j - 1, k};
                c[3] = {i - 1, j - 1, k};
            }

            return c;
        };

        auto add_one_edge = [&](FieldBlock &E,
                                FieldBlock &dr,
                                int edge_axis,
                                int di,
                                int dj,
                                int dk)
        {
            if (!E.is_allocated() || !dr.is_allocated())
                return;

            Int3 lo = E.inner_lo();
            Int3 hi = E.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double xm = 0.5 * (x(i, j, k) + x(i + di, j + dj, k + dk));
                        const double ym = 0.5 * (y(i, j, k) + y(i + di, j + dj, k + dk));
                        const double zm = 0.5 * (z(i, j, k) + z(i + di, j + dj, k + dk));
                        const double r = std::sqrt(xm * xm + ym * ym + zm * zm);

                        const double yita0 = yita0_of_r(r);
                        if (yita0 == 0.0)
                            continue;

                        const auto cells = get_edge_cells(edge_axis, i, j, k);
                        double Jx = 0.0;
                        double Jy = 0.0;
                        double Jz = 0.0;

                        for (int q = 0; q < 4; ++q)
                        {
                            const int ic = cells[q].i;
                            const int jc = cells[q].j;
                            const int kc = cells[q].k;
                            Jx += Jcell(ic, jc, kc, 0);
                            Jy += Jcell(ic, jc, kc, 1);
                            Jz += Jcell(ic, jc, kc, 2);
                        }

                        Jx *= 0.25;
                        Jy *= 0.25;
                        Jz *= 0.25;

                        const double Jedge =
                            Jx * dr(i, j, k, 0) +
                            Jy * dr(i, j, k, 1) +
                            Jz * dr(i, j, k, 2);

                        E(i, j, k, 0) += (inver_Rem * yita0) * Jedge;
                    }
        };

        // --- EdgeXi: use the edge midpoint of node (i,j,k) -> (i+1,j,k) ---
        {
            add_one_edge(Exi, dr_xi, 0, 1, 0, 0);
        }

        // --- EdgeEt: node (i,j,k) -> (i,j+1,k) ---
        {
            add_one_edge(Eeta, dr_eta, 1, 0, 1, 0);
        }

        // --- EdgeZe: node (i,j,k) -> (i,j,k+1) ---
        {
            add_one_edge(Eze, dr_zeta, 2, 0, 0, 1);
        }
    }
}
