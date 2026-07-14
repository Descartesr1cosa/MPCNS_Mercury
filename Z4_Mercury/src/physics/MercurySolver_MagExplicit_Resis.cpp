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

void MercurySolver::AddLocalArtificialResistivityToEdgeEMF_()
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

void MercurySolver::AddResistiveEdgeEMF_To_(const IdTriplet &fid_Etarget)
{
    if (!resist_control.is_Mercury_resistance)
        return;

    const int nb = fld_->num_blocks();

    // Add the radially windowed Mercury resistive EMF: E += invRem*yita0*J.
    for (int ib = 0; ib < nb; ++ib)
    {
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

                        const double yita0 = MercuryResistivityShape_(r);
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

    if (singular_edges_ && !singular_edges_->empty())
    {
        auto contribution=[&](const METRIC::SingularPhysicalEdge &edge,
                              const METRIC::WeightedIncidentEntity &inc)->double
        {
            const auto &c=inc.entity;
            FieldBlock &jc=fld_->field(fid_.fid_Jcell,c.block);
            if(!jc.is_allocated()) return 0.0;
            auto &cx=grd_->grids(c.block).dual_x; auto &cy=grd_->grids(c.block).dual_y; auto &cz=grd_->grids(c.block).dual_z;
            const double x=cx(c.i+1,c.j+1,c.k+1), y=cy(c.i+1,c.j+1,c.k+1), z=cz(c.i+1,c.j+1,c.k+1);
            const double eta=MercuryResistivityShape_(std::sqrt(x*x+y*y+z*z));
            const double jedge=jc(c.i,c.j,c.k,0)*edge.canonical_edge_vector[0]+
                               jc(c.i,c.j,c.k,1)*edge.canonical_edge_vector[1]+
                               jc(c.i,c.j,c.k,2)*edge.canonical_edge_vector[2];
            return inver_Rem*eta*jedge;
        };
        singular_edges_->assemble_cell_field_to_local_owners(*fld_,fld_->descriptor(fid_Etarget.xi).name,contribution);
        singular_edges_->assemble_cell_field_to_local_owners(*fld_,fld_->descriptor(fid_Etarget.eta).name,contribution);
        singular_edges_->assemble_cell_field_to_local_owners(*fld_,fld_->descriptor(fid_Etarget.zeta).name,contribution);
    }
}
