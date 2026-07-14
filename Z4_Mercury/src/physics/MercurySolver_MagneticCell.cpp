#include "MercurySolver.h"

void MercurySolver::calc_Bcell()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &Bcell = fld_->field(fid_.fid_Bcell, ib);
        auto &Bindcell = fld_->field(fid_.fid_Bindcell, ib);

        auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
        auto &Beta = fld_->field(fid_.fid_B.eta, ib);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);

        auto &Baddxi = fld_->field(fid_.fid_Badd.xi, ib);
        auto &Baddeta = fld_->field(fid_.fid_Badd.eta, ib);
        auto &Baddzeta = fld_->field(fid_.fid_Badd.zeta, ib);

        auto &W = fld_->field(fid_.fid_Bcell_from_Bface_w, ib);

        if (!Bcell.is_allocated() || !Bindcell.is_allocated() ||
            !Bxi.is_allocated() || !Beta.is_allocated() || !Bzeta.is_allocated())
            continue;

        Int3 lo = Bcell.inner_lo();
        Int3 hi = Bcell.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    // -------- total B = B + Badd --------
                    const double phi_tot[6] = {
                        -Bxi(i, j, k, 0) - Baddxi(i, j, k, 0),           // xi-
                        Bxi(i + 1, j, k, 0) + Baddxi(i + 1, j, k, 0),    // xi+
                        -Beta(i, j, k, 0) - Baddeta(i, j, k, 0),         // eta-
                        Beta(i, j + 1, k, 0) + Baddeta(i, j + 1, k, 0),  // eta+
                        -Bzeta(i, j, k, 0) - Baddzeta(i, j, k, 0),       // zeta-
                        Bzeta(i, j, k + 1, 0) + Baddzeta(i, j, k + 1, 0) // zeta+
                    };

                    // -------- induced B = B only --------
                    const double phi_ind[6] = {
                        -Bxi(i, j, k, 0),
                        Bxi(i + 1, j, k, 0),
                        -Beta(i, j, k, 0),
                        Beta(i, j + 1, k, 0),
                        -Bzeta(i, j, k, 0),
                        Bzeta(i, j, k + 1, 0)};

                    double Bx_tot = 0.0, By_tot = 0.0, Bz_tot = 0.0;
                    double Bx_ind = 0.0, By_ind = 0.0, Bz_ind = 0.0;

                    for (int n = 0; n < 6; ++n)
                    {
                        const double wt_x = W(i, j, k, n);
                        const double wt_y = W(i, j, k, 6 + n);
                        const double wt_z = W(i, j, k, 12 + n);

                        Bx_tot += wt_x * phi_tot[n];
                        By_tot += wt_y * phi_tot[n];
                        Bz_tot += wt_z * phi_tot[n];

                        Bx_ind += wt_x * phi_ind[n];
                        By_ind += wt_y * phi_ind[n];
                        Bz_ind += wt_z * phi_ind[n];
                    }

                    Bcell(i, j, k, 0) = Bx_tot;
                    Bcell(i, j, k, 1) = By_tot;
                    Bcell(i, j, k, 2) = Bz_tot;

                    Bindcell(i, j, k, 0) = Bx_ind;
                    Bindcell(i, j, k, 1) = By_ind;
                    Bindcell(i, j, k, 2) = Bz_ind;
                }
    }

    mercury_bound_.Sync("B_cell");

    if (topo_)
    {
        for (const auto &p : topo_->physical_patches)
        {
            if (p.bc_name != "Pole")
                continue;

            const int dir = std::abs(p.direction);
            if (dir != 1 && dir != 2)
                continue;

            const int ib = p.this_block;
            if (ib < 0 || ib >= nblock)
                continue;

            auto &Bindcell = fld_->field(fid_.fid_Bindcell, ib);
            auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);
            // auto &Baddzeta = fld_->field(fid_.fid_Badd.zeta, ib);
            auto &Aze = fld_->field(fid_.fid_metric.zeta, ib);

            if (!Bzeta.is_allocated() || !Aze.is_allocated())
                continue;

            const int sgn = (p.direction > 0) ? +1 : -1;

            auto write_zeta_flux = [&](int i, int j, int k)
            {
                const int kc = k;
                const double phi =
                    Bindcell(i, j, kc, 0) * Aze(i, j, k, 0) +
                    Bindcell(i, j, kc, 1) * Aze(i, j, k, 1) +
                    Bindcell(i, j, kc, 2) * Aze(i, j, k, 2);

                Bzeta(i, j, k, 0) = phi; //- Baddzeta(i, j, k, 0);
            };

            const Box3 &node = p.this_box_node;

            if (dir == 1)
            {
                const int icell = (sgn < 0) ? node.lo.i : (node.lo.i - 1);

                const int jface_lo = node.lo.j;
                const int jface_hi = node.hi.j - 1;
                const int kface_lo = node.lo.k;
                const int kface_hi = node.hi.k;

                for (int j = jface_lo; j < jface_hi; ++j)
                    for (int k = kface_lo; k < kface_hi; ++k)
                        write_zeta_flux(icell, j, k);
            }
            else
            {
                const int jcell = (sgn < 0) ? node.lo.j : (node.lo.j - 1);

                const int iface_lo = node.lo.i;
                const int iface_hi = node.hi.i - 1;
                const int kface_lo = node.lo.k;
                const int kface_hi = node.hi.k;

                for (int i = iface_lo; i < iface_hi; ++i)
                    for (int k = kface_lo; k < kface_hi; ++k)
                        write_zeta_flux(i, jcell, k);
            }
        }
    }

}

void MercurySolver::calc_divB()
{
    const int nblock = fld_->num_blocks();

    auto get_comp = [](const Int3 &a, int ax) -> int
    {
        if (ax == 0)
            return a.i;
        if (ax == 1)
            return a.j;
        return a.k;
    };

    auto set_comp = [](Int3 &a, int ax, int v)
    {
        if (ax == 0)
            a.i = v;
        else if (ax == 1)
            a.j = v;
        else
            a.k = v;
    };

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &divB = fld_->field(fid_.fid_divB, ib);
        auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
        auto &Beta = fld_->field(fid_.fid_B.eta, ib);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);
        auto &Jac = fld_->field(fid_.fid_Jac, ib);

        Int3 lo = divB.inner_lo();
        Int3 hi = divB.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    divB(i, j, k, 0) = (Bxi(i + 1, j, k, 0) - Bxi(i, j, k, 0) +
                                        Beta(i, j + 1, k, 0) - Beta(i, j, k, 0) +
                                        Bzeta(i, j, k + 1, 0) - Bzeta(i, j, k, 0)) /
                                       Jac(i, j, k, 0);
                }

        // Pole collapse intentionally identifies degenerate face fluxes, so the
        // adjacent first off-axis shell is not a meaningful divB diagnostic.
        for (const auto &p : topo_->physical_patches)
        {
            if (p.this_block != ib)
                continue;
            if (p.bc_name != "Pole")
                continue;

            const int ax = std::abs(p.direction) - 1;
            const bool high_side = (p.direction > 0);

            Int3 plo = lo;
            Int3 phi = hi;

            for (int d = 0; d < 3; ++d)
            {
                if (d == ax)
                    continue;

                const int tlo = std::max(get_comp(lo, d), get_comp(p.this_box_node.lo, d));
                const int thi = std::min(get_comp(hi, d), get_comp(p.this_box_node.hi, d) - 1);
                set_comp(plo, d, tlo);
                set_comp(phi, d, thi);
            }

            if (high_side)
            {
                set_comp(plo, ax, get_comp(hi, ax) - 1);
                set_comp(phi, ax, get_comp(hi, ax));
            }
            else
            {
                set_comp(plo, ax, get_comp(lo, ax));
                set_comp(phi, ax, get_comp(lo, ax) + 1);
            }

            if (!(plo.i < phi.i && plo.j < phi.j && plo.k < phi.k))
                continue;

            for (int i = plo.i; i < phi.i; ++i)
                for (int j = plo.j; j < phi.j; ++j)
                    for (int k = plo.k; k < phi.k; ++k)
                        divB(i, j, k, 0) = 0.0;
        }
    }
}

