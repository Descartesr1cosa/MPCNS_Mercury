#include "MercurySolver.h"

#include "mpcns/operators/reconstruction_operators.hpp"

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
        auto &weights = fld_->field(fid_.fid_Bcell_from_Bface_w, ib);

        const mpcns::state::PrimalFaceFormView induced{Bxi, Beta, Bzeta};
        const mpcns::state::PrimalFaceFormView background{Baddxi, Baddeta, Baddzeta};
        if (!mpcns::operators::WeightedCellReconstruction::sum_face_2forms_to_cell_vector(
                induced, background, weights, Bcell))
            continue;

        mpcns::operators::WeightedCellReconstruction::face_2form_to_cell_vector(
            induced, weights, Bindcell);
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

