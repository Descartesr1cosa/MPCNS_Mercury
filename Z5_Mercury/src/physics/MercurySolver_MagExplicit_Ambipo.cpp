#include "MercurySolver.h"

void MercurySolver::AddAmbipolarEdgeEMF_()
{
    double a8 = ambi_coef;

    const double ne_floor = 1e-30;
    const double inv23 = M_H / M_Na; // consistent with your calc_Uplus()

    auto dd = [](double a, double b, double c,
                 double fp_i, double fm_i,
                 double fp_j, double fm_j,
                 double fp_k, double fm_k) -> double
    {
        return 0.5 * (a * (fp_i - fm_i) + b * (fp_j - fm_j) + c * (fp_k - fm_k));
    };

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);
        FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);
        FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib);

        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UNa = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &PVH = fld_->field(fid_.fid_PV_H, ib); // p in comp=3
        FieldBlock &PVN = fld_->field(fid_.fid_PV_Na, ib);

        FieldBlock &Exi = fld_->field(fid_.fid_E.xi, ib); // edge xi EMF (line integral)
        FieldBlock &Eet = fld_->field(fid_.fid_E.eta, ib);
        FieldBlock &Eze = fld_->field(fid_.fid_E.zeta, ib);

        if (!Jac.is_allocated() || !Axi.is_allocated() || !Aet.is_allocated() || !Aze.is_allocated())
            continue;
        if (!UH.is_allocated() || !UNa.is_allocated() || !PVH.is_allocated() || !PVN.is_allocated())
            continue;
        if (!Exi.is_allocated() || !Eet.is_allocated() || !Eze.is_allocated())
            continue;

        Block &blk = fld_->grd->grids(ib);

        // ---- local helpers, identical logic to AddSourceToRHS_B() ----
        auto derv_at = [&](int i, int j, int k,
                           double &ax, double &ay, double &az,
                           double &bx, double &by, double &bz,
                           double &cx, double &cy, double &cz)
        {
            const double V = std::abs(Jac(i, j, k, 0));
            if (V <= 0.0)
            {
                ax = ay = az = bx = by = bz = cx = cy = cz = 0.0;
                return;
            }

            ax = 0.5 * (Axi(i - 1, j, k, 0) + Axi(i, j, k, 0)) / V;
            ay = 0.5 * (Axi(i - 1, j, k, 1) + Axi(i, j, k, 1)) / V;
            az = 0.5 * (Axi(i - 1, j, k, 2) + Axi(i, j, k, 2)) / V;

            bx = 0.5 * (Aet(i, j - 1, k, 0) + Aet(i, j, k, 0)) / V;
            by = 0.5 * (Aet(i, j - 1, k, 1) + Aet(i, j, k, 1)) / V;
            bz = 0.5 * (Aet(i, j - 1, k, 2) + Aet(i, j, k, 2)) / V;

            cx = 0.5 * (Aze(i, j, k - 1, 0) + Aze(i, j, k, 0)) / V;
            cy = 0.5 * (Aze(i, j, k - 1, 1) + Aze(i, j, k, 1)) / V;
            cz = 0.5 * (Aze(i, j, k - 1, 2) + Aze(i, j, k, 2)) / V;
        };

        auto ne_at = [&](int i, int j, int k) -> double
        {
            const double rhoH = std::max(UH(i, j, k, 0), 0.0);
            const double rhoNa = std::max(UNa(i, j, k, 0), 0.0);
            return rhoH + rhoNa * inv23;
        };

        auto pe_at = [&](int i, int j, int k) -> double
        {
            return PVH(i, j, k, 3) + PVN(i, j, k, 3);
        };

        // G = a8*(1/ne)*grad(pe)  (physical vector at cell center)
        auto Gcell_at = [&](int i, int j, int k, double G[3])
        {
            const double ne = std::max(ne_at(i, j, k), ne_floor);
            const double invne = 1.0 / ne;

            double ax, ay, az, bx, by, bz, cx, cy, cz;
            derv_at(i, j, k, ax, ay, az, bx, by, bz, cx, cy, cz);

            const double pe_p_i = pe_at(i + 1, j, k), pe_m_i = pe_at(i - 1, j, k);
            const double pe_p_j = pe_at(i, j + 1, k), pe_m_j = pe_at(i, j - 1, k);
            const double pe_p_k = pe_at(i, j, k + 1), pe_m_k = pe_at(i, j, k - 1);

            const double dpex = dd(ax, bx, cx, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
            const double dpey = dd(ay, by, cy, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
            const double dpez = dd(az, bz, cz, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);

            G[0] = a8 * invne * dpex;
            G[1] = a8 * invne * dpey;
            G[2] = a8 * invne * dpez;
        };

        auto avg4 = [](const double a[3], const double b[3], const double c[3], const double d[3], double out[3])
        {
            out[0] = 0.25 * (a[0] + b[0] + c[0] + d[0]);
            out[1] = 0.25 * (a[1] + b[1] + c[1] + d[1]);
            out[2] = 0.25 * (a[2] + b[2] + c[2] + d[2]);
        };

        // ===== EdgeXi: use 4 surrounding cells averaged in (j,k) =====
        {
            Int3 elo = Exi.inner_lo();
            Int3 ehi = Exi.inner_hi();
            for (int i = elo.i; i < ehi.i; ++i)
                for (int j = elo.j; j < ehi.j; ++j)
                    for (int k = elo.k; k < ehi.k; ++k)
                    {
                        // cells around EdgeXi(i,j,k)  ~ (i, j-1/2, k-1/2)
                        double g00[3], g10[3], g01[3], g11[3], gavg[3];
                        Gcell_at(i, j, k, g00);
                        Gcell_at(i, j - 1, k, g10);
                        Gcell_at(i, j, k - 1, g01);
                        Gcell_at(i, j - 1, k - 1, g11);
                        avg4(g00, g10, g01, g11, gavg);

                        // EdgeXi uses dual node shift (j+1,k+1) per your edge sup rule:
                        const double x1 = blk.dual_x(i + 1, j + 1, k + 1);
                        const double y1 = blk.dual_y(i + 1, j + 1, k + 1);
                        const double z1 = blk.dual_z(i + 1, j + 1, k + 1);

                        const double x0 = blk.dual_x(i, j + 1, k + 1);
                        const double y0 = blk.dual_y(i, j + 1, k + 1);
                        const double z0 = blk.dual_z(i, j + 1, k + 1);

                        const double dlx = x1 - x0, dly = y1 - y0, dlz = z1 - z0;

                        // line integral contribution
                        Exi(i, j, k, 0) += (gavg[0] * dlx + gavg[1] * dly + gavg[2] * dlz);
                    }
        }

        // ===== EdgeEta: average in (i,k) =====
        {
            Int3 elo = Eet.inner_lo();
            Int3 ehi = Eet.inner_hi();
            for (int i = elo.i; i < ehi.i; ++i)
                for (int j = elo.j; j < ehi.j; ++j)
                    for (int k = elo.k; k < ehi.k; ++k)
                    {
                        // cells around EdgeEta(i,j,k) ~ (i-1/2, j, k-1/2)
                        double g00[3], g10[3], g01[3], g11[3], gavg[3];
                        Gcell_at(i, j, k, g00);
                        Gcell_at(i - 1, j, k, g10);
                        Gcell_at(i, j, k - 1, g01);
                        Gcell_at(i - 1, j, k - 1, g11);
                        avg4(g00, g10, g01, g11, gavg);

                        // EdgeEta sup shift (i+1,k+1)
                        const double x1 = blk.dual_x(i + 1, j + 1, k + 1);
                        const double y1 = blk.dual_y(i + 1, j + 1, k + 1);
                        const double z1 = blk.dual_z(i + 1, j + 1, k + 1);

                        const double x0 = blk.dual_x(i + 1, j, k + 1);
                        const double y0 = blk.dual_y(i + 1, j, k + 1);
                        const double z0 = blk.dual_z(i + 1, j, k + 1);

                        const double dlx = x1 - x0, dly = y1 - y0, dlz = z1 - z0;

                        Eet(i, j, k, 0) += (gavg[0] * dlx + gavg[1] * dly + gavg[2] * dlz);
                    }
        }

        // ===== EdgeZeta: average in (i,j) =====
        {
            Int3 elo = Eze.inner_lo();
            Int3 ehi = Eze.inner_hi();
            for (int i = elo.i; i < ehi.i; ++i)
                for (int j = elo.j; j < ehi.j; ++j)
                    for (int k = elo.k; k < ehi.k; ++k)
                    {
                        // cells around EdgeZeta(i,j,k) ~ (i-1/2, j-1/2, k)
                        double g00[3], g10[3], g01[3], g11[3], gavg[3];
                        Gcell_at(i, j, k, g00);
                        Gcell_at(i - 1, j, k, g10);
                        Gcell_at(i, j - 1, k, g01);
                        Gcell_at(i - 1, j - 1, k, g11);
                        avg4(g00, g10, g01, g11, gavg);

                        // EdgeZeta sup shift (i+1,j+1)
                        const double x1 = blk.dual_x(i + 1, j + 1, k + 1);
                        const double y1 = blk.dual_y(i + 1, j + 1, k + 1);
                        const double z1 = blk.dual_z(i + 1, j + 1, k + 1);

                        const double x0 = blk.dual_x(i + 1, j + 1, k);
                        const double y0 = blk.dual_y(i + 1, j + 1, k);
                        const double z0 = blk.dual_z(i + 1, j + 1, k);

                        const double dlx = x1 - x0, dly = y1 - y0, dlz = z1 - z0;

                        Eze(i, j, k, 0) += (gavg[0] * dlx + gavg[1] * dly + gavg[2] * dlz);
                    }
        }
    }
}