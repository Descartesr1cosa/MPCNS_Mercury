#pragma once

#include <vector>

#include "0_basic/1_MPCNS_Parameter.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"

#include "operators/CTOperators.h"
#include "operators/Vector.h"

class DipoleField
{
private:
    struct Dipole
    {
        Vec3 r0;  // Center of Diople
        Vec3 m;   // Diople Moment (the constant coefficient is absorbed in m)
        double a; // softening length（core）
    };

    // Dioples
    std::vector<Dipole> dips;

    // Magnetic vector potential A of dipoles at location
    void A_of_dipole(Vec3 &location, Vec3 &A)
    {
        A[0] = 0.0;
        A[1] = 0.0;
        A[2] = 0.0;
        for (auto &dip : dips)
        {
            Vec3 dr = location - dip.r0;
            double factor = dr.norm2() + dip.a * dip.a;
            factor = sqrt(factor) * factor;
            Vec3 A_temp = (dip.m ^ dr) / factor;
            A += A_temp;
        }
    }

public:
    DipoleField() {};

    void load_from_param(Param *par)
    {
        dips.clear();

        const int n = par->GetInt("num_of_dipole");
        dips.reserve(n);

        const double L_ref = par->GetDou_List("REF").data["L_ref"];
        const double B_ref = par->GetDou_List("REF").data["B_ref"];
        double mu_mag = par->GetDou_List("constant").data["mu_mag"];

        // nondimensional parameter for dioples
        double inver_m_mag_ref = 4 * 3.1415926535 * B_ref * L_ref * L_ref * L_ref / mu_mag;
        inver_m_mag_ref = 1.0 / inver_m_mag_ref;

        // Build list of dioples
        for (int i = 1; i <= n; ++i)
        {
            const std::string key = "Dipole" + std::to_string(i);

            auto dl = par->GetDou_List(key);

            const double x0 = dl.data["x0"];
            const double y0 = dl.data["y0"];
            const double z0 = dl.data["z0"];
            const double mx = dl.data["mx"];
            const double my = dl.data["my"];
            const double mz = dl.data["mz"];
            const double a = dl.data["a"];

            // --- map the read parameters to Dipole structure ---
            Dipole d;

            d.r0 = Vec3(x0, y0, z0); // dimensional length (position) [m]
            d.m = Vec3(mx, my, mz);  // dimensional dipole moment with [A.m.m]
            d.a = a;                 // dimensional softening length [m]

            d.r0 /= L_ref;          // nondimensionalization
            d.a /= L_ref;           // nondimensionalization
            d.m *= inver_m_mag_ref; // nondimensionalization

            dips.push_back(d);
        }
    };

    void Build_Badd_FaceFlux(Grid *grd_, Field *fld_, Param *par_,
                             int Badd_xi_id, int Badd_eta_id, int Badd_zeta_id)
    {
        // 0) Initialized as ZERO
        //-----------------------------------------------------------------------------------------
        // Badd_ initialized as 0.0
        for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
        {
            // Badd_xi
            {
                auto &Badd_xi = fld_->field(Badd_xi_id, iblock);
                const Int3 &sub = Badd_xi.get_lo();
                const Int3 &sup = Badd_xi.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_xi(i, j, k, 0) = 0.0;
            }
            // Badd_eta
            {
                auto &Badd_eta = fld_->field(Badd_eta_id, iblock);
                const Int3 &sub = Badd_eta.get_lo();
                const Int3 &sup = Badd_eta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_eta(i, j, k, 0) = 0.0;
            }
            // Badd_zeta
            {
                auto &Badd_zeta = fld_->field(Badd_zeta_id, iblock);
                const Int3 &sub = Badd_zeta.get_lo();
                const Int3 &sup = Badd_zeta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_zeta(i, j, k, 0) = 0.0;
            }
        }

        // 1) Build edge 1-form buffer：A·dr（only existed and saved here）
        std::vector<FieldBlock> Aadd_xi;
        std::vector<FieldBlock> Aadd_eta;
        std::vector<FieldBlock> Aadd_zeta;
        Aadd_xi.resize(fld_->num_blocks());
        Aadd_eta.resize(fld_->num_blocks());
        Aadd_zeta.resize(fld_->num_blocks());
        for (int b = 0; b < fld_->num_blocks(); ++b)
        {
            Aadd_xi[b].allocate(grd_->grids(b), {"Aadd_xi", StaggerLocation::EdgeXi, 1, 1});
            Aadd_eta[b].allocate(grd_->grids(b), {"Aadd_eta", StaggerLocation::EdgeEt, 1, 1});
            Aadd_zeta[b].allocate(grd_->grids(b), {"Aadd_zeta", StaggerLocation::EdgeZe, 1, 1});
        }

        for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
        {
            // Badd_xi
            {
                auto &A_xi = Aadd_xi[iblock];
                auto &x = grd_->grids(iblock).x;
                auto &y = grd_->grids(iblock).y;
                auto &z = grd_->grids(iblock).z;

                const Int3 &sub = A_xi.inner_lo();
                const Int3 &sup = A_xi.inner_hi();
                for (int i = sub.i - 1; i <= sup.i; i++)
                    for (int j = sub.j - 1; j <= sup.j; j++)
                        for (int k = sub.k - 1; k <= sup.k; k++)
                        {
                            Vec3 x1 = {x(i, j, k), y(i, j, k), z(i, j, k)};
                            Vec3 x2 = {x(i + 1, j, k), y(i + 1, j, k), z(i + 1, j, k)};
                            Vec3 xm = 0.5 * (x1 + x2);
                            Vec3 dr = x2 - x1;
                            Vec3 A;
                            A_of_dipole(xm, A);        // A is nondimensional
                            A_xi(i, j, k, 0) = A * dr; // A·dr，to corresponding edge
                        }
            }
            // Badd_eta
            {
                auto &A_eta = Aadd_eta[iblock];
                auto &x = grd_->grids(iblock).x;
                auto &y = grd_->grids(iblock).y;
                auto &z = grd_->grids(iblock).z;

                const Int3 &sub = A_eta.inner_lo();
                const Int3 &sup = A_eta.inner_hi();
                for (int i = sub.i - 1; i <= sup.i; i++)
                    for (int j = sub.j - 1; j <= sup.j; j++)
                        for (int k = sub.k - 1; k <= sup.k; k++)
                        {
                            Vec3 x1 = {x(i, j, k), y(i, j, k), z(i, j, k)};
                            Vec3 x2 = {x(i, j + 1, k), y(i, j + 1, k), z(i, j + 1, k)};
                            Vec3 xm = 0.5 * (x1 + x2);
                            Vec3 dr = x2 - x1;
                            Vec3 A;
                            A_of_dipole(xm, A);         // A is nondimensional
                            A_eta(i, j, k, 0) = A * dr; // A·dr，to corresponding edge
                        }
            }
            // Badd_zeta
            {
                auto &A_zeta = Aadd_zeta[iblock];
                auto &x = grd_->grids(iblock).x;
                auto &y = grd_->grids(iblock).y;
                auto &z = grd_->grids(iblock).z;

                const Int3 &sub = A_zeta.inner_lo();
                const Int3 &sup = A_zeta.inner_hi();
                for (int i = sub.i - 1; i <= sup.i; i++)
                    for (int j = sub.j - 1; j <= sup.j; j++)
                        for (int k = sub.k - 1; k <= sup.k; k++)
                        {
                            Vec3 x1 = {x(i, j, k), y(i, j, k), z(i, j, k)};
                            Vec3 x2 = {x(i, j, k + 1), y(i, j, k + 1), z(i, j, k + 1)};
                            Vec3 xm = 0.5 * (x1 + x2);
                            Vec3 dr = x2 - x1;
                            Vec3 A;
                            A_of_dipole(xm, A);          // A is nondimensional
                            A_zeta(i, j, k, 0) = A * dr; // A·dr，to corresponding edge
                        }
            }
        }

        // 2) curl: edge -> face The dopoles' parts of Badd_face，added to IMF
        for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
        {
            auto &Face_xi = fld_->field(Badd_xi_id, iblk);
            auto &Face_eta = fld_->field(Badd_eta_id, iblk);
            auto &Face_zeta = fld_->field(Badd_zeta_id, iblk);

            auto &Edge_xi = Aadd_xi[iblk];
            auto &Edge_eta = Aadd_eta[iblk];
            auto &Edge_zeta = Aadd_zeta[iblk];

            // CTOperators::CurlEdgeToFace(iblk, A_xi, A_eta, A_zeta, Badd_xi, Badd_eta, Badd_zeta, /*multiper=*/1.0); //  B = curl A

            double multiper = 1.0;

            // =====================================================
            // Face_xi = multiper * [ (d/dzeta)Edge_eta - (d/deta)Edge_zeta ]
            // curl_xi = (Edge_eta(i,j,k) - Edge_eta(i,j,k+1))
            //         + (Edge_zeta(i,j+1,k) - Edge_zeta(i,j,k))
            // needs: Edge_eta(...,k+1), Edge_zeta(...,j+1,...)
            // =====================================================
            {
                Int3 lo = Face_xi.get_lo();
                Int3 hi = Face_xi.get_hi();

                const Int3 lo_eeta = Edge_eta.get_lo();
                const Int3 hi_eeta = Edge_eta.get_hi();
                const Int3 lo_eze = Edge_zeta.get_lo();
                const Int3 hi_eze = Edge_zeta.get_hi();

                lo.i = std::max(lo.i, std::max(lo_eeta.i, lo_eze.i));
                lo.j = std::max(lo.j, std::max(lo_eeta.j, lo_eze.j));
                lo.k = std::max(lo.k, std::max(lo_eeta.k, lo_eze.k));

                hi.i = std::min(hi.i, std::min(hi_eeta.i, hi_eze.i));
                hi.j = std::min(hi.j, std::min(hi_eeta.j, hi_eze.j - 1)); // j+1 on Edge_zeta
                hi.k = std::min(hi.k, std::min(hi_eeta.k - 1, hi_eze.k)); // k+1 on Edge_eta

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double curl_xi =
                                (Edge_eta(i, j, k, 0) - Edge_eta(i, j, k + 1, 0)) +
                                (Edge_zeta(i, j + 1, k, 0) - Edge_zeta(i, j, k, 0));

                            Face_xi(i, j, k, 0) = multiper * curl_xi;
                        }
            }

            // =====================================================
            // Face_eta = multiper * [ (d/dxi)Edge_zeta - (d/dzeta)Edge_xi ]
            // curl_eta = (Edge_xi(i,j,k+1) - Edge_xi(i,j,k))
            //          + (Edge_zeta(i,j,k) - Edge_zeta(i+1,j,k))
            // needs: Edge_xi(...,k+1), Edge_zeta(i+1,...)
            // =====================================================
            {
                Int3 lo = Face_eta.get_lo();
                Int3 hi = Face_eta.get_hi();

                const Int3 lo_exi = Edge_xi.get_lo();
                const Int3 hi_exi = Edge_xi.get_hi();
                const Int3 lo_eze = Edge_zeta.get_lo();
                const Int3 hi_eze = Edge_zeta.get_hi();

                lo.i = std::max(lo.i, std::max(lo_exi.i, lo_eze.i));
                lo.j = std::max(lo.j, std::max(lo_exi.j, lo_eze.j));
                lo.k = std::max(lo.k, std::max(lo_exi.k, lo_eze.k));

                hi.i = std::min(hi.i, std::min(hi_exi.i, hi_eze.i - 1)); // i+1 on Edge_zeta
                hi.j = std::min(hi.j, std::min(hi_exi.j, hi_eze.j));
                hi.k = std::min(hi.k, std::min(hi_exi.k - 1, hi_eze.k)); // k+1 on Edge_xi

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double curl_eta =
                                (Edge_xi(i, j, k + 1, 0) - Edge_xi(i, j, k, 0)) +
                                (Edge_zeta(i, j, k, 0) - Edge_zeta(i + 1, j, k, 0));

                            Face_eta(i, j, k, 0) = multiper * curl_eta;
                        }
            }

            // =====================================================
            // Face_zeta = multiper * [ (d/deta)Edge_xi - (d/dxi)Edge_eta ]
            // curl_zeta = (Edge_xi(i,j,k) - Edge_xi(i,j+1,k))
            //           + (Edge_eta(i+1,j,k) - Edge_eta(i,j,k))
            // needs: Edge_xi(...,j+1,...), Edge_eta(i+1,...)
            // =====================================================
            {
                Int3 lo = Face_zeta.get_lo();
                Int3 hi = Face_zeta.get_hi();

                const Int3 lo_exi = Edge_xi.get_lo();
                const Int3 hi_exi = Edge_xi.get_hi();
                const Int3 lo_eeta = Edge_eta.get_lo();
                const Int3 hi_eeta = Edge_eta.get_hi();

                lo.i = std::max(lo.i, std::max(lo_exi.i, lo_eeta.i));
                lo.j = std::max(lo.j, std::max(lo_exi.j, lo_eeta.j));
                lo.k = std::max(lo.k, std::max(lo_exi.k, lo_eeta.k));

                hi.i = std::min(hi.i, std::min(hi_exi.i, hi_eeta.i - 1)); // i+1 on Edge_eta
                hi.j = std::min(hi.j, std::min(hi_exi.j - 1, hi_eeta.j)); // j+1 on Edge_xi
                hi.k = std::min(hi.k, std::min(hi_exi.k, hi_eeta.k));

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double curl_zeta =
                                (Edge_xi(i, j, k, 0) - Edge_xi(i, j + 1, k, 0)) +
                                (Edge_eta(i + 1, j, k, 0) - Edge_eta(i, j, k, 0));

                            Face_zeta(i, j, k, 0) = multiper * curl_zeta;
                        }
            }
        }

        // 3) IMF can be directly added to face magnetic flux：Phi = B0 · S_face
        Vec3 B0(par_->GetDou_List("INITIAL").data["Bx"] / par_->GetDou_List("REF").data["B_ref"],
                par_->GetDou_List("INITIAL").data["By"] / par_->GetDou_List("REF").data["B_ref"],
                par_->GetDou_List("INITIAL").data["Bz"] / par_->GetDou_List("REF").data["B_ref"]);
        for (int iblock = 0; iblock < fld_->num_blocks(); iblock++)
        {
            // Badd_xi
            {
                auto &Badd_xi = fld_->field(Badd_xi_id, iblock);
                auto &Aera = fld_->field("JDxi", iblock);
                const Int3 &sub = Badd_xi.get_lo();
                const Int3 &sup = Badd_xi.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_xi(i, j, k, 0) += (B0[0] * Aera(i, j, k, 0) + B0[1] * Aera(i, j, k, 1) + B0[2] * Aera(i, j, k, 2));
            }
            // Badd_eta
            {
                auto &Badd_eta = fld_->field(Badd_eta_id, iblock);
                auto &Aera = fld_->field("JDet", iblock);
                const Int3 &sub = Badd_eta.get_lo();
                const Int3 &sup = Badd_eta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_eta(i, j, k, 0) += (B0[0] * Aera(i, j, k, 0) + B0[1] * Aera(i, j, k, 1) + B0[2] * Aera(i, j, k, 2));
            }
            // Badd_zeta
            {
                auto &Badd_zeta = fld_->field(Badd_zeta_id, iblock);
                auto &Aera = fld_->field("JDze", iblock);

                const Int3 &sub = Badd_zeta.get_lo();
                const Int3 &sup = Badd_zeta.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                            Badd_zeta(i, j, k, 0) += (B0[0] * Aera(i, j, k, 0) + B0[1] * Aera(i, j, k, 1) + B0[2] * Aera(i, j, k, 2));
            }
        }
    };
};