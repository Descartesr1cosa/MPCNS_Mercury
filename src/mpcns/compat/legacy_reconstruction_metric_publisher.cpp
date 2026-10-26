#include "mpcns/compat/legacy_reconstruction_metric_publisher.hpp"

#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "mpcns/geometry/reconstruction_kernels.hpp"

#include <array>
#include <cmath>

namespace mpcns::compat
{
namespace
{
std::array<double, 3> minus(const std::array<double, 3> &a, const std::array<double, 3> &b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
std::array<double, 3> plus(const std::array<double, 3> &a, const std::array<double, 3> &b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
std::array<double, 3> scale(const std::array<double, 3> &a, double s) { return {s * a[0], s * a[1], s * a[2]}; }
} // namespace

void populate_legacy_reconstruction_cache(Field &fields, Grid &grid)
{
    Grid *grd = &grid;
    {
        // Calc_ covariant basis vectors a_xi, a_eta, a_zeta at cell centers
        auto &axi_ = fields.field("a_xi");
        auto &aeta_ = fields.field("a_eta");
        auto &azeta_ = fields.field("a_zeta");

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &axi = axi_[ib];
            auto &aeta = aeta_[ib];
            auto &azeta = azeta_[ib];

            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            Int3 lo = axi.inner_lo();
            Int3 hi = axi.inner_hi();

            auto node = [&](int i, int j, int k) -> std::array<double, 3>
            {
                return {x(i, j, k), y(i, j, k), z(i, j, k)};
            };

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const std::array<mpcns::geometry::CartesianVector, 8>
                            corners{{
                                node(i, j, k),
                                node(i + 1, j, k),
                                node(i, j + 1, k),
                                node(i + 1, j + 1, k),
                                node(i, j, k + 1),
                                node(i + 1, j, k + 1),
                                node(i, j + 1, k + 1),
                                node(i + 1, j + 1, k + 1)}};
                        const mpcns::geometry::CellCovariantBasis basis =
                            mpcns::geometry::cell_covariant_basis_from_corners(
                                corners);

                        for (int component = 0; component < 3; ++component)
                        {
                            axi(i, j, k, component) = basis.xi[component];
                            aeta(i, j, k, component) = basis.eta[component];
                            azeta(i, j, k, component) = basis.zeta[component];
                        }
                    }
        }
    }

    {
        // ===== Build cell cache: pinvGT_cell =====
        // pinvGT_cell approximates (G^T)^(-1) with Tikhonov regularization:
        //   pinv = (G G^T + reg I)^(-1) G
        // so that for covariant components w = [v·g_xi, v·g_eta, v·g_zeta],
        // the physical vector v is reconstructed by: v = pinvGT_cell * w.

        auto &pinvGTc_ = fields.field("pinvGT_cell");
        auto &axi_ = fields.field("a_xi");
        auto &aeta_ = fields.field("a_eta");
        auto &azeta_ = fields.field("a_zeta");

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &pinvGTc = pinvGTc_[ib];
            auto &axi = axi_[ib];
            auto &aeta = aeta_[ib];
            auto &azeta = azeta_[ib];

            Int3 lo = pinvGTc.inner_lo();
            Int3 hi = pinvGTc.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        std::array<double, 3> g_xi = {axi(i, j, k, 0), axi(i, j, k, 1), axi(i, j, k, 2)};
                        std::array<double, 3> g_eta = {aeta(i, j, k, 0), aeta(i, j, k, 1), aeta(i, j, k, 2)};
                        std::array<double, 3> g_ze = {azeta(i, j, k, 0), azeta(i, j, k, 1), azeta(i, j, k, 2)};

                        // G columns: g_xi, g_eta, g_zeta
                        double Gm[9] = {
                            g_xi[0], g_eta[0], g_ze[0],
                            g_xi[1], g_eta[1], g_ze[1],
                            g_xi[2], g_eta[2], g_ze[2]};

                        const auto pinvG =
                                mpcns::geometry::regularized_pseudoinverse_transpose(
                                    mpcns::geometry::matrix3_from_row_major(Gm));

                        for (int m = 0; m < 9; ++m)
                            pinvGTc(i, j, k, m) = pinvG[m];
                    }
        }
    }

    {
        // ===== Build edge cache: pinvGT_{xi,eta,ze} and pinvAT_{xi,eta,ze} =====
        auto &pinvGT_xi_ = fields.field("pinvGT_xi");
        auto &pinvAT_xi_ = fields.field("pinvAT_xi");
        auto &pinvGT_et_ = fields.field("pinvGT_eta");
        auto &pinvAT_et_ = fields.field("pinvAT_eta");
        auto &pinvGT_ze_ = fields.field("pinvGT_zeta");
        auto &pinvAT_ze_ = fields.field("pinvAT_zeta");

        auto &Axi_ = fields.field("JDxi");
        auto &Aeta_ = fields.field("JDet");
        auto &Azeta_ = fields.field("JDze");

        auto scale = [&](const std::array<double, 3> &a, double s) -> std::array<double, 3>
        {
            return {s * a[0], s * a[1], s * a[2]};
        };

        auto getvec3 = [&](FieldBlock &fb, int i, int j, int k) -> std::array<double, 3>
        {
            return {fb(i, j, k, 0), fb(i, j, k, 1), fb(i, j, k, 2)};
        };

        auto in_range = [&](const Int3 &lo, const Int3 &hi, int i, int j, int k) -> bool
        {
            return (i >= lo.i && i < hi.i &&
                    j >= lo.j && j < hi.j &&
                    k >= lo.k && k < hi.k);
        };

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            auto node = [&](int i, int j, int k) -> std::array<double, 3>
            {
                return {x(i, j, k), y(i, j, k), z(i, j, k)};
            };

            auto &A_xi = Axi_[ib];
            auto &A_et = Aeta_[ib];
            auto &A_ze = Azeta_[ib];

            const Int3 loAx = A_xi.inner_lo(), hiAx = A_xi.inner_hi();
            const Int3 loAe = A_et.inner_lo(), hiAe = A_et.inner_hi();
            const Int3 loAz = A_ze.inner_lo(), hiAz = A_ze.inner_hi();

            auto d_node_xi = [&](int ii, int jj, int kk) -> std::array<double, 3>
            {
                return scale(minus(node(ii + 1, jj, kk), node(ii - 1, jj, kk)), 0.5);
            };

            auto d_node_eta = [&](int ii, int jj, int kk) -> std::array<double, 3>
            {
                return scale(minus(node(ii, jj + 1, kk), node(ii, jj - 1, kk)), 0.5);
            };

            auto d_node_ze = [&](int ii, int jj, int kk) -> std::array<double, 3>
            {
                return scale(minus(node(ii, jj, kk + 1), node(ii, jj, kk - 1)), 0.5);
            };

            // -------------------------
            // EdgeXi: pinvGT_xi / pinvAT_xi
            // -------------------------
            {
                auto &pinvGT = pinvGT_xi_[ib];
                auto &pinvAT = pinvAT_xi_[ib];
                Int3 lo = pinvGT.inner_lo();
                Int3 hi = pinvGT.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            // // need j-1, k-1 for symmetric stencils
                            // bool okA =
                            //     in_range(loAe, hiAe, i, j, k) && in_range(loAe, hiAe, i, j, k - 1) &&
                            //     in_range(loAz, hiAz, i, j, k) && in_range(loAz, hiAz, i, j - 1, k);

                            // // 8-point xi-face average around edge
                            // for (int di : {0, 1})
                            //     for (int dj : {0, -1})
                            //         for (int dk : {0, -1})
                            //             okA = okA && in_range(loAx, hiAx, i + di, j + dj, k + dk);

                            // if (!okA)
                            // {
                            //     for (int m = 0; m < 9; ++m)
                            //     {
                            //         pinvGT(i, j, k, m) = 0.0;
                            //         pinvAT(i, j, k, m) = 0.0;
                            //     }
                            //     continue;
                            // }

                            // G columns: g_xi, g_eta, g_zeta (all at the xi-edge center)
                            auto r000 = node(i, j, k);
                            auto r100 = node(i + 1, j, k);

                            auto g_xi = minus(r100, r000);
                            auto g_eta = scale(plus(d_node_eta(i, j, k), d_node_eta(i + 1, j, k)), 0.5); // avg endpoints
                            auto g_ze = scale(plus(d_node_ze(i, j, k), d_node_ze(i + 1, j, k)), 0.5);    // avg endpoints

                            double Gm[9] = {
                                g_xi[0], g_eta[0], g_ze[0],
                                g_xi[1], g_eta[1], g_ze[1],
                                g_xi[2], g_eta[2], g_ze[2]};

                            const auto pinvG =
                                mpcns::geometry::regularized_pseudoinverse_transpose(
                                    mpcns::geometry::matrix3_from_row_major(Gm));
                            for (int m = 0; m < 9; ++m)
                                pinvGT(i, j, k, m) = pinvG[m];

                            // A columns: A_xi, A_eta, A_zeta (all co-located at the xi-edge center)
                            auto A_eta = scale(plus(getvec3(A_et, i, j, k), getvec3(A_et, i, j, k - 1)), 0.5);
                            auto A_zeta = scale(plus(getvec3(A_ze, i, j, k), getvec3(A_ze, i, j - 1, k)), 0.5);

                            std::array<double, 3> A_xi_edge = {0.0, 0.0, 0.0};
                            for (int di : {0, 1})
                                for (int dj : {0, -1})
                                    for (int dk : {0, -1})
                                        A_xi_edge = plus(A_xi_edge, getvec3(A_xi, i + di, j + dj, k + dk));
                            A_xi_edge = scale(A_xi_edge, 0.125);

                            double Am[9] = {
                                A_xi_edge[0], A_eta[0], A_zeta[0],
                                A_xi_edge[1], A_eta[1], A_zeta[1],
                                A_xi_edge[2], A_eta[2], A_zeta[2]};

                            const auto pinvA =
                                mpcns::geometry::regularized_pseudoinverse_transpose(
                                    mpcns::geometry::matrix3_from_row_major(Am));
                            for (int m = 0; m < 9; ++m)
                                pinvAT(i, j, k, m) = pinvA[m];
                        }
            }

            // -------------------------
            // EdgeEt: pinvGT_eta / pinvAT_eta
            // -------------------------
            {
                auto &pinvGT = pinvGT_et_[ib];
                auto &pinvAT = pinvAT_et_[ib];
                Int3 lo = pinvGT.inner_lo();
                Int3 hi = pinvGT.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            // // need i-1, k-1, and j+1 (edge direction)
                            // bool okA =
                            //     in_range(loAx, hiAx, i, j, k) && in_range(loAx, hiAx, i, j, k - 1) &&
                            //     in_range(loAz, hiAz, i, j, k) && in_range(loAz, hiAz, i - 1, j, k);

                            // // 8-point eta-face average around edge
                            // for (int di : {0, -1})
                            //     for (int dj : {0, 1})
                            //         for (int dk : {0, -1})
                            //             okA = okA && in_range(loAe, hiAe, i + di, j + dj, k + dk);

                            // if (!okA)
                            // {
                            //     for (int m = 0; m < 9; ++m)
                            //     {
                            //         pinvGT(i, j, k, m) = 0.0;
                            //         pinvAT(i, j, k, m) = 0.0;
                            //     }
                            //     continue;
                            // }

                            auto r000 = node(i, j, k);
                            auto r010 = node(i, j + 1, k);

                            auto g_eta = minus(r010, r000);
                            auto g_xi = scale(plus(d_node_xi(i, j, k), d_node_xi(i, j + 1, k)), 0.5);
                            auto g_ze = scale(plus(d_node_ze(i, j, k), d_node_ze(i, j + 1, k)), 0.5); // avg endpoints

                            double Gm[9] = {
                                g_xi[0], g_eta[0], g_ze[0],
                                g_xi[1], g_eta[1], g_ze[1],
                                g_xi[2], g_eta[2], g_ze[2]};

                            const auto pinvG =
                                mpcns::geometry::regularized_pseudoinverse_transpose(
                                    mpcns::geometry::matrix3_from_row_major(Gm));
                            for (int m = 0; m < 9; ++m)
                                pinvGT(i, j, k, m) = pinvG[m];

                            // A columns co-located at eta-edge center
                            auto A_xi_edge = scale(plus(getvec3(A_xi, i, j, k), getvec3(A_xi, i, j, k - 1)), 0.5);
                            auto A_zeta = scale(plus(getvec3(A_ze, i, j, k), getvec3(A_ze, i - 1, j, k)), 0.5);

                            std::array<double, 3> A_eta_edge = {0.0, 0.0, 0.0};
                            for (int di : {0, -1})
                                for (int dj : {0, 1})
                                    for (int dk : {0, -1})
                                        A_eta_edge = plus(A_eta_edge, getvec3(A_et, i + di, j + dj, k + dk));
                            A_eta_edge = scale(A_eta_edge, 0.125);

                            double Am[9] = {
                                A_xi_edge[0], A_eta_edge[0], A_zeta[0],
                                A_xi_edge[1], A_eta_edge[1], A_zeta[1],
                                A_xi_edge[2], A_eta_edge[2], A_zeta[2]};

                            const auto pinvA =
                                mpcns::geometry::regularized_pseudoinverse_transpose(
                                    mpcns::geometry::matrix3_from_row_major(Am));
                            for (int m = 0; m < 9; ++m)
                                pinvAT(i, j, k, m) = pinvA[m];
                        }
            }
            // -------------------------
            // EdgeZe: pinvGT_ze / pinvAT_ze
            // -------------------------
            {
                auto &pinvGT = pinvGT_ze_[ib];
                auto &pinvAT = pinvAT_ze_[ib];
                Int3 lo = pinvGT.inner_lo();
                Int3 hi = pinvGT.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            // // need i-1, j-1, and k+1 (edge direction)
                            // bool okA =
                            //     in_range(loAx, hiAx, i, j, k) && in_range(loAx, hiAx, i, j - 1, k) &&
                            //     in_range(loAe, hiAe, i, j, k) && in_range(loAe, hiAe, i - 1, j, k);

                            // // 8-point zeta-face average around edge
                            // for (int di : {0, -1})
                            //     for (int dj : {0, -1})
                            //         for (int dk : {0, 1})
                            //             okA = okA && in_range(loAz, hiAz, i + di, j + dj, k + dk);

                            // if (!okA)
                            // {
                            //     for (int m = 0; m < 9; ++m)
                            //     {
                            //         pinvGT(i, j, k, m) = 0.0;
                            //         pinvAT(i, j, k, m) = 0.0;
                            //     }
                            //     continue;
                            // }

                            auto r000 = node(i, j, k);
                            auto r001 = node(i, j, k + 1);

                            auto g_ze = minus(r001, r000);
                            auto g_xi = scale(plus(d_node_xi(i, j, k), d_node_xi(i, j, k + 1)), 0.5);
                            auto g_eta = scale(plus(d_node_eta(i, j, k), d_node_eta(i, j, k + 1)), 0.5);

                            double Gm[9] = {
                                g_xi[0], g_eta[0], g_ze[0],
                                g_xi[1], g_eta[1], g_ze[1],
                                g_xi[2], g_eta[2], g_ze[2]};

                            const auto pinvG =
                                mpcns::geometry::regularized_pseudoinverse_transpose(
                                    mpcns::geometry::matrix3_from_row_major(Gm));
                            for (int m = 0; m < 9; ++m)
                                pinvGT(i, j, k, m) = pinvG[m];

                            // A columns co-located at zeta-edge center
                            auto A_xi_edge = scale(plus(getvec3(A_xi, i, j, k), getvec3(A_xi, i, j - 1, k)), 0.5);
                            auto A_eta = scale(plus(getvec3(A_et, i, j, k), getvec3(A_et, i - 1, j, k)), 0.5);

                            std::array<double, 3> A_ze_edge = {0.0, 0.0, 0.0};
                            for (int di : {0, -1})
                                for (int dj : {0, -1})
                                    for (int dk : {0, 1})
                                        A_ze_edge = plus(A_ze_edge, getvec3(A_ze, i + di, j + dj, k + dk));
                            A_ze_edge = scale(A_ze_edge, 0.125);

                            double Am[9] = {
                                A_xi_edge[0], A_eta[0], A_ze_edge[0],
                                A_xi_edge[1], A_eta[1], A_ze_edge[1],
                                A_xi_edge[2], A_eta[2], A_ze_edge[2]};

                            const auto pinvA =
                                mpcns::geometry::regularized_pseudoinverse_transpose(
                                    mpcns::geometry::matrix3_from_row_major(Am));
                            for (int m = 0; m < 9; ++m)
                                pinvAT(i, j, k, m) = pinvA[m];
                        }
            }
        }
    }
}

void populate_legacy_diagonal_hodge_closure(Field &fields, Grid &grid)
{
    Grid *grd = &grid;
        auto &Jac_ = fields.field("Jac");
        auto &Hodge_star_0form_to_3form_cell_ = fields.field("Hodge_star_0form_to_3form_cell_lumped");
        auto &Hodge_star_3form_to_0form_cell_ = fields.field("Hodge_star_3form_to_0form_cell_lumped");

        auto &Hodge_star_1form_to_2form_edge_xi_ = fields.field("Hodge_star_1form_to_2form_edge_xi_lumped");
        auto &Hodge_star_1form_to_2form_edge_eta_ = fields.field("Hodge_star_1form_to_2form_edge_eta_lumped");
        auto &Hodge_star_1form_to_2form_edge_ze_ = fields.field("Hodge_star_1form_to_2form_edge_zeta_lumped");

        auto &Hodge_star_inverse_2form_to_1form_edge_xi_ = fields.field("Hodge_star_inverse_2form_to_1form_edge_xi_lumped");
        auto &Hodge_star_inverse_2form_to_1form_edge_eta_ = fields.field("Hodge_star_inverse_2form_to_1form_edge_eta_lumped");
        auto &Hodge_star_inverse_2form_to_1form_edge_ze_ = fields.field("Hodge_star_inverse_2form_to_1form_edge_zeta_lumped");

        constexpr double eps = 1e-300;

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &Jac = Jac_[ib];
            auto &Hodge_star_0form_to_3form_cell = Hodge_star_0form_to_3form_cell_[ib];
            auto &Hodge_star_3form_to_0form_cell = Hodge_star_3form_to_0form_cell_[ib];
            Int3 lo = Jac.get_lo();
            Int3 hi = Jac.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double V = Jac(i, j, k, 0);
                        Hodge_star_0form_to_3form_cell(i, j, k, 0) = V;
                        Hodge_star_3form_to_0form_cell(i, j, k, 0) = (std::fabs(V) < eps) ? 0.0 : 1.0 / V;
                    }

            auto fill_Hodge_star_1form_to_2form = [](FieldBlock &dst, FieldBlock &Hodge_star_inverse)
            {
                Int3 lo = dst.get_lo();
                Int3 hi = dst.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double Hodge_star_inverse_scale = Hodge_star_inverse(i, j, k, 0);
                            dst(i, j, k, 0) = (std::fabs(Hodge_star_inverse_scale) < 1e-300)
                                                  ? 0.0
                                                  : 1.0 / Hodge_star_inverse_scale;
                        }
            };

            fill_Hodge_star_1form_to_2form(Hodge_star_1form_to_2form_edge_xi_[ib],
                                           Hodge_star_inverse_2form_to_1form_edge_xi_[ib]);
            fill_Hodge_star_1form_to_2form(Hodge_star_1form_to_2form_edge_eta_[ib],
                                           Hodge_star_inverse_2form_to_1form_edge_eta_[ib]);
            fill_Hodge_star_1form_to_2form(Hodge_star_1form_to_2form_edge_ze_[ib],
                                           Hodge_star_inverse_2form_to_1form_edge_ze_[ib]);
        }
}

} // namespace mpcns::compat
