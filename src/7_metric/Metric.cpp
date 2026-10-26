#include "7_metric/Metric.h"

#include "3_field/Field.h"
#include "mpcns/compat/legacy_geometry_field_registration.hpp"
#include "mpcns/geometry/geometry_primitives.hpp"
#include "mpcns/discrete/dual_face_inverse_hodge.hpp"

#include <array>
#include <limits>

namespace METRIC
{
void register_metric_fields(Field &fields, int geomtry_ghost_)
{
    mpcns::compat::register_legacy_geometry_fields(fields, geomtry_ghost_);
}

void compute_metric_fields(Field &fields, Grid &grid)
{
    Grid *grd = &grid;

    auto minus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };
    auto plus = [&](const std::array<double, 3> &a, const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    };
    auto scale = [&](const std::array<double, 3> &a, double s) -> std::array<double, 3>
    {
        return {s * a[0], s * a[1], s * a[2]};
    };

    {
        // Calc_ JDxi  = Jac\nabla\xi = Area
        auto &data_ = fields.field("JDxi");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            // Int3 inner_range_lo = data.inner_lo();
            // Int3 inner_range_hi = data.inner_hi();
            Int3 inner_range_lo = data.get_lo();
            Int3 inner_range_hi = data.get_hi();
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        const auto area = mpcns::geometry::quadrilateral_area_vector(
                            {x(i, j, k), y(i, j, k), z(i, j, k)},
                            {x(i, j + 1, k), y(i, j + 1, k), z(i, j + 1, k)},
                            {x(i, j, k + 1), y(i, j, k + 1), z(i, j, k + 1)},
                            {x(i, j + 1, k + 1), y(i, j + 1, k + 1), z(i, j + 1, k + 1)});
                        data(i, j, k, 0) = area[0];
                        data(i, j, k, 1) = area[1];
                        data(i, j, k, 2) = area[2];
                    }
        }
    }

    {
        // Calc_ JDet  = Jac\nabla\eta = Area
        auto &data_ = fields.field("JDet");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            // Int3 inner_range_lo = data.inner_lo();
            // Int3 inner_range_hi = data.inner_hi();
            Int3 inner_range_lo = data.get_lo();
            Int3 inner_range_hi = data.get_hi();
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        const auto area = mpcns::geometry::quadrilateral_area_vector(
                            {x(i, j, k), y(i, j, k), z(i, j, k)},
                            {x(i, j, k + 1), y(i, j, k + 1), z(i, j, k + 1)},
                            {x(i + 1, j, k), y(i + 1, j, k), z(i + 1, j, k)},
                            {x(i + 1, j, k + 1), y(i + 1, j, k + 1), z(i + 1, j, k + 1)});
                        data(i, j, k, 0) = area[0];
                        data(i, j, k, 1) = area[1];
                        data(i, j, k, 2) = area[2];
                    }
        }
    }

    {
        // Calc_ JDze  = Jac\nabla\zeta = Area
        auto &data_ = fields.field("JDze");
        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &data = data_[ib];
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            // Int3 inner_range_lo = data.inner_lo();
            // Int3 inner_range_hi = data.inner_hi();
            Int3 inner_range_lo = data.get_lo();
            Int3 inner_range_hi = data.get_hi();
            for (int i = inner_range_lo.i; i < inner_range_hi.i; i++)
                for (int j = inner_range_lo.j; j < inner_range_hi.j; j++)
                    for (int k = inner_range_lo.k; k < inner_range_hi.k; k++)
                    {
                        const auto area = mpcns::geometry::quadrilateral_area_vector(
                            {x(i, j, k), y(i, j, k), z(i, j, k)},
                            {x(i + 1, j, k), y(i + 1, j, k), z(i + 1, j, k)},
                            {x(i, j + 1, k), y(i, j + 1, k), z(i, j + 1, k)},
                            {x(i + 1, j + 1, k), y(i + 1, j + 1, k), z(i + 1, j + 1, k)});
                        data(i, j, k, 0) = area[0];
                        data(i, j, k, 1) = area[1];
                        data(i, j, k, 2) = area[2];
                    }
        }
    }

    {
        // Calc_ Jac = 1 / 3* \SUM Area \cdot dr
        auto &Jac_ = fields.field("Jac");
        auto &Axi_ = fields.field("JDxi");
        auto &Aeta_ = fields.field("JDet");
        auto &Azeta_ = fields.field("JDze");

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &Jac = Jac_[ib];
            auto &A_xi = Axi_[ib];
            auto &A_et = Aeta_[ib];
            auto &A_ze = Azeta_[ib];

            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            auto &cx = grd->grids(ib).dual_x; // cell center
            auto &cy = grd->grids(ib).dual_y;
            auto &cz = grd->grids(ib).dual_z;

            // Int3 lo = Jac.inner_lo();
            // Int3 hi = Jac.inner_hi();
            Int3 lo = Jac.get_lo();
            Int3 hi = Jac.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        std::array<double, 3> xc = {cx(i + 1, j + 1, k + 1),
                                                    cy(i + 1, j + 1, k + 1),
                                                    cz(i + 1, j + 1, k + 1)};


                        // 下面要小心：对 minus 面，把 A 反号变成“对 cell 的外法向”
                        // 并且 face center 要用对应那 4 个顶点的平均

                        // --- xi- face at i ---
                        std::array<double, 3> Af_xm = {-A_xi(i, j, k, 0),
                                                       -A_xi(i, j, k, 1),
                                                       -A_xi(i, j, k, 2)};
                        std::array<double, 3> xf_xm = {
                            0.25 * (x(i, j, k) + x(i, j + 1, k) + x(i, j, k + 1) + x(i, j + 1, k + 1)),
                            0.25 * (y(i, j, k) + y(i, j + 1, k) + y(i, j, k + 1) + y(i, j + 1, k + 1)),
                            0.25 * (z(i, j, k) + z(i, j + 1, k) + z(i, j, k + 1) + z(i, j + 1, k + 1))};

                        // --- xi+ face at i+1 ---
                        std::array<double, 3> Af_xp = {A_xi(i + 1, j, k, 0),
                                                       A_xi(i + 1, j, k, 1),
                                                       A_xi(i + 1, j, k, 2)};
                        std::array<double, 3> xf_xp = {
                            0.25 * (x(i + 1, j, k) + x(i + 1, j + 1, k) + x(i + 1, j, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i + 1, j, k) + y(i + 1, j + 1, k) + y(i + 1, j, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i + 1, j, k) + z(i + 1, j + 1, k) + z(i + 1, j, k + 1) + z(i + 1, j + 1, k + 1))};

                        // --- eta- face at j ---
                        std::array<double, 3> Af_em = {-A_et(i, j, k, 0),
                                                       -A_et(i, j, k, 1),
                                                       -A_et(i, j, k, 2)};
                        std::array<double, 3> xf_em = {
                            0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j, k + 1) + x(i + 1, j, k + 1)),
                            0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j, k + 1) + y(i + 1, j, k + 1)),
                            0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j, k + 1) + z(i + 1, j, k + 1))};

                        // --- eta+ face at j+1 ---
                        std::array<double, 3> Af_ep = {A_et(i, j + 1, k, 0),
                                                       A_et(i, j + 1, k, 1),
                                                       A_et(i, j + 1, k, 2)};
                        std::array<double, 3> xf_ep = {
                            0.25 * (x(i, j + 1, k) + x(i + 1, j + 1, k) + x(i, j + 1, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i, j + 1, k) + y(i + 1, j + 1, k) + y(i, j + 1, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i, j + 1, k) + z(i + 1, j + 1, k) + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1))};

                        // --- zeta- face at k ---
                        std::array<double, 3> Af_zm = {-A_ze(i, j, k, 0),
                                                       -A_ze(i, j, k, 1),
                                                       -A_ze(i, j, k, 2)};
                        std::array<double, 3> xf_zm = {
                            0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j + 1, k) + x(i + 1, j + 1, k)),
                            0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j + 1, k) + y(i + 1, j + 1, k)),
                            0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k))};

                        // --- zeta+ face at k+1 ---
                        std::array<double, 3> Af_zp = {A_ze(i, j, k + 1, 0),
                                                       A_ze(i, j, k + 1, 1),
                                                       A_ze(i, j, k + 1, 2)};
                        std::array<double, 3> xf_zp = {
                            0.25 * (x(i, j, k + 1) + x(i + 1, j, k + 1) + x(i, j + 1, k + 1) + x(i + 1, j + 1, k + 1)),
                            0.25 * (y(i, j, k + 1) + y(i + 1, j, k + 1) + y(i, j + 1, k + 1) + y(i + 1, j + 1, k + 1)),
                            0.25 * (z(i, j, k + 1) + z(i + 1, j, k + 1) + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1))};
                        const std::array<mpcns::geometry::OrientedFaceGeometry, 6> faces = {{
                            {Af_xm, xf_xm}, {Af_xp, xf_xp},
                            {Af_em, xf_em}, {Af_ep, xf_ep},
                            {Af_zm, xf_zm}, {Af_zp, xf_zp}}};
                        Jac(i, j, k, 0) =
                            mpcns::geometry::cell_volume_from_oriented_faces(xc, faces);
                    }
        }
    }

    {
        constexpr double eps = 1e-25;

        // Pull fields
        auto &JDxi_ = fields.field("JDxi");
        auto &JDet_ = fields.field("JDet");
        auto &JDze_ = fields.field("JDze");

        auto &Area_xi_ = fields.field("Area_xi");
        auto &Area_eta_ = fields.field("Area_eta");
        auto &Area_ze_ = fields.field("Area_zeta");

        auto &dlstar_xi_ = fields.field("dlstar_xi");
        auto &dlstar_eta_ = fields.field("dlstar_eta");
        auto &dlstar_ze_ = fields.field("dlstar_zeta");

        auto &Hodge_star_2form_to_1form_face_xi_ = fields.field("Hodge_star_2form_to_1form_face_xi_lumped");
        auto &Hodge_star_2form_to_1form_face_eta_ = fields.field("Hodge_star_2form_to_1form_face_eta_lumped");
        auto &Hodge_star_2form_to_1form_face_ze_ = fields.field("Hodge_star_2form_to_1form_face_zeta_lumped");

        auto &dl_xi_ = fields.field("dl_xi");
        auto &dl_eta_ = fields.field("dl_eta");
        auto &dl_ze_ = fields.field("dl_zeta");

        auto &Sstar_xi_ = fields.field("Sstar_xi");
        auto &Sstar_eta_ = fields.field("Sstar_eta");
        auto &Sstar_ze_ = fields.field("Sstar_zeta");

        auto &Astar_xi_ = fields.field("Astar_xi");
        auto &Astar_eta_ = fields.field("Astar_eta");
        auto &Astar_ze_ = fields.field("Astar_zeta");

        auto &Hodge_star_inverse_2form_to_1form_edge_xi_ = fields.field("Hodge_star_inverse_2form_to_1form_edge_xi_lumped");
        auto &Hodge_star_inverse_2form_to_1form_edge_eta_ = fields.field("Hodge_star_inverse_2form_to_1form_edge_eta_lumped");
        auto &Hodge_star_inverse_2form_to_1form_edge_ze_ = fields.field("Hodge_star_inverse_2form_to_1form_edge_zeta_lumped");

        auto get_cellc = [&](double3D &cx, double3D &cy, double3D &cz, int ii, int jj, int kk) -> std::array<double, 3>
        {
            return {cx(ii, jj, kk), cy(ii, jj, kk), cz(ii, jj, kk)};
        };

        auto get_node = [&](double3D &x, double3D &y, double3D &z, int i, int j, int k) -> std::array<double, 3>
        {
            return {x(i, j, k), y(i, j, k), z(i, j, k)};
        };


        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            auto &cx = grd->grids(ib).dual_x; // cell centers
            auto &cy = grd->grids(ib).dual_y;
            auto &cz = grd->grids(ib).dual_z;

            // -------------------------
            // Face: Area magnitude, dlstar, Hodge_star_2form_to_1form
            // -------------------------

            // FaceXi: face at (i,j,k) uses JDxi(i,j,k,:)
            {
                auto &JDxi = JDxi_[ib];
                auto &Area = Area_xi_[ib];
                auto &dlst = dlstar_xi_[ib];
                auto &Hodge_star_2form_to_1form_face = Hodge_star_2form_to_1form_face_xi_[ib];

                Int3 lo = Area.get_lo();
                Int3 hi = Area.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            std::array<double, 3> S = {JDxi(i, j, k, 0), JDxi(i, j, k, 1), JDxi(i, j, k, 2)};
                            double Smag = mpcns::geometry::euclidean_norm(S);
                            Area(i, j, k, 0) = Smag; // |S|

                            // dual-edge length across this face: distance between adjacent cell centers
                            // FaceXi(i, j, k) separates Cell(i-1,j,k) and Cell(i,j,k)
                            // Cell(i,j,k) center is dual(i+1,j+1,k+1); so:
                            // left  cell (i-1) => dual(i,  j+1,k+1)
                            // right cell (i)   => dual(i+1,j+1,k+1)
                            std::array<double, 3> cL = get_cellc(cx, cy, cz, i, j + 1, k + 1);
                            std::array<double, 3> cR = get_cellc(cx, cy, cz, i + 1, j + 1, k + 1);
                            double Lstar = mpcns::geometry::euclidean_norm(minus(cR, cL));
                            dlst(i, j, k, 0) = Lstar; // |l*|

                            Hodge_star_2form_to_1form_face(i, j, k, 0) = (Smag < eps) ? 0.0 : Lstar / Smag;
                        }
            }

            // FaceEt
            {
                auto &JDet = JDet_[ib];
                auto &Area = Area_eta_[ib];
                auto &dlst = dlstar_eta_[ib];
                auto &Hodge_star_2form_to_1form_face = Hodge_star_2form_to_1form_face_eta_[ib];

                Int3 lo = Area.get_lo();
                Int3 hi = Area.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            std::array<double, 3> S = {JDet(i, j, k, 0), JDet(i, j, k, 1), JDet(i, j, k, 2)};
                            double Smag = mpcns::geometry::euclidean_norm(S);
                            Area(i, j, k, 0) = Smag;

                            // FaceEt(i,j,k) separates Cell(i,j-1,k) and Cell(i,j,k)
                            std::array<double, 3> cL = get_cellc(cx, cy, cz, i + 1, j, k + 1);
                            std::array<double, 3> cR = get_cellc(cx, cy, cz, i + 1, j + 1, k + 1);
                            double Lstar = mpcns::geometry::euclidean_norm(minus(cR, cL));
                            dlst(i, j, k, 0) = Lstar;

                            Hodge_star_2form_to_1form_face(i, j, k, 0) = (Smag < eps) ? 0.0 : Lstar / Smag;
                        }
            }

            // FaceZe
            {
                auto &JDze = JDze_[ib];
                auto &Area = Area_ze_[ib];
                auto &dlst = dlstar_ze_[ib];
                auto &Hodge_star_2form_to_1form_face = Hodge_star_2form_to_1form_face_ze_[ib];

                Int3 lo = Area.get_lo();
                Int3 hi = Area.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            std::array<double, 3> S = {JDze(i, j, k, 0), JDze(i, j, k, 1), JDze(i, j, k, 2)};
                            double Smag = mpcns::geometry::euclidean_norm(S);
                            Area(i, j, k, 0) = Smag;

                            // FaceZe(i,j,k) separates Cell(i,j,k-1) and Cell(i,j,k)
                            std::array<double, 3> cL = get_cellc(cx, cy, cz, i + 1, j + 1, k);
                            std::array<double, 3> cR = get_cellc(cx, cy, cz, i + 1, j + 1, k + 1);
                            double Lstar = mpcns::geometry::euclidean_norm(minus(cR, cL));
                            dlst(i, j, k, 0) = Lstar;

                            Hodge_star_2form_to_1form_face(i, j, k, 0) = (Smag < eps) ? 0.0 : Lstar / Smag;
                        }
            }

            // -------------------------
            // Edge: dl, Sstar, Astar, alpha
            // -------------------------

            // Keep the accepted axis-degeneracy treatment explicit and independent
            // from legacy Grid/Field storage.
            const mpcns::discrete::AxisSafeInverseStar1Policy
                axis_safe_inverse_star_1{};

            // EdgeXi: edge from node(i,j,k) to node(i+1,j,k)
            {
                auto &dl = dl_xi_[ib];
                auto &Sstar = Sstar_xi_[ib];
                auto &Astar = Astar_xi_[ib];
                auto &Hodge_star_inverse_2form_to_1form_edge = Hodge_star_inverse_2form_to_1form_edge_xi_[ib];

                Int3 lo = dl.get_lo();
                Int3 hi = dl.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            auto r0 = get_node(x, y, z, i, j, k);
                            auto r1 = get_node(x, y, z, i + 1, j, k);
                            double L = mpcns::geometry::euclidean_norm(minus(r1, r0));
                            dl(i, j, k, 0) = L; // |e|

                            // Dual face around this xi-edge uses 4 surrounding cell centers:
                            // cells: (i, j-1, k-1), (i, j, k-1), (i, j-1, k), (i, j, k)
                            // mapping to dual indices: cell(i,j,k) -> dual(i+1,j+1,k+1)
                            auto p00 = get_cellc(cx, cy, cz, i + 1, j, k);
                            auto p10 = get_cellc(cx, cy, cz, i + 1, j + 1, k);
                            auto p01 = get_cellc(cx, cy, cz, i + 1, j, k + 1);
                            auto p11 = get_cellc(cx, cy, cz, i + 1, j + 1, k + 1);

                            const auto Avec = mpcns::geometry::quadrilateral_area_vector(
                                p00, p10, p01, p11);
                            Sstar(i, j, k, 0) = Avec[0];
                            Sstar(i, j, k, 1) = Avec[1];
                            Sstar(i, j, k, 2) = Avec[2];

                            double Amag = mpcns::geometry::euclidean_norm(Avec);
                            Astar(i, j, k, 0) = Amag;

                            // Hodge_star_inverse_2form_to_1form_edge(i, j, k, 0) = (Amag < eps) ? 0.0 : L / Amag; // Hodge_star_inverse_2form_to_1form = |e|/|S*|
                            Hodge_star_inverse_2form_to_1form_edge(i, j, k, 0) = axis_safe_inverse_star_1.apply(L, {p00, p10, p01, p11}, Amag, eps);
                        }
            }

            // EdgeEt: edge from node(i,j,k) to node(i,j+1,k)
            {
                auto &dl = dl_eta_[ib];
                auto &Sstar = Sstar_eta_[ib];
                auto &Astar = Astar_eta_[ib];
                auto &Hodge_star_inverse_2form_to_1form_edge = Hodge_star_inverse_2form_to_1form_edge_eta_[ib];

                Int3 lo = dl.get_lo();
                Int3 hi = dl.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            auto r0 = get_node(x, y, z, i, j, k);
                            auto r1 = get_node(x, y, z, i, j + 1, k);
                            double L = mpcns::geometry::euclidean_norm(minus(r1, r0));
                            dl(i, j, k, 0) = L;

                            // surrounding cells: (i-1,j,k-1),(i,j,k-1),(i-1,j,k),(i,j,k)
                            auto p00 = get_cellc(cx, cy, cz, i, j + 1, k);
                            auto p10 = get_cellc(cx, cy, cz, i + 1, j + 1, k);
                            auto p01 = get_cellc(cx, cy, cz, i, j + 1, k + 1);
                            auto p11 = get_cellc(cx, cy, cz, i + 1, j + 1, k + 1);

                            const auto Avec = mpcns::geometry::quadrilateral_area_vector(
                                p00, p10, p01, p11);
                            Sstar(i, j, k, 0) = Avec[0];
                            Sstar(i, j, k, 1) = Avec[1];
                            Sstar(i, j, k, 2) = Avec[2];

                            double Amag = mpcns::geometry::euclidean_norm(Avec);
                            Astar(i, j, k, 0) = Amag;

                            // Hodge_star_inverse_2form_to_1form_edge(i, j, k, 0) = (Amag < eps) ? 0.0 : L / Amag;
                            Hodge_star_inverse_2form_to_1form_edge(i, j, k, 0) = axis_safe_inverse_star_1.apply(L, {p00, p10, p01, p11}, Amag, eps);
                        }
            }

            // EdgeZe: edge from node(i,j,k) to node(i,j,k+1)
            {
                auto &dl = dl_ze_[ib];
                auto &Sstar = Sstar_ze_[ib];
                auto &Astar = Astar_ze_[ib];
                auto &Hodge_star_inverse_2form_to_1form_edge = Hodge_star_inverse_2form_to_1form_edge_ze_[ib];

                Int3 lo = dl.get_lo();
                Int3 hi = dl.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            auto r0 = get_node(x, y, z, i, j, k);
                            auto r1 = get_node(x, y, z, i, j, k + 1);
                            double L = mpcns::geometry::euclidean_norm(minus(r1, r0));
                            dl(i, j, k, 0) = L;

                            // surrounding cells: (i-1,j-1,k),(i,j-1,k),(i-1,j,k),(i,j,k)
                            auto p00 = get_cellc(cx, cy, cz, i, j, k + 1);
                            auto p10 = get_cellc(cx, cy, cz, i + 1, j, k + 1);
                            auto p01 = get_cellc(cx, cy, cz, i, j + 1, k + 1);
                            auto p11 = get_cellc(cx, cy, cz, i + 1, j + 1, k + 1);

                            const auto Avec = mpcns::geometry::quadrilateral_area_vector(
                                p00, p10, p01, p11);
                            Sstar(i, j, k, 0) = Avec[0];
                            Sstar(i, j, k, 1) = Avec[1];
                            Sstar(i, j, k, 2) = Avec[2];

                            double Amag = mpcns::geometry::euclidean_norm(Avec);
                            Astar(i, j, k, 0) = Amag;

                            // Hodge_star_inverse_2form_to_1form_edge(i, j, k, 0) = (Amag < eps) ? 0.0 : L / Amag;
                            Hodge_star_inverse_2form_to_1form_edge(i, j, k, 0) = axis_safe_inverse_star_1.apply(L, {p00, p10, p01, p11}, Amag, eps);
                        }
            }
        }
    }

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
                        // 8 corner nodes of cell (i,j,k)
                        auto r000 = node(i, j, k);
                        auto r100 = node(i + 1, j, k);
                        auto r010 = node(i, j + 1, k);
                        auto r110 = node(i + 1, j + 1, k);

                        auto r001 = node(i, j, k + 1);
                        auto r101 = node(i + 1, j, k + 1);
                        auto r011 = node(i, j + 1, k + 1);
                        auto r111 = node(i + 1, j + 1, k + 1);

                        // ---- a_xi: average of 4 xi-directed edges ----
                        auto ex0 = minus(r100, r000);
                        auto ex1 = minus(r110, r010);
                        auto ex2 = minus(r101, r001);
                        auto ex3 = minus(r111, r011);
                        auto exs = plus(plus(ex0, ex1), plus(ex2, ex3));

                        axi(i, j, k, 0) = 0.25 * exs[0];
                        axi(i, j, k, 1) = 0.25 * exs[1];
                        axi(i, j, k, 2) = 0.25 * exs[2];

                        // ---- a_eta: average of 4 eta-directed edges ----
                        auto ey0 = minus(r010, r000);
                        auto ey1 = minus(r110, r100);
                        auto ey2 = minus(r011, r001);
                        auto ey3 = minus(r111, r101);
                        auto eys = plus(plus(ey0, ey1), plus(ey2, ey3));

                        aeta(i, j, k, 0) = 0.25 * eys[0];
                        aeta(i, j, k, 1) = 0.25 * eys[1];
                        aeta(i, j, k, 2) = 0.25 * eys[2];

                        // ---- a_zeta: average of 4 zeta-directed edges ----
                        auto ez0 = minus(r001, r000);
                        auto ez1 = minus(r101, r100);
                        auto ez2 = minus(r011, r010);
                        auto ez3 = minus(r111, r110);
                        auto ezs = plus(plus(ez0, ez1), plus(ez2, ez3));

                        azeta(i, j, k, 0) = 0.25 * ezs[0];
                        azeta(i, j, k, 1) = 0.25 * ezs[1];
                        azeta(i, j, k, 2) = 0.25 * ezs[2];
                    }
        }
    }

    {
        // ===== Build cell cache: pinvGT_cell =====
        // pinvGT_cell approximates (G^T)^(-1) with Tikhonov regularization:
        //   pinv = (G G^T + reg I)^(-1) G
        // so that for covariant components w = [v·g_xi, v·g_eta, v·g_zeta],
        // the physical vector v is reconstructed by: v = pinvGT_cell * w.

        // invert 3x3 matrix (row-major). return false if singular.
        auto inv3x3 = [&](const double M[9], double invM[9]) -> bool
        {
            const double a00 = M[0], a01 = M[1], a02 = M[2];
            const double a10 = M[3], a11 = M[4], a12 = M[5];
            const double a20 = M[6], a21 = M[7], a22 = M[8];

            const double c00 = a11 * a22 - a12 * a21;
            const double c01 = a02 * a21 - a01 * a22;
            const double c02 = a01 * a12 - a02 * a11;

            const double c10 = a12 * a20 - a10 * a22;
            const double c11 = a00 * a22 - a02 * a20;
            const double c12 = a02 * a10 - a00 * a12;

            const double c20 = a10 * a21 - a11 * a20;
            const double c21 = a01 * a20 - a00 * a21;
            const double c22 = a00 * a11 - a01 * a10;

            const double det = a00 * c00 + a01 * c10 + a02 * c20;
            if (std::fabs(det) < 1e-300)
                return false;

            const double invdet = 1.0 / det;

            invM[0] = c00 * invdet;
            invM[1] = c01 * invdet;
            invM[2] = c02 * invdet;
            invM[3] = c10 * invdet;
            invM[4] = c11 * invdet;
            invM[5] = c12 * invdet;
            invM[6] = c20 * invdet;
            invM[7] = c21 * invdet;
            invM[8] = c22 * invdet;
            return true;
        };

        // pinv(X^T) with Tikhonov: (X X^T + reg I)^-1 X
        auto build_pinv_transpose = [&](const double X[9], double pinv[9])
        {
            // S = X X^T (3x3)
            double S[9] = {0.0};
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                {
                    double s = 0.0;
                    for (int kk = 0; kk < 3; ++kk)
                        s += X[r * 3 + kk] * X[c * 3 + kk];
                    S[r * 3 + c] = s;
                }

            const double tr = S[0] + S[4] + S[8];
            double reg = 1e-14 * (tr + 1e-300); // scale-aware regularization

            double Sinv[9];
            bool ok = false;

            for (int it = 0; it < 8; ++it)
            {
                double Stmp[9] = {
                    S[0] + reg, S[1], S[2],
                    S[3], S[4] + reg, S[5],
                    S[6], S[7], S[8] + reg};

                if (inv3x3(Stmp, Sinv))
                {
                    ok = true;
                    break;
                }
                reg *= 10.0;
            }

            if (!ok)
            {
                for (int i = 0; i < 9; ++i)
                    pinv[i] = 0.0;
                return;
            }

            // pinv = Sinv * X
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                {
                    double s = 0.0;
                    for (int kk = 0; kk < 3; ++kk)
                        s += Sinv[r * 3 + kk] * X[kk * 3 + c];
                    pinv[r * 3 + c] = s;
                }
        };

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

                        double pinvG[9];
                        build_pinv_transpose(Gm, pinvG);

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

        // invert 3x3 matrix (row-major). return false if singular.
        auto inv3x3 = [&](const double M[9], double invM[9]) -> bool
        {
            const double a00 = M[0], a01 = M[1], a02 = M[2];
            const double a10 = M[3], a11 = M[4], a12 = M[5];
            const double a20 = M[6], a21 = M[7], a22 = M[8];

            const double c00 = a11 * a22 - a12 * a21;
            const double c01 = a02 * a21 - a01 * a22;
            const double c02 = a01 * a12 - a02 * a11;

            const double c10 = a12 * a20 - a10 * a22;
            const double c11 = a00 * a22 - a02 * a20;
            const double c12 = a02 * a10 - a00 * a12;

            const double c20 = a10 * a21 - a11 * a20;
            const double c21 = a01 * a20 - a00 * a21;
            const double c22 = a00 * a11 - a01 * a10;

            const double det = a00 * c00 + a01 * c10 + a02 * c20;
            if (std::fabs(det) < 1e-300)
                return false;

            const double invdet = 1.0 / det;

            invM[0] = c00 * invdet;
            invM[1] = c01 * invdet;
            invM[2] = c02 * invdet;
            invM[3] = c10 * invdet;
            invM[4] = c11 * invdet;
            invM[5] = c12 * invdet;
            invM[6] = c20 * invdet;
            invM[7] = c21 * invdet;
            invM[8] = c22 * invdet;
            return true;
        };

        // pinv(X^T) with Tikhonov: (X X^T + reg I)^-1 X
        auto build_pinv_transpose = [&](const double X[9], double pinv[9])
        {
            // S = X X^T (3x3)
            double S[9] = {0.0};
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                {
                    double s = 0.0;
                    for (int kk = 0; kk < 3; ++kk)
                        s += X[r * 3 + kk] * X[c * 3 + kk];
                    S[r * 3 + c] = s;
                }

            const double tr = S[0] + S[4] + S[8];
            double reg = 1e-14 * (tr + 1e-300); // scale-aware regularization

            double Sinv[9];
            bool ok = false;

            for (int it = 0; it < 8; ++it)
            {
                double Stmp[9] = {
                    S[0] + reg, S[1], S[2],
                    S[3], S[4] + reg, S[5],
                    S[6], S[7], S[8] + reg};

                if (inv3x3(Stmp, Sinv))
                {
                    ok = true;
                    break;
                }
                reg *= 10.0;
            }

            if (!ok)
            {
                for (int i = 0; i < 9; ++i)
                    pinv[i] = 0.0;
                return;
            }

            // pinv = Sinv * X
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                {
                    double s = 0.0;
                    for (int kk = 0; kk < 3; ++kk)
                        s += Sinv[r * 3 + kk] * X[kk * 3 + c];
                    pinv[r * 3 + c] = s;
                }
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

                            double pinvG[9];
                            build_pinv_transpose(Gm, pinvG);
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

                            double pinvA[9];
                            build_pinv_transpose(Am, pinvA);
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

                            double pinvG[9];
                            build_pinv_transpose(Gm, pinvG);
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

                            double pinvA[9];
                            build_pinv_transpose(Am, pinvA);
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

                            double pinvG[9];
                            build_pinv_transpose(Gm, pinvG);
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

                            double pinvA[9];
                            build_pinv_transpose(Am, pinvA);
                            for (int m = 0; m < 9; ++m)
                                pinvAT(i, j, k, m) = pinvA[m];
                        }
            }
        }
    }

    {
        auto &dr_xi_ = fields.field("dr_xi");
        auto &dr_eta_ = fields.field("dr_eta");
        auto &dr_zeta_ = fields.field("dr_zeta");

        for (int ib = 0; ib < grd->nblock; ++ib)
        {
            auto &x = grd->grids(ib).x;
            auto &y = grd->grids(ib).y;
            auto &z = grd->grids(ib).z;

            // -------------------------
            // EdgeXi: node(i,j,k) -> node(i+1,j,k)
            // -------------------------
            {
                auto &dr = dr_xi_[ib];
                Int3 lo = dr.get_lo();
                Int3 hi = dr.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            dr(i, j, k, 0) = x(i + 1, j, k) - x(i, j, k);
                            dr(i, j, k, 1) = y(i + 1, j, k) - y(i, j, k);
                            dr(i, j, k, 2) = z(i + 1, j, k) - z(i, j, k);
                        }
            }

            // -------------------------
            // EdgeEta: node(i,j,k) -> node(i,j+1,k)
            // -------------------------
            {
                auto &dr = dr_eta_[ib];
                Int3 lo = dr.get_lo();
                Int3 hi = dr.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            dr(i, j, k, 0) = x(i, j + 1, k) - x(i, j, k);
                            dr(i, j, k, 1) = y(i, j + 1, k) - y(i, j, k);
                            dr(i, j, k, 2) = z(i, j + 1, k) - z(i, j, k);
                        }
            }

            // -------------------------
            // EdgeZeta: node(i,j,k) -> node(i,j,k+1)
            // -------------------------
            {
                auto &dr = dr_zeta_[ib];
                Int3 lo = dr.get_lo();
                Int3 hi = dr.get_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            dr(i, j, k, 0) = x(i, j, k + 1) - x(i, j, k);
                            dr(i, j, k, 1) = y(i, j, k + 1) - y(i, j, k);
                            dr(i, j, k, 2) = z(i, j, k + 1) - z(i, j, k);
                        }
            }
        }
    }

    {
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

    // GCL Test
    // {
    //     auto &Jac_ = fields.field("Jac");
    //     auto &Axi_ = fields.field("JDxi");
    //     auto &Aeta_ = fields.field("JDet");
    //     auto &Azeta_ = fields.field("JDze");
    //     for (int ib = 0; ib < grd->nblock; ++ib)
    //     {
    //         auto &Jac = Jac_[ib];
    //         auto &A_xi = Axi_[ib];
    //         auto &A_et = Aeta_[ib];
    //         auto &A_ze = Azeta_[ib];
    //         Int3 lo = Jac.inner_lo();
    //         Int3 hi = Jac.inner_hi();
    //         for (int i = lo.i; i < hi.i; ++i)
    //             for (int j = lo.j; j < hi.j; ++j)
    //                 for (int k = lo.k; k < hi.k; ++k)
    //                 {
    //                     double error0 = A_xi(i + 1, j, k, 0) - A_xi(i, j, k, 0);
    //                     double error1 = A_xi(i + 1, j, k, 1) - A_xi(i, j, k, 1);
    //                     double error2 = A_xi(i + 1, j, k, 2) - A_xi(i, j, k, 2);
    //                     error0 += A_et(i, j + 1, k, 0) - A_et(i, j, k, 0);
    //                     error1 += A_et(i, j + 1, k, 1) - A_et(i, j, k, 1);
    //                     error2 += A_et(i, j + 1, k, 2) - A_et(i, j, k, 2);
    //                     error0 += A_ze(i, j, k + 1, 0) - A_ze(i, j, k, 0);
    //                     error1 += A_ze(i, j, k + 1, 1) - A_ze(i, j, k, 1);
    //                     error2 += A_ze(i, j, k + 1, 2) - A_ze(i, j, k, 2);
    //                     std::cout << sqrt(error0 * error0 + error1 * error1 + error2 * error2) << "\n";
    //                 }
    //     }
    // }
}

void build_field_geometry(Field &fields, Grid &grid, int geometry_ghost)
{
    register_metric_fields(fields, geometry_ghost);
    compute_metric_fields(fields, grid);
}

MetricDiagnostics diagnose_metric_fields(Field &fields)
{
    MetricDiagnostics diag;

    auto scan_positive = [](FieldBlock &fb, std::int64_t &counter)
    {
        Int3 lo = fb.get_lo();
        Int3 hi = fb.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (!(fb(i, j, k, 0) > 0.0))
                        ++counter;
    };

    auto scan_finite = [](FieldBlock &fb, std::int64_t &counter)
    {
        Int3 lo = fb.get_lo();
        Int3 hi = fb.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (!std::isfinite(fb(i, j, k, 0)))
                        ++counter;
    };

    auto scan_zero = [](FieldBlock &fb, std::int64_t &counter)
    {
        Int3 lo = fb.get_lo();
        Int3 hi = fb.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (fb(i, j, k, 0) == 0.0)
                        ++counter;
    };

    for (int ib = 0; ib < fields.num_blocks(); ++ib)
    {
        scan_positive(fields.field("Jac", ib), diag.jac_nonpositive);

        scan_positive(fields.field("Area_xi", ib), diag.area_nonpositive);
        scan_positive(fields.field("Area_eta", ib), diag.area_nonpositive);
        scan_positive(fields.field("Area_zeta", ib), diag.area_nonpositive);

        scan_positive(fields.field("dl_xi", ib), diag.dl_nonpositive);
        scan_positive(fields.field("dl_eta", ib), diag.dl_nonpositive);
        scan_positive(fields.field("dl_zeta", ib), diag.dl_nonpositive);

        scan_finite(fields.field("Hodge_star_inverse_2form_to_1form_edge_xi_lumped", ib),
                    diag.Hodge_star_inverse_2form_to_1form_nonfinite);
        scan_finite(fields.field("Hodge_star_inverse_2form_to_1form_edge_eta_lumped", ib),
                    diag.Hodge_star_inverse_2form_to_1form_nonfinite);
        scan_finite(fields.field("Hodge_star_inverse_2form_to_1form_edge_zeta_lumped", ib),
                    diag.Hodge_star_inverse_2form_to_1form_nonfinite);

        scan_finite(fields.field("Hodge_star_2form_to_1form_face_xi_lumped", ib),
                    diag.Hodge_star_2form_to_1form_nonfinite);
        scan_finite(fields.field("Hodge_star_2form_to_1form_face_eta_lumped", ib),
                    diag.Hodge_star_2form_to_1form_nonfinite);
        scan_finite(fields.field("Hodge_star_2form_to_1form_face_zeta_lumped", ib),
                    diag.Hodge_star_2form_to_1form_nonfinite);

        scan_zero(fields.field("Astar_xi", ib), diag.near_axis_singular);
        scan_zero(fields.field("Astar_eta", ib), diag.near_axis_singular);
        scan_zero(fields.field("Astar_zeta", ib), diag.near_axis_singular);

        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_xi_lumped", ib), diag.near_axis_capped);
        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_eta_lumped", ib), diag.near_axis_capped);
        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_zeta_lumped", ib), diag.near_axis_capped);
    }

    return diag;
}

namespace
{
bool metric_in_range(const FieldBlock &fb, int i, int j, int k)
{
    const Int3 lo = fb.get_lo();
    const Int3 hi = fb.get_hi();
    return i >= lo.i && i < hi.i &&
           j >= lo.j && j < hi.j &&
           k >= lo.k && k < hi.k;
}

void metric_zero(double out_xyz[3])
{
    out_xyz[0] = 0.0;
    out_xyz[1] = 0.0;
    out_xyz[2] = 0.0;
}

bool metric_solve3x3(double A[3][3], double b[3], double x[3])
{
    double M[3][4] = {
        {A[0][0], A[0][1], A[0][2], b[0]},
        {A[1][0], A[1][1], A[1][2], b[1]},
        {A[2][0], A[2][1], A[2][2], b[2]}};

    for (int c = 0; c < 3; ++c)
    {
        int piv = c;
        double amax = std::fabs(M[c][c]);
        for (int r = c + 1; r < 3; ++r)
        {
            const double av = std::fabs(M[r][c]);
            if (av > amax)
            {
                amax = av;
                piv = r;
            }
        }
        if (amax < 1.0e-300)
            return false;
        if (piv != c)
            for (int q = c; q < 4; ++q)
                std::swap(M[c][q], M[piv][q]);

        const double inv = 1.0 / M[c][c];
        for (int q = c; q < 4; ++q)
            M[c][q] *= inv;

        for (int r = 0; r < 3; ++r)
        {
            if (r == c)
                continue;
            const double f = M[r][c];
            for (int q = c; q < 4; ++q)
                M[r][q] -= f * M[c][q];
        }
    }

    x[0] = M[0][3];
    x[1] = M[1][3];
    x[2] = M[2][3];
    return true;
}
}

bool reconstruct_face_2form_to_cell(Field &fields,
                                    int fid_xi, int fid_eta, int fid_zeta,
                                    int iblock, int i, int j, int k,
                                    double out_xyz[3])
{
    metric_zero(out_xyz);

    FieldBlock &Fxi = fields.field(fid_xi, iblock);
    FieldBlock &Feta = fields.field(fid_eta, iblock);
    FieldBlock &Fzeta = fields.field(fid_zeta, iblock);

    if (!Fxi.is_allocated() || !Feta.is_allocated() || !Fzeta.is_allocated())
        return false;

    if (!metric_in_range(Fxi, i, j, k) ||
        !metric_in_range(Fxi, i + 1, j, k) ||
        !metric_in_range(Feta, i, j, k) ||
        !metric_in_range(Feta, i, j + 1, k) ||
        !metric_in_range(Fzeta, i, j, k) ||
        !metric_in_range(Fzeta, i, j, k + 1))
        return false;

    const double phi[6] = {
        -Fxi(i, j, k, 0),
        Fxi(i + 1, j, k, 0),
        -Feta(i, j, k, 0),
        Feta(i, j + 1, k, 0),
        -Fzeta(i, j, k, 0),
        Fzeta(i, j, k + 1, 0)};

    if (fields.has_field("Bcell_from_Bface_w"))
    {
        FieldBlock &W = fields.field("Bcell_from_Bface_w", iblock);
        if (W.is_allocated() && metric_in_range(W, i, j, k))
        {
            for (int n = 0; n < 6; ++n)
            {
                out_xyz[0] += W(i, j, k, n) * phi[n];
                out_xyz[1] += W(i, j, k, 6 + n) * phi[n];
                out_xyz[2] += W(i, j, k, 12 + n) * phi[n];
            }
            return true;
        }
    }

    if (!fields.has_field("JDxi") || !fields.has_field("JDet") || !fields.has_field("JDze"))
        return false;

    FieldBlock &JDxi = fields.field("JDxi", iblock);
    FieldBlock &JDet = fields.field("JDet", iblock);
    FieldBlock &JDze = fields.field("JDze", iblock);
    if (!JDxi.is_allocated() || !JDet.is_allocated() || !JDze.is_allocated())
        return false;
    if (!metric_in_range(JDxi, i, j, k) || !metric_in_range(JDxi, i + 1, j, k) ||
        !metric_in_range(JDet, i, j, k) || !metric_in_range(JDet, i, j + 1, k) ||
        !metric_in_range(JDze, i, j, k) || !metric_in_range(JDze, i, j, k + 1))
        return false;

    double A[3][3] = {
        {0.5 * (JDxi(i, j, k, 0) + JDxi(i + 1, j, k, 0)),
         0.5 * (JDxi(i, j, k, 1) + JDxi(i + 1, j, k, 1)),
         0.5 * (JDxi(i, j, k, 2) + JDxi(i + 1, j, k, 2))},
        {0.5 * (JDet(i, j, k, 0) + JDet(i, j + 1, k, 0)),
         0.5 * (JDet(i, j, k, 1) + JDet(i, j + 1, k, 1)),
         0.5 * (JDet(i, j, k, 2) + JDet(i, j + 1, k, 2))},
        {0.5 * (JDze(i, j, k, 0) + JDze(i, j, k + 1, 0)),
         0.5 * (JDze(i, j, k, 1) + JDze(i, j, k + 1, 1)),
         0.5 * (JDze(i, j, k, 2) + JDze(i, j, k + 1, 2))}};

    double rhs[3] = {
        0.5 * (Fxi(i, j, k, 0) + Fxi(i + 1, j, k, 0)),
        0.5 * (Feta(i, j, k, 0) + Feta(i, j + 1, k, 0)),
        0.5 * (Fzeta(i, j, k, 0) + Fzeta(i, j, k + 1, 0))};

    return metric_solve3x3(A, rhs, out_xyz);
}

bool reconstruct_edge_1form_to_cell(Field &fields,
                                    int fid_xi, int fid_eta, int fid_zeta,
                                    int iblock, int i, int j, int k,
                                    double out_xyz[3])
{
    metric_zero(out_xyz);

    FieldBlock &Exi = fields.field(fid_xi, iblock);
    FieldBlock &Eeta = fields.field(fid_eta, iblock);
    FieldBlock &Ezeta = fields.field(fid_zeta, iblock);

    if (!Exi.is_allocated() || !Eeta.is_allocated() || !Ezeta.is_allocated())
        return false;

    if (!metric_in_range(Exi, i, j, k) ||
        !metric_in_range(Exi, i, j + 1, k) ||
        !metric_in_range(Exi, i, j, k + 1) ||
        !metric_in_range(Exi, i, j + 1, k + 1) ||
        !metric_in_range(Eeta, i, j, k) ||
        !metric_in_range(Eeta, i + 1, j, k) ||
        !metric_in_range(Eeta, i, j, k + 1) ||
        !metric_in_range(Eeta, i + 1, j, k + 1) ||
        !metric_in_range(Ezeta, i, j, k) ||
        !metric_in_range(Ezeta, i + 1, j, k) ||
        !metric_in_range(Ezeta, i, j + 1, k) ||
        !metric_in_range(Ezeta, i + 1, j + 1, k))
        return false;

    const double s[12] = {
        Exi(i, j, k, 0),
        Exi(i, j + 1, k, 0),
        Exi(i, j, k + 1, 0),
        Exi(i, j + 1, k + 1, 0),
        Eeta(i, j, k, 0),
        Eeta(i + 1, j, k, 0),
        Eeta(i, j, k + 1, 0),
        Eeta(i + 1, j, k + 1, 0),
        Ezeta(i, j, k, 0),
        Ezeta(i + 1, j, k, 0),
        Ezeta(i, j + 1, k, 0),
        Ezeta(i + 1, j + 1, k, 0)};

    if (fields.has_field("Jcell_from_Jedge_w"))
    {
        FieldBlock &W = fields.field("Jcell_from_Jedge_w", iblock);
        if (W.is_allocated() && metric_in_range(W, i, j, k))
        {
            for (int n = 0; n < 12; ++n)
            {
                out_xyz[0] += W(i, j, k, n) * s[n];
                out_xyz[1] += W(i, j, k, 12 + n) * s[n];
                out_xyz[2] += W(i, j, k, 24 + n) * s[n];
            }
            return true;
        }
    }

    if (!fields.has_field("pinvGT_cell"))
        return false;
    FieldBlock &pinv = fields.field("pinvGT_cell", iblock);
    if (!pinv.is_allocated() || !metric_in_range(pinv, i, j, k))
        return false;

    const double w[3] = {
        0.25 * (s[0] + s[1] + s[2] + s[3]),
        0.25 * (s[4] + s[5] + s[6] + s[7]),
        0.25 * (s[8] + s[9] + s[10] + s[11])};

    for (int r = 0; r < 3; ++r)
        out_xyz[r] = pinv(i, j, k, 3 * r + 0) * w[0] +
                     pinv(i, j, k, 3 * r + 1) * w[1] +
                     pinv(i, j, k, 3 * r + 2) * w[2];

    return true;
}

bool reconstruct_face_2form_to_node(Field &fields,
                                    int fid_xi, int fid_eta, int fid_zeta,
                                    int iblock, int i, int j, int k,
                                    double out_xyz[3])
{
    metric_zero(out_xyz);

    const FieldBlock &ref = fields.field(fid_xi, iblock);
    const Block &blk = ref.get_block();
    const int dim = blk.dimension;

    double sum[3] = {0.0, 0.0, 0.0};
    int count = 0;
    for (int di = -1; di <= 0; ++di)
        for (int dj = -1; dj <= 0; ++dj)
            for (int dk = (dim == 3 ? -1 : 0); dk <= 0; ++dk)
            {
                if (dim == 2 && dk != 0)
                    continue;
                double cell_xyz[3];
                if (!reconstruct_face_2form_to_cell(fields, fid_xi, fid_eta, fid_zeta,
                                                     iblock, i + di, j + dj, k + dk, cell_xyz))
                    continue;
                sum[0] += cell_xyz[0];
                sum[1] += cell_xyz[1];
                sum[2] += cell_xyz[2];
                ++count;
            }

    if (count == 0)
        return false;

    out_xyz[0] = sum[0] / count;
    out_xyz[1] = sum[1] / count;
    out_xyz[2] = sum[2] / count;
    return true;
}

bool reconstruct_edge_1form_to_node(Field &fields,
                                    int fid_xi, int fid_eta, int fid_zeta,
                                    int iblock, int i, int j, int k,
                                    double out_xyz[3])
{
    metric_zero(out_xyz);

    const FieldBlock &ref = fields.field(fid_xi, iblock);
    const Block &blk = ref.get_block();
    const int dim = blk.dimension;

    double sum[3] = {0.0, 0.0, 0.0};
    int count = 0;
    for (int di = -1; di <= 0; ++di)
        for (int dj = -1; dj <= 0; ++dj)
            for (int dk = (dim == 3 ? -1 : 0); dk <= 0; ++dk)
            {
                if (dim == 2 && dk != 0)
                    continue;
                double cell_xyz[3];
                if (!reconstruct_edge_1form_to_cell(fields, fid_xi, fid_eta, fid_zeta,
                                                     iblock, i + di, j + dj, k + dk, cell_xyz))
                    continue;
                sum[0] += cell_xyz[0];
                sum[1] += cell_xyz[1];
                sum[2] += cell_xyz[2];
                ++count;
            }

    if (count == 0)
        return false;

    out_xyz[0] = sum[0] / count;
    out_xyz[1] = sum[1] / count;
    out_xyz[2] = sum[2] / count;
    return true;
}
}
