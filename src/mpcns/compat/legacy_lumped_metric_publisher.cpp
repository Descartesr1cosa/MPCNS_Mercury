#include "mpcns/compat/legacy_lumped_metric_publisher.hpp"

#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "mpcns/discrete/dual_face_inverse_hodge.hpp"
#include "mpcns/geometry/geometry_primitives.hpp"

#include <array>

namespace mpcns::compat
{
namespace
{

std::array<double, 3> minus(const std::array<double, 3> &lhs,
                               const std::array<double, 3> &rhs)
{
    return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
}

} // namespace

void populate_legacy_lumped_face_and_edge_metrics(Field &fields, Grid &grid)
{
    Grid *grd = &grid;
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

} // namespace mpcns::compat
