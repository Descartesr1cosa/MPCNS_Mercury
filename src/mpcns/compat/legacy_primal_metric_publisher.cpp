#include "mpcns/compat/legacy_primal_metric_publisher.hpp"

#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "mpcns/geometry/geometry_primitives.hpp"

#include <array>

namespace mpcns::compat
{
void populate_legacy_primal_face_and_cell_metrics(Field &fields, Grid &grid)
{
    Grid *grd = &grid;
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
}

void populate_legacy_primal_edge_vectors(Field &fields, Grid &grid)
{
    Grid *grd = &grid;
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
}
} // namespace mpcns::compat
