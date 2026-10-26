#include "Z0_Tests.h"

#include "Z0_Boundary.h"
#include "Z0_TestCommon.h"
#include "7_metric/Metric.h"
#include "8_dec/DecOps.h"
#include "mpcns/compat/legacy_geometry_field_ids.hpp"
#include "mpcns/compat/legacy_geometry_field_registration.hpp"
#include "mpcns/compat/legacy_geometry_snapshot.hpp"
#include "mpcns/discrete/incidence_operator.hpp"
#include "mpcns/discrete/lumped_hodge.hpp"
#include "mpcns/discrete/dual_face_inverse_hodge.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

namespace
{
    void fill_phi(Field &field)
    {
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &phi = field.field("phi", ib);
            if (!phi.is_allocated())
                continue;
            const Int3 lo = phi.get_lo();
            const Int3 hi = phi.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        phi(i, j, k, 0) =
                            0.125 * i * i - 0.25 * j * k + 0.5 * k + 3.0;
        }
    }

    void fill_edge_form(Field &field)
    {
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            auto fill = [&](const char *name, double bias)
            {
                FieldBlock &fb = field.field(name, ib);
                if (!fb.is_allocated())
                    return;
                const Int3 lo = fb.get_lo();
                const Int3 hi = fb.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            fb(i, j, k, 0) =
                                bias + 0.375 * i - 0.5 * j + 0.625 * k +
                                0.0625 * i * j - 0.03125 * j * k;
            };
            fill("E_xi", 1.0);
            fill("E_eta", -2.0);
            fill("E_zeta", 0.75);
        }
    }

    void fill_face_form(Field &field)
    {
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            auto fill = [&](const char *name, double bias)
            {
                FieldBlock &fb = field.field(name, ib);
                if (!fb.is_allocated())
                    return;
                const Int3 lo = fb.get_lo();
                const Int3 hi = fb.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            fb(i, j, k, 0) =
                                bias + 0.25 * i - 0.375 * j + 0.5 * k +
                                0.03125 * i * k - 0.0625 * j * k;
            };
            fill("B_xi", 1.25);
            fill("B_eta", -0.75);
            fill("B_zeta", 2.0);
        }
    }

    double max_abs_owned(Field &field, const std::string &name)
    {
        double max_abs = 0.0;
        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &fb = field.field(fid, ib);
            if (!fb.is_allocated())
                continue;
            const Box3 owned = Z0_TEST::owned_box(field, ib, desc.location);
            for (int i = owned.lo.i; i < owned.hi.i; ++i)
                for (int j = owned.lo.j; j < owned.hi.j; ++j)
                    for (int k = owned.lo.k; k < owned.hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                            max_abs = std::max(max_abs, std::abs(fb(i, j, k, m)));
        }
        return Z0_TEST::global_max(max_abs);
    }

    // Regression oracle for the production CT curl stencil as it existed
    // before the d1 implementation was centralized.  Keep this formula in
    // the test only: numerical code must use mpcns::discrete::apply_d1.
    double max_difference_from_legacy_ct_d1(Field &field,
                                            const DEC::EdgeFormNames &edge_names,
                                            const DEC::FaceFormNames &face_names)
    {
        double max_difference = 0.0;
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &edge_xi = field.field(edge_names.xi, ib);
            FieldBlock &edge_eta = field.field(edge_names.eta, ib);
            FieldBlock &edge_zeta = field.field(edge_names.zeta, ib);
            FieldBlock &face_xi = field.field(face_names.xi, ib);
            FieldBlock &face_eta = field.field(face_names.eta, ib);
            FieldBlock &face_zeta = field.field(face_names.zeta, ib);

            if (!edge_xi.is_allocated() || !edge_eta.is_allocated() ||
                !edge_zeta.is_allocated() || !face_xi.is_allocated() ||
                !face_eta.is_allocated() || !face_zeta.is_allocated())
                continue;

            const auto update_max = [&](const FieldBlock &face, int i, int j, int k, double expected)
            {
                max_difference = std::max(
                    max_difference,
                    std::abs(face(i, j, k, 0) - expected));
            };

            {
                const Int3 lo = face_xi.inner_lo();
                const Int3 hi = face_xi.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            update_max(
                                face_xi, i, j, k,
                                (edge_eta(i, j, k, 0) - edge_eta(i, j, k + 1, 0)) +
                                    (edge_zeta(i, j + 1, k, 0) - edge_zeta(i, j, k, 0)));
            }

            {
                const Int3 lo = face_eta.inner_lo();
                const Int3 hi = face_eta.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            update_max(
                                face_eta, i, j, k,
                                (edge_xi(i, j, k + 1, 0) - edge_xi(i, j, k, 0)) +
                                    (edge_zeta(i, j, k, 0) - edge_zeta(i + 1, j, k, 0)));
            }

            {
                const Int3 lo = face_zeta.inner_lo();
                const Int3 hi = face_zeta.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            update_max(
                                face_zeta, i, j, k,
                                (edge_xi(i, j, k, 0) - edge_xi(i, j + 1, k, 0)) +
                                    (edge_eta(i + 1, j, k, 0) - edge_eta(i, j, k, 0)));
            }
        }
        return Z0_TEST::global_max(max_difference);
    }

    void apply_legacy_current_from_magnetic_flux(Field &field,
                                                  const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &b_xi = field.field("B_xi", ib);
            FieldBlock &b_eta = field.field("B_eta", ib);
            FieldBlock &b_zeta = field.field("B_zeta", ib);
            FieldBlock &j_xi = field.field("J_legacy_xi", ib);
            FieldBlock &j_eta = field.field("J_legacy_eta", ib);
            FieldBlock &j_zeta = field.field("J_legacy_zeta", ib);
            FieldBlock &beta_xi = field.field(geometry_ids.lumped_star_primal_2_to_dual_1.xi, ib);
            FieldBlock &beta_eta = field.field(geometry_ids.lumped_star_primal_2_to_dual_1.eta, ib);
            FieldBlock &beta_zeta = field.field(geometry_ids.lumped_star_primal_2_to_dual_1.zeta, ib);
            FieldBlock &alpha_xi = field.field(geometry_ids.lumped_inverse_star_primal_1.xi, ib);
            FieldBlock &alpha_eta = field.field(geometry_ids.lumped_inverse_star_primal_1.eta, ib);
            FieldBlock &alpha_zeta = field.field(geometry_ids.lumped_inverse_star_primal_1.zeta, ib);

            if (!b_xi.is_allocated() || !b_eta.is_allocated() || !b_zeta.is_allocated() ||
                !j_xi.is_allocated() || !j_eta.is_allocated() || !j_zeta.is_allocated())
                continue;

            {
                const Int3 lo = j_xi.inner_lo();
                const Int3 hi = j_xi.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double circulation =
                                (beta_eta(i, j, k - 1, 0) * b_eta(i, j, k - 1, 0) - beta_eta(i, j, k, 0) * b_eta(i, j, k, 0)) +
                                (beta_zeta(i, j, k, 0) * b_zeta(i, j, k, 0) - beta_zeta(i, j - 1, k, 0) * b_zeta(i, j - 1, k, 0));
                            j_xi(i, j, k, 0) = circulation;
                            j_xi(i, j, k, 0) *= alpha_xi(i, j, k, 0);
                        }
            }

            {
                const Int3 lo = j_eta.inner_lo();
                const Int3 hi = j_eta.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double circulation =
                                (beta_xi(i, j, k, 0) * b_xi(i, j, k, 0) - beta_xi(i, j, k - 1, 0) * b_xi(i, j, k - 1, 0)) +
                                (beta_zeta(i - 1, j, k, 0) * b_zeta(i - 1, j, k, 0) - beta_zeta(i, j, k, 0) * b_zeta(i, j, k, 0));
                            j_eta(i, j, k, 0) = circulation;
                            j_eta(i, j, k, 0) *= alpha_eta(i, j, k, 0);
                        }
            }

            {
                const Int3 lo = j_zeta.inner_lo();
                const Int3 hi = j_zeta.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double circulation =
                                (beta_xi(i, j - 1, k, 0) * b_xi(i, j - 1, k, 0) - beta_xi(i, j, k, 0) * b_xi(i, j, k, 0)) +
                                (beta_eta(i, j, k, 0) * b_eta(i, j, k, 0) - beta_eta(i - 1, j, k, 0) * b_eta(i - 1, j, k, 0));
                            j_zeta(i, j, k, 0) = circulation;
                            j_zeta(i, j, k, 0) *= alpha_zeta(i, j, k, 0);
                        }
            }
        }
    }

    double max_difference_inner(Field &field,
                                const std::string &lhs_name,
                                const std::string &rhs_name)
    {
        double max_difference = 0.0;
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &lhs = field.field(lhs_name, ib);
            FieldBlock &rhs = field.field(rhs_name, ib);
            if (!lhs.is_allocated() || !rhs.is_allocated())
                continue;

            const Int3 lo = lhs.inner_lo();
            const Int3 hi = lhs.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double lhs_value = lhs(i, j, k, 0);
                        const double rhs_value = rhs(i, j, k, 0);
                        if (!std::isfinite(lhs_value) || !std::isfinite(rhs_value))
                            max_difference = std::numeric_limits<double>::infinity();
                        else
                            max_difference = std::max(
                                max_difference,
                                std::abs(lhs_value - rhs_value));
                    }
        }
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_legacy_primal_face_area_stencil(
        Grid &grid,
        Field &field,
        const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
        using Point = std::array<double, 3>;
        using Offset = std::array<int, 3>;

        const auto legacy_area = [](const Point &p00,
                                    const Point &p10,
                                    const Point &p01,
                                    const Point &p11)
        {
            const auto subtract = [](const Point &lhs, const Point &rhs)
            {
                return Point{lhs[0] - rhs[0],
                             lhs[1] - rhs[1],
                             lhs[2] - rhs[2]};
            };
            const auto cross = [](const Point &lhs, const Point &rhs)
            {
                return Point{lhs[1] * rhs[2] - lhs[2] * rhs[1],
                             lhs[2] * rhs[0] - lhs[0] * rhs[2],
                             lhs[0] * rhs[1] - lhs[1] * rhs[0]};
            };

            const Point first_triangle =
                cross(subtract(p10, p00), subtract(p01, p00));
            const Point second_triangle =
                cross(subtract(p01, p11), subtract(p10, p11));
            return Point{0.5 * (first_triangle[0] + second_triangle[0]),
                         0.5 * (first_triangle[1] + second_triangle[1]),
                         0.5 * (first_triangle[2] + second_triangle[2])};
        };

        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };

        const auto compare_direction = [&](const Block &block,
                                           FieldBlock &actual,
                                           const Offset &p10_offset,
                                           const Offset &p01_offset,
                                           const Offset &p11_offset)
        {
            if (!actual.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                return;
            }

            const auto point = [&](int i, int j, int k)
            {
                return Point{block.x(i, j, k),
                             block.y(i, j, k),
                             block.z(i, j, k)};
            };
            const Int3 lo = actual.get_lo();
            const Int3 hi = actual.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const Point expected = legacy_area(
                            point(i, j, k),
                            point(i + p10_offset[0],
                                  j + p10_offset[1],
                                  k + p10_offset[2]),
                            point(i + p01_offset[0],
                                  j + p01_offset[1],
                                  k + p01_offset[2]),
                            point(i + p11_offset[0],
                                  j + p11_offset[1],
                                  k + p11_offset[2]));
                        for (int component = 0; component < 3; ++component)
                            update(expected[component],
                                   actual(i, j, k, component));
                    }
        };

        const Offset xi_p10{{0, 1, 0}};
        const Offset xi_p01{{0, 0, 1}};
        const Offset xi_p11{{0, 1, 1}};
        const Offset eta_p10{{0, 0, 1}};
        const Offset eta_p01{{1, 0, 0}};
        const Offset eta_p11{{1, 0, 1}};
        const Offset zeta_p10{{1, 0, 0}};
        const Offset zeta_p01{{0, 1, 0}};
        const Offset zeta_p11{{1, 1, 0}};

        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const Block &block = grid.grids(ib);
            compare_direction(block,
                              field.field(geometry_ids.primal_face_area_vector.xi, ib),
                              xi_p10, xi_p01, xi_p11);
            compare_direction(block,
                              field.field(geometry_ids.primal_face_area_vector.eta, ib),
                              eta_p10, eta_p01, eta_p11);
            compare_direction(block,
                              field.field(geometry_ids.primal_face_area_vector.zeta, ib),
                              zeta_p10, zeta_p01, zeta_p11);
        }
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_legacy_cell_volume_stencil(
        Grid &grid,
        Field &field,
        const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
        using Point = std::array<double, 3>;

        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };
        const auto subtract = [](const Point &lhs, const Point &rhs)
        {
            return Point{lhs[0] - rhs[0],
                         lhs[1] - rhs[1],
                         lhs[2] - rhs[2]};
        };
        const auto dot = [](const Point &lhs, const Point &rhs)
        {
            return lhs[0] * rhs[0] +
                   lhs[1] * rhs[1] +
                   lhs[2] * rhs[2];
        };
        const auto face_center = [](const Point &p00,
                                    const Point &p10,
                                    const Point &p01,
                                    const Point &p11)
        {
            return Point{0.25 * (p00[0] + p10[0] + p01[0] + p11[0]),
                         0.25 * (p00[1] + p10[1] + p01[1] + p11[1]),
                         0.25 * (p00[2] + p10[2] + p01[2] + p11[2])};
        };

        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            FieldBlock &jac = field.field(geometry_ids.primal_cell_volume, ib);
            FieldBlock &area_xi =
                field.field(geometry_ids.primal_face_area_vector.xi, ib);
            FieldBlock &area_eta =
                field.field(geometry_ids.primal_face_area_vector.eta, ib);
            FieldBlock &area_zeta =
                field.field(geometry_ids.primal_face_area_vector.zeta, ib);
            if (!jac.is_allocated() || !area_xi.is_allocated() ||
                !area_eta.is_allocated() || !area_zeta.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                continue;
            }

            const Block &block = grid.grids(ib);
            const auto node = [&](int i, int j, int k)
            {
                return Point{block.x(i, j, k),
                             block.y(i, j, k),
                             block.z(i, j, k)};
            };
            const auto cell_center = [&](int i, int j, int k)
            {
                return Point{block.dual_x(i + 1, j + 1, k + 1),
                             block.dual_y(i + 1, j + 1, k + 1),
                             block.dual_z(i + 1, j + 1, k + 1)};
            };

            const Int3 lo = jac.get_lo();
            const Int3 hi = jac.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const Point center = cell_center(i, j, k);
                        double volume_sum = 0.0;

                        volume_sum += dot(
                            Point{-area_xi(i, j, k, 0),
                                  -area_xi(i, j, k, 1),
                                  -area_xi(i, j, k, 2)},
                            subtract(face_center(node(i, j, k),
                                                 node(i, j + 1, k),
                                                 node(i, j, k + 1),
                                                 node(i, j + 1, k + 1)),
                                     center));
                        volume_sum += dot(
                            Point{area_xi(i + 1, j, k, 0),
                                  area_xi(i + 1, j, k, 1),
                                  area_xi(i + 1, j, k, 2)},
                            subtract(face_center(node(i + 1, j, k),
                                                 node(i + 1, j + 1, k),
                                                 node(i + 1, j, k + 1),
                                                 node(i + 1, j + 1, k + 1)),
                                     center));
                        volume_sum += dot(
                            Point{-area_eta(i, j, k, 0),
                                  -area_eta(i, j, k, 1),
                                  -area_eta(i, j, k, 2)},
                            subtract(face_center(node(i, j, k),
                                                 node(i + 1, j, k),
                                                 node(i, j, k + 1),
                                                 node(i + 1, j, k + 1)),
                                     center));
                        volume_sum += dot(
                            Point{area_eta(i, j + 1, k, 0),
                                  area_eta(i, j + 1, k, 1),
                                  area_eta(i, j + 1, k, 2)},
                            subtract(face_center(node(i, j + 1, k),
                                                 node(i + 1, j + 1, k),
                                                 node(i, j + 1, k + 1),
                                                 node(i + 1, j + 1, k + 1)),
                                     center));
                        volume_sum += dot(
                            Point{-area_zeta(i, j, k, 0),
                                  -area_zeta(i, j, k, 1),
                                  -area_zeta(i, j, k, 2)},
                            subtract(face_center(node(i, j, k),
                                                 node(i + 1, j, k),
                                                 node(i, j + 1, k),
                                                 node(i + 1, j + 1, k)),
                                     center));
                        volume_sum += dot(
                            Point{area_zeta(i, j, k + 1, 0),
                                  area_zeta(i, j, k + 1, 1),
                                  area_zeta(i, j, k + 1, 2)},
                            subtract(face_center(node(i, j, k + 1),
                                                 node(i + 1, j, k + 1),
                                                 node(i, j + 1, k + 1),
                                                 node(i + 1, j + 1, k + 1)),
                                     center));
                        update(volume_sum / 3.0, jac(i, j, k, 0));
                    }
        }
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_legacy_face_hodge_metric_stencil(
        Grid &grid,
        Field &field,
        const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
        constexpr double eps = 1e-25;
        using Point = std::array<double, 3>;
        using Offset = std::array<int, 3>;

        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };
        const auto subtract = [](const Point &lhs, const Point &rhs)
        {
            return Point{lhs[0] - rhs[0],
                         lhs[1] - rhs[1],
                         lhs[2] - rhs[2]};
        };
        const auto norm = [](const Point &vector)
        {
            return std::sqrt(vector[0] * vector[0] +
                             vector[1] * vector[1] +
                             vector[2] * vector[2]);
        };

        const auto compare_direction = [&](const Block &block,
                                           FieldBlock &face_vector,
                                           FieldBlock &face_area,
                                           FieldBlock &dual_edge_length,
                                           FieldBlock &star_2,
                                           const Offset &left_offset,
                                           const Offset &right_offset)
        {
            if (!face_vector.is_allocated() || !face_area.is_allocated() ||
                !dual_edge_length.is_allocated() || !star_2.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                return;
            }

            const auto cell_center = [&](int i, int j, int k)
            {
                return Point{block.dual_x(i, j, k),
                             block.dual_y(i, j, k),
                             block.dual_z(i, j, k)};
            };
            const Int3 lo = face_area.get_lo();
            const Int3 hi = face_area.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const Point area_vector{
                            face_vector(i, j, k, 0),
                            face_vector(i, j, k, 1),
                            face_vector(i, j, k, 2)};
                        const double area = norm(area_vector);
                        const Point left = cell_center(
                            i + left_offset[0],
                            j + left_offset[1],
                            k + left_offset[2]);
                        const Point right = cell_center(
                            i + right_offset[0],
                            j + right_offset[1],
                            k + right_offset[2]);
                        const double dual_length = norm(subtract(right, left));
                        const double hodge_star_2 =
                            (area < eps) ? 0.0 : dual_length / area;

                        update(area, face_area(i, j, k, 0));
                        update(dual_length, dual_edge_length(i, j, k, 0));
                        update(hodge_star_2, star_2(i, j, k, 0));
                    }
        };

        const Offset xi_left{{0, 1, 1}};
        const Offset xi_right{{1, 1, 1}};
        const Offset eta_left{{1, 0, 1}};
        const Offset eta_right{{1, 1, 1}};
        const Offset zeta_left{{1, 1, 0}};
        const Offset zeta_right{{1, 1, 1}};

        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const Block &block = grid.grids(ib);
            compare_direction(
                block,
                field.field(geometry_ids.primal_face_area_vector.xi, ib),
                field.field(geometry_ids.primal_face_area.xi, ib),
                field.field(geometry_ids.dual_edge_length_at_face.xi, ib),
                field.field(geometry_ids.lumped_star_primal_2_to_dual_1.xi, ib),
                xi_left, xi_right);
            compare_direction(
                block,
                field.field(geometry_ids.primal_face_area_vector.eta, ib),
                field.field(geometry_ids.primal_face_area.eta, ib),
                field.field(geometry_ids.dual_edge_length_at_face.eta, ib),
                field.field(geometry_ids.lumped_star_primal_2_to_dual_1.eta, ib),
                eta_left, eta_right);
            compare_direction(
                block,
                field.field(geometry_ids.primal_face_area_vector.zeta, ib),
                field.field(geometry_ids.primal_face_area.zeta, ib),
                field.field(geometry_ids.dual_edge_length_at_face.zeta, ib),
                field.field(geometry_ids.lumped_star_primal_2_to_dual_1.zeta, ib),
                zeta_left, zeta_right);
        }
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_legacy_dual_face_area_stencil(
        Grid &grid,
        Field &field,
        const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
        using Point = std::array<double, 3>;
        using Offset = std::array<int, 3>;

        const auto legacy_area = [](const Point &p00,
                                    const Point &p10,
                                    const Point &p01,
                                    const Point &p11)
        {
            const auto subtract = [](const Point &lhs, const Point &rhs)
            {
                return Point{lhs[0] - rhs[0],
                             lhs[1] - rhs[1],
                             lhs[2] - rhs[2]};
            };
            const auto cross = [](const Point &lhs, const Point &rhs)
            {
                return Point{lhs[1] * rhs[2] - lhs[2] * rhs[1],
                             lhs[2] * rhs[0] - lhs[0] * rhs[2],
                             lhs[0] * rhs[1] - lhs[1] * rhs[0]};
            };

            const Point first_triangle =
                cross(subtract(p10, p00), subtract(p01, p00));
            const Point second_triangle =
                cross(subtract(p01, p11), subtract(p10, p11));
            return Point{0.5 * (first_triangle[0] + second_triangle[0]),
                         0.5 * (first_triangle[1] + second_triangle[1]),
                         0.5 * (first_triangle[2] + second_triangle[2])};
        };

        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };
        const auto compare_direction = [&](const Block &block,
                                           FieldBlock &actual,
                                           const Offset &p00_offset,
                                           const Offset &p10_offset,
                                           const Offset &p01_offset,
                                           const Offset &p11_offset)
        {
            if (!actual.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                return;
            }

            const auto dual_cell_center = [&](int i, int j, int k)
            {
                return Point{block.dual_x(i, j, k),
                             block.dual_y(i, j, k),
                             block.dual_z(i, j, k)};
            };
            const Int3 lo = actual.get_lo();
            const Int3 hi = actual.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const Point expected = legacy_area(
                            dual_cell_center(i + p00_offset[0],
                                             j + p00_offset[1],
                                             k + p00_offset[2]),
                            dual_cell_center(i + p10_offset[0],
                                             j + p10_offset[1],
                                             k + p10_offset[2]),
                            dual_cell_center(i + p01_offset[0],
                                             j + p01_offset[1],
                                             k + p01_offset[2]),
                            dual_cell_center(i + p11_offset[0],
                                             j + p11_offset[1],
                                             k + p11_offset[2]));
                        for (int component = 0; component < 3; ++component)
                            update(expected[component],
                                   actual(i, j, k, component));
                    }
        };

        const Offset xi_p00{{1, 0, 0}};
        const Offset xi_p10{{1, 1, 0}};
        const Offset xi_p01{{1, 0, 1}};
        const Offset xi_p11{{1, 1, 1}};
        const Offset eta_p00{{0, 1, 0}};
        const Offset eta_p10{{1, 1, 0}};
        const Offset eta_p01{{0, 1, 1}};
        const Offset eta_p11{{1, 1, 1}};
        const Offset zeta_p00{{0, 0, 1}};
        const Offset zeta_p10{{1, 0, 1}};
        const Offset zeta_p01{{0, 1, 1}};
        const Offset zeta_p11{{1, 1, 1}};

        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const Block &block = grid.grids(ib);
            compare_direction(
                block,
                field.field(geometry_ids.dual_face_area_vector_at_edge.xi, ib),
                xi_p00, xi_p10, xi_p01, xi_p11);
            compare_direction(
                block,
                field.field(geometry_ids.dual_face_area_vector_at_edge.eta, ib),
                eta_p00, eta_p10, eta_p01, eta_p11);
            compare_direction(
                block,
                field.field(geometry_ids.dual_face_area_vector_at_edge.zeta, ib),
                zeta_p00, zeta_p10, zeta_p01, zeta_p11);
        }
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_legacy_edge_dual_measure_stencil(
        Grid &grid,
        Field &field,
        const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
        using Point = std::array<double, 3>;
        using Offset = std::array<int, 3>;

        const auto subtract = [](const Point &lhs, const Point &rhs)
        {
            return Point{lhs[0] - rhs[0],
                         lhs[1] - rhs[1],
                         lhs[2] - rhs[2]};
        };
        const auto legacy_norm = [](const Point &vector)
        {
            return std::sqrt(vector[0] * vector[0] +
                             vector[1] * vector[1] +
                             vector[2] * vector[2]);
        };

        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };
        const auto compare_direction = [&](const Block &block,
                                           FieldBlock &edge_length,
                                           FieldBlock &dual_face_vector,
                                           FieldBlock &dual_face_area,
                                           const Offset &endpoint_offset)
        {
            if (!edge_length.is_allocated() || !dual_face_vector.is_allocated() ||
                !dual_face_area.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                return;
            }

            const auto node = [&](int i, int j, int k)
            {
                return Point{block.x(i, j, k),
                             block.y(i, j, k),
                             block.z(i, j, k)};
            };
            const Int3 lo = edge_length.get_lo();
            const Int3 hi = edge_length.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const Point edge = subtract(
                            node(i + endpoint_offset[0],
                                 j + endpoint_offset[1],
                                 k + endpoint_offset[2]),
                            node(i, j, k));
                        const Point face_vector{
                            dual_face_vector(i, j, k, 0),
                            dual_face_vector(i, j, k, 1),
                            dual_face_vector(i, j, k, 2)};
                        update(legacy_norm(edge), edge_length(i, j, k, 0));
                        update(legacy_norm(face_vector),
                               dual_face_area(i, j, k, 0));
                    }
        };

        const Offset xi_endpoint{{1, 0, 0}};
        const Offset eta_endpoint{{0, 1, 0}};
        const Offset zeta_endpoint{{0, 0, 1}};
        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const Block &block = grid.grids(ib);
            compare_direction(
                block,
                field.field(geometry_ids.primal_edge_length.xi, ib),
                field.field(geometry_ids.dual_face_area_vector_at_edge.xi, ib),
                field.field(geometry_ids.dual_face_area_at_edge.xi, ib),
                xi_endpoint);
            compare_direction(
                block,
                field.field(geometry_ids.primal_edge_length.eta, ib),
                field.field(geometry_ids.dual_face_area_vector_at_edge.eta, ib),
                field.field(geometry_ids.dual_face_area_at_edge.eta, ib),
                eta_endpoint);
            compare_direction(
                block,
                field.field(geometry_ids.primal_edge_length.zeta, ib),
                field.field(geometry_ids.dual_face_area_vector_at_edge.zeta, ib),
                field.field(geometry_ids.dual_face_area_at_edge.zeta, ib),
                zeta_endpoint);
        }
        return Z0_TEST::global_max(max_difference);
    }

    double legacy_axis_safe_inverse_star_1(
        double primal_edge_length,
        const std::array<double, 3> &p00,
        const std::array<double, 3> &p10,
        const std::array<double, 3> &p01,
        const std::array<double, 3> &p11,
        double dual_face_area,
        double zero_area_tolerance)
    {
        constexpr double dual_collapse_ratio = 1.0e-3;
        constexpr double primal_edge_absolute_tolerance = 1.0e-12;
        const auto subtract = [](const std::array<double, 3> &lhs,
                                 const std::array<double, 3> &rhs)
        {
            return std::array<double, 3>{lhs[0] - rhs[0],
                                         lhs[1] - rhs[1],
                                         lhs[2] - rhs[2]};
        };
        const auto norm = [](const std::array<double, 3> &vector)
        {
            return std::sqrt(vector[0] * vector[0] +
                             vector[1] * vector[1] +
                             vector[2] * vector[2]);
        };

        const double inverse_star_raw =
            (dual_face_area < zero_area_tolerance)
                ? 0.0
                : (primal_edge_length / dual_face_area);
        const double d00_10 = norm(subtract(p10, p00));
        const double d00_01 = norm(subtract(p01, p00));
        const double d10_11 = norm(subtract(p11, p10));
        const double d01_11 = norm(subtract(p11, p01));
        const double d00_11 = norm(subtract(p11, p00));
        const double d10_01 = norm(subtract(p01, p10));
        const double href = std::max(
            std::max(std::max(d00_10, d00_01), std::max(d10_11, d01_11)),
            std::max(std::max(d00_11, d10_01),
                     std::max(primal_edge_length, zero_area_tolerance)));
        const double side_min = std::min(
            std::min(d00_10, d00_01),
            std::min(d10_11, d01_11));
        const double diagonal_min = std::min(d00_11, d10_01);
        const bool dual_touches_axis =
            (primal_edge_length < primal_edge_absolute_tolerance) ||
            (side_min < dual_collapse_ratio * href) ||
            (diagonal_min < dual_collapse_ratio * href) ||
            (dual_face_area < dual_collapse_ratio * href * href);
        return dual_touches_axis ? 0.0 : inverse_star_raw;
    }

    double max_difference_from_axis_safe_inverse_star_policy_cases()
    {
        using Point = std::array<double, 3>;
        constexpr double zero_area_tolerance = 1.0e-25;
        const mpcns::discrete::AxisSafeInverseStar1Policy policy;
        const Point p00{0.0, 0.0, 0.0};
        const Point p10{1.0, 0.0, 0.0};
        const Point p01{0.0, 1.0, 0.0};
        const Point p11{1.0, 1.0, 0.0};
        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };
        const auto compare = [&](double primal_edge_length,
                                 const Point &case_p00,
                                 const Point &case_p10,
                                 const Point &case_p01,
                                 const Point &case_p11,
                                 double dual_face_area)
        {
            const double expected = legacy_axis_safe_inverse_star_1(
                primal_edge_length, case_p00, case_p10, case_p01, case_p11,
                dual_face_area, zero_area_tolerance);
            const double actual = policy.apply(
                primal_edge_length,
                {case_p00, case_p10, case_p01, case_p11},
                dual_face_area,
                zero_area_tolerance);
            update(expected, actual);
        };

        compare(1.0, p00, p10, p01, p11, 1.0);
        compare(0.0, p00, p10, p01, p11, 1.0);
        compare(1.0, p00, p00, p01, p11, 1.0);
        compare(1.0, p00, p10, Point{0.0, 1.0e-4, 0.0},
                Point{1.0, 1.0e-4, 0.0}, 1.0e-4);
        compare(1.0, p00, p10, p01, p11, 1.0e-30);
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_legacy_axis_safe_inverse_star_stencil(
        Grid &grid,
        Field &field,
        const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
        using Point = std::array<double, 3>;
        using Offset = std::array<int, 3>;
        constexpr double zero_area_tolerance = 1.0e-25;
        const auto subtract = [](const Point &lhs, const Point &rhs)
        {
            return Point{lhs[0] - rhs[0],
                         lhs[1] - rhs[1],
                         lhs[2] - rhs[2]};
        };
        const auto norm = [](const Point &vector)
        {
            return std::sqrt(vector[0] * vector[0] +
                             vector[1] * vector[1] +
                             vector[2] * vector[2]);
        };

        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };
        const auto compare_direction = [&](const Block &block,
                                           FieldBlock &inverse_star,
                                           FieldBlock &dual_face_vector,
                                           const Offset &endpoint_offset,
                                           const Offset &p00_offset,
                                           const Offset &p10_offset,
                                           const Offset &p01_offset,
                                           const Offset &p11_offset)
        {
            if (!inverse_star.is_allocated() || !dual_face_vector.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                return;
            }

            const auto node = [&](int i, int j, int k)
            {
                return Point{block.x(i, j, k),
                             block.y(i, j, k),
                             block.z(i, j, k)};
            };
            const auto dual_cell_center = [&](int i, int j, int k)
            {
                return Point{block.dual_x(i, j, k),
                             block.dual_y(i, j, k),
                             block.dual_z(i, j, k)};
            };
            const Int3 lo = inverse_star.get_lo();
            const Int3 hi = inverse_star.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double primal_edge_length = norm(subtract(
                            node(i + endpoint_offset[0],
                                 j + endpoint_offset[1],
                                 k + endpoint_offset[2]),
                            node(i, j, k)));
                        const Point p00 = dual_cell_center(
                            i + p00_offset[0], j + p00_offset[1],
                            k + p00_offset[2]);
                        const Point p10 = dual_cell_center(
                            i + p10_offset[0], j + p10_offset[1],
                            k + p10_offset[2]);
                        const Point p01 = dual_cell_center(
                            i + p01_offset[0], j + p01_offset[1],
                            k + p01_offset[2]);
                        const Point p11 = dual_cell_center(
                            i + p11_offset[0], j + p11_offset[1],
                            k + p11_offset[2]);
                        const Point dual_face{
                            dual_face_vector(i, j, k, 0),
                            dual_face_vector(i, j, k, 1),
                            dual_face_vector(i, j, k, 2)};
                        const double expected = legacy_axis_safe_inverse_star_1(
                            primal_edge_length, p00, p10, p01, p11,
                            norm(dual_face), zero_area_tolerance);
                        update(expected, inverse_star(i, j, k, 0));
                    }
        };

        const Offset xi_endpoint{{1, 0, 0}};
        const Offset eta_endpoint{{0, 1, 0}};
        const Offset zeta_endpoint{{0, 0, 1}};
        const Offset xi_p00{{1, 0, 0}};
        const Offset xi_p10{{1, 1, 0}};
        const Offset xi_p01{{1, 0, 1}};
        const Offset xi_p11{{1, 1, 1}};
        const Offset eta_p00{{0, 1, 0}};
        const Offset eta_p10{{1, 1, 0}};
        const Offset eta_p01{{0, 1, 1}};
        const Offset eta_p11{{1, 1, 1}};
        const Offset zeta_p00{{0, 0, 1}};
        const Offset zeta_p10{{1, 0, 1}};
        const Offset zeta_p01{{0, 1, 1}};
        const Offset zeta_p11{{1, 1, 1}};
        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const Block &block = grid.grids(ib);
            compare_direction(
                block,
                field.field(geometry_ids.lumped_inverse_star_primal_1.xi, ib),
                field.field(geometry_ids.dual_face_area_vector_at_edge.xi, ib),
                xi_endpoint, xi_p00, xi_p10, xi_p01, xi_p11);
            compare_direction(
                block,
                field.field(geometry_ids.lumped_inverse_star_primal_1.eta, ib),
                field.field(geometry_ids.dual_face_area_vector_at_edge.eta, ib),
                eta_endpoint, eta_p00, eta_p10, eta_p01, eta_p11);
            compare_direction(
                block,
                field.field(geometry_ids.lumped_inverse_star_primal_1.zeta, ib),
                field.field(geometry_ids.dual_face_area_vector_at_edge.zeta, ib),
                zeta_endpoint, zeta_p00, zeta_p10, zeta_p01, zeta_p11);
        }
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_geometry_snapshot(
        Grid &grid,
        Field &field,
        const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
        const mpcns::geometry::GeometrySnapshot snapshot =
            mpcns::compat::borrow_legacy_geometry_snapshot(grid, field);
        double max_difference = snapshot.is_complete() ? 0.0 : std::numeric_limits<double>::infinity();

        const auto update = [&](double lhs, double rhs)
        {
            if (!std::isfinite(lhs) || !std::isfinite(rhs))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference, std::abs(lhs - rhs));
        };

        const auto compare_field = [&](const mpcns::geometry::ConstFieldView &view,
                                       const FieldBlock &legacy,
                                       int components)
        {
            if (!legacy.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                return;
            }
            const Int3 index = legacy.get_lo();
            const mpcns::core::LogicalIndex geometry_index{index.i, index.j, index.k};
            for (int component = 0; component < components; ++component)
                update(view.value(geometry_index, component), legacy(index.i, index.j, index.k, component));
        };

        const auto compare_directional = [&](const mpcns::geometry::DirectionalFieldViews &views,
                                             const mpcns::compat::DirectionalFieldIds &ids,
                                             int components,
                                             int block)
        {
            compare_field(views.xi, field.field(ids.xi, block), components);
            compare_field(views.eta, field.field(ids.eta, block), components);
            compare_field(views.zeta, field.field(ids.zeta, block), components);
        };

        if (snapshot.block_count() != static_cast<std::size_t>(grid.nblock))
            max_difference = std::numeric_limits<double>::infinity();

        for (int ib = 0; ib < grid.nblock; ++ib)
        {
            const mpcns::geometry::GeometryBlock &geometry = snapshot.block(ib);
            const Block &legacy_block = grid.grids(ib);
            if (geometry.name != legacy_block.block_name ||
                geometry.spatial_dimension != legacy_block.dimension)
                max_difference = std::numeric_limits<double>::infinity();

            const mpcns::core::LogicalIndex node_index{0, 0, 0};
            const mpcns::geometry::CartesianPoint point = geometry.coordinates.position(node_index);
            update(point.x, legacy_block.x(0, 0, 0));
            update(point.y, legacy_block.y(0, 0, 0));
            update(point.z, legacy_block.z(0, 0, 0));

            compare_field(geometry.primal.cell_volume,
                          field.field(geometry_ids.primal_cell_volume, ib), 1);
            compare_directional(geometry.primal.face_area_vector,
                                geometry_ids.primal_face_area_vector, 3, ib);
            compare_directional(geometry.primal.edge_vector,
                                geometry_ids.primal_edge_vector, 3, ib);
            compare_directional(geometry.reconstruction.cell_covariant_basis,
                                geometry_ids.cell_covariant_basis, 3, ib);
            compare_directional(geometry.primal.face_area,
                                geometry_ids.primal_face_area, 1, ib);
            compare_directional(geometry.primal.edge_length,
                                geometry_ids.primal_edge_length, 1, ib);
            compare_directional(geometry.dual.edge_length_at_primal_face,
                                geometry_ids.dual_edge_length_at_face, 1, ib);
            compare_directional(geometry.dual.face_area_vector_at_primal_edge,
                                geometry_ids.dual_face_area_vector_at_edge, 3, ib);
            compare_directional(geometry.dual.face_area_at_primal_edge,
                                geometry_ids.dual_face_area_at_edge, 1, ib);
            compare_field(geometry.reconstruction.cell_gradient_pseudoinverse,
                          field.field(geometry_ids.cell_reconstruction_gradient_pseudoinverse, ib), 9);
            compare_directional(geometry.reconstruction.edge_gradient_pseudoinverse,
                                geometry_ids.edge_reconstruction_gradient_pseudoinverse, 9, ib);
            compare_directional(geometry.reconstruction.edge_area_pseudoinverse,
                                geometry_ids.edge_reconstruction_area_pseudoinverse, 9, ib);
        }
        return Z0_TEST::global_max(max_difference);
    }

    double legacy_geometry_catalog_layout_difference(
        const Field &field,
        int geometry_ghost)
    {
        double max_difference = 0.0;
        for (const auto &registration :
             mpcns::compat::legacy_geometry_field_registrations())
        {
            if (!field.has_field(registration.legacy_name))
            {
                max_difference = std::numeric_limits<double>::infinity();
                continue;
            }

            const FieldDescriptor &descriptor =
                field.descriptor(registration.legacy_name);
            const bool matches =
                descriptor.name == registration.legacy_name &&
                descriptor.location == registration.location &&
                descriptor.ncomp == registration.component_count &&
                descriptor.nghost == registration.ghost_layers(geometry_ghost);
            if (!matches)
                max_difference = std::numeric_limits<double>::infinity();
        }
        return Z0_TEST::global_max(max_difference);
    }
}

namespace Z0_TEST
{
    bool RunDecChainTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param)
    {
        (void)halo;
        (void)param;
        if (rank() == 0)
            std::cout << "Z0 DEC chain tests\n";

        bool passed = true;
        const DEC::EdgeFormNames e{"E_xi", "E_eta", "E_zeta"};
        const DEC::FaceFormNames b{"B_xi", "B_eta", "B_zeta"};

        fill_all_with_nan(field);
        const double geometry_catalog_difference =
            legacy_geometry_catalog_layout_difference(
                field, param.GetInt("ngg"));
        passed &= print_result(
            "Legacy geometry field layout matches catalog exactly",
            geometry_catalog_difference == 0.0,
            max_detail(geometry_catalog_difference));

        fill_phi(field);
        DEC::d0_node_to_edge(field, "phi", e);
        DEC::d1_edge_to_face(field, e, b);
        const double d1d0 = std::max(max_abs_owned(field, "B_xi"),
                                     std::max(max_abs_owned(field, "B_eta"),
                                              max_abs_owned(field, "B_zeta")));
        passed &= print_result("D1 * D0 == 0", d1d0 < 1.0e-12, max_detail(d1d0));

        fill_all_with_nan(field);
        fill_edge_form(field);
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            auto &edge_xi = field.field(e.xi, ib);
            auto &edge_eta = field.field(e.eta, ib);
            auto &edge_zeta = field.field(e.zeta, ib);
            auto &face_xi = field.field(b.xi, ib);
            auto &face_eta = field.field(b.eta, ib);
            auto &face_zeta = field.field(b.zeta, ib);

            if (!edge_xi.is_allocated() || !edge_eta.is_allocated() ||
                !edge_zeta.is_allocated() || !face_xi.is_allocated() ||
                !face_eta.is_allocated() || !face_zeta.is_allocated())
                continue;

            const mpcns::state::PrimalEdgeFormView edge{
                edge_xi, edge_eta, edge_zeta};
            const mpcns::state::PrimalFaceFormView face{
                face_xi, face_eta, face_zeta};
            mpcns::discrete::apply_d1(edge, face);
        }
        const double d1_legacy_difference =
            max_difference_from_legacy_ct_d1(field, e, b);
        passed &= print_result("MPCNS d1 matches legacy CT stencil exactly",
                               d1_legacy_difference == 0.0,
                               max_detail(d1_legacy_difference));
        DEC::d2_face_to_cell(field, b, "divB");
        const double d2d1 = max_abs_owned(field, "divB");
        passed &= print_result("MPCNS D2 * D1 == 0", d2d1 < 1.0e-12, max_detail(d2d1));

        // Rebuild the immutable metric fields after the NaN-based DEC tests,
        // then compare the new Hodge composition with the exact legacy fused
        // current stencil.  No halo exchange is required: this test fills the
        // full local ghost range deterministically.
        fill_all_with_nan(field);
        METRIC::compute_metric_fields(field, *field.grd);
        fill_face_form(field);
        const auto geometry_ids = mpcns::compat::bind_legacy_geometry_field_ids(field);
        const double primal_face_area_difference =
            max_difference_from_legacy_primal_face_area_stencil(
                *field.grd, field, geometry_ids);
        passed &= print_result("MPCNS quadrilateral area primitive matches legacy JD stencil exactly",
                               primal_face_area_difference == 0.0,
                               max_detail(primal_face_area_difference));
        const double cell_volume_difference =
            max_difference_from_legacy_cell_volume_stencil(
                *field.grd, field, geometry_ids);
        passed &= print_result("MPCNS cell-volume primitive matches legacy Jac stencil exactly",
                               cell_volume_difference == 0.0,
                               max_detail(cell_volume_difference));
        const double face_hodge_metric_difference =
            max_difference_from_legacy_face_hodge_metric_stencil(
                *field.grd, field, geometry_ids);
        passed &= print_result(
            "MPCNS face-area and star-2 metrics match legacy stencil exactly",
            face_hodge_metric_difference == 0.0,
            max_detail(face_hodge_metric_difference));
        const double dual_face_area_difference =
            max_difference_from_legacy_dual_face_area_stencil(
                *field.grd, field, geometry_ids);
        passed &= print_result(
            "MPCNS dual-face area primitive matches legacy Sstar stencil exactly",
            dual_face_area_difference == 0.0,
            max_detail(dual_face_area_difference));
        const double edge_dual_measure_difference =
            max_difference_from_legacy_edge_dual_measure_stencil(
                *field.grd, field, geometry_ids);
        passed &= print_result(
            "MPCNS edge length and dual-face magnitude match legacy stencil exactly",
            edge_dual_measure_difference == 0.0,
            max_detail(edge_dual_measure_difference));
        const double axis_safe_policy_cases_difference =
            max_difference_from_axis_safe_inverse_star_policy_cases();
        passed &= print_result(
            "AxisSafeInverseStar1Policy matches legacy branch cases exactly",
            axis_safe_policy_cases_difference == 0.0,
            max_detail(axis_safe_policy_cases_difference));
        const double inverse_star_policy_difference =
            max_difference_from_legacy_axis_safe_inverse_star_stencil(
                *field.grd, field, geometry_ids);
        passed &= print_result(
            "MPCNS inverse star-1 axis-safety policy matches legacy stencil exactly",
            inverse_star_policy_difference == 0.0,
            max_detail(inverse_star_policy_difference));


        const double geometry_snapshot_difference =
            max_difference_from_geometry_snapshot(*field.grd, field, geometry_ids);
        passed &= print_result("GeometrySnapshot reads exact legacy geometry",
                               geometry_snapshot_difference == 0.0,
                               max_detail(geometry_snapshot_difference));

        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &b_xi = field.field("B_xi", ib);
            FieldBlock &b_eta = field.field("B_eta", ib);
            FieldBlock &b_zeta = field.field("B_zeta", ib);
            FieldBlock &work_xi = field.field("dual_edge_work_xi", ib);
            FieldBlock &work_eta = field.field("dual_edge_work_eta", ib);
            FieldBlock &work_zeta = field.field("dual_edge_work_zeta", ib);
            FieldBlock &j_xi = field.field("E_xi", ib);
            FieldBlock &j_eta = field.field("E_eta", ib);
            FieldBlock &j_zeta = field.field("E_zeta", ib);

            if (!b_xi.is_allocated() || !b_eta.is_allocated() || !b_zeta.is_allocated() ||
                !work_xi.is_allocated() || !work_eta.is_allocated() || !work_zeta.is_allocated() ||
                !j_xi.is_allocated() || !j_eta.is_allocated() || !j_zeta.is_allocated())
                continue;

            const mpcns::discrete::LumpedHodge hodge{
                mpcns::compat::bind_legacy_lumped_hodge_coefficients(field, ib)};
            const mpcns::state::PrimalFaceFormView magnetic_flux{
                b_xi, b_eta, b_zeta};
            const mpcns::state::DualEdgeFormView dual_edge_work{
                work_xi, work_eta, work_zeta};
            const mpcns::state::DualFaceFormView dual_current{
                j_xi, j_eta, j_zeta};
            const mpcns::state::PrimalEdgeFormView electric_current{
                j_xi, j_eta, j_zeta};

            hodge.apply_star_2(magnetic_flux, dual_edge_work);
            mpcns::discrete::apply_d1_transpose(dual_edge_work, dual_current);
            hodge.apply_inverse_star_1(dual_current, electric_current);
        }
        apply_legacy_current_from_magnetic_flux(field, geometry_ids);
        const double hodge_current_difference = std::max(
            max_difference_inner(field, "E_xi", "J_legacy_xi"),
            std::max(max_difference_inner(field, "E_eta", "J_legacy_eta"),
                     max_difference_inner(field, "E_zeta", "J_legacy_zeta")));
        passed &= print_result("MPCNS lumped Hodge and D1 transpose match legacy current CT exactly",
                               hodge_current_difference == 0.0,
                               max_detail(hodge_current_difference));

        return passed;
    }
}
