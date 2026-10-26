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

#include <algorithm>
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
