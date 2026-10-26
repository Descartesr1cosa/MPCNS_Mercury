#include "Z0_Tests.h"

#include "Z0_Boundary.h"
#include "Z0_TestCommon.h"
#include "7_metric/Metric.h"
#include "8_dec/DecOps.h"
#include "mpcns/compat/legacy_geometry_field_ids.hpp"
#include "mpcns/compat/legacy_geometry_field_registration.hpp"
#include "mpcns/compat/legacy_geometry_build_pipeline.hpp"
#include "mpcns/compat/legacy_geometry_snapshot.hpp"
#include "mpcns/compat/legacy_field_reconstruction.hpp"
#include "mpcns/compat/legacy_metric_diagnostics.hpp"
#include "mpcns/compat/legacy_state_layout.hpp"
#include "mpcns/compat/legacy_state_view_adapter.hpp"
#include "mpcns/discrete/incidence_operator.hpp"
#include "mpcns/operators/electromagnetic_operators.hpp"
#include "mpcns/operators/faraday_operator.hpp"
#include "mpcns/state/simulation_state.hpp"

#include "mpcns/discrete/lumped_hodge.hpp"
#include "mpcns/discrete/dual_face_inverse_hodge.hpp"
#include "mpcns/geometry/singular_edge_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>
#include <utility>
#include <type_traits>
namespace Z0_METRIC_REFERENCE
{
    struct MetricPublisherSnapshot
    {
        const char *name;
        int block_index;
        Int3 lo;
        Int3 hi;
        int component_count;
        std::vector<double> values;
    };

    std::vector<MetricPublisherSnapshot>
    snapshot_metric_publisher_outputs(Field &field)
    {
        static constexpr std::array<const char *, 15> field_names{{
            "a_xi",
            "a_eta",
            "a_zeta",
            "pinvGT_cell",
            "pinvGT_xi",
            "pinvAT_xi",
            "pinvGT_eta",
            "pinvAT_eta",
            "pinvGT_zeta",
            "pinvAT_zeta",
            "Hodge_star_0form_to_3form_cell_lumped",
            "Hodge_star_3form_to_0form_cell_lumped",
            "Hodge_star_1form_to_2form_edge_xi_lumped",
            "Hodge_star_1form_to_2form_edge_eta_lumped",
            "Hodge_star_1form_to_2form_edge_zeta_lumped"}};

        std::vector<MetricPublisherSnapshot> snapshots;
        for (const char *name : field_names)
        {
            const int field_id = field.field_id(name);
            const FieldDescriptor &descriptor = field.descriptor(field_id);
            for (int ib = 0; ib < field.num_blocks(); ++ib)
            {
                FieldBlock &block = field.field(field_id, ib);
                MetricPublisherSnapshot snapshot{
                    name, ib, block.get_lo(), block.get_hi(), descriptor.ncomp, {}};
                if (!block.is_allocated())
                {
                    snapshots.push_back(std::move(snapshot));
                    continue;
                }

                const std::size_t count =
                    static_cast<std::size_t>(snapshot.hi.i - snapshot.lo.i) *
                    static_cast<std::size_t>(snapshot.hi.j - snapshot.lo.j) *
                    static_cast<std::size_t>(snapshot.hi.k - snapshot.lo.k) *
                    static_cast<std::size_t>(snapshot.component_count);
                snapshot.values.reserve(count);
                for (int i = snapshot.lo.i; i < snapshot.hi.i; ++i)
                    for (int j = snapshot.lo.j; j < snapshot.hi.j; ++j)
                        for (int k = snapshot.lo.k; k < snapshot.hi.k; ++k)
                            for (int m = 0; m < snapshot.component_count; ++m)
                                snapshot.values.push_back(block(i, j, k, m));
                snapshots.push_back(std::move(snapshot));
            }
        }
        return snapshots;
    }

    double max_difference_from_metric_publisher_snapshot(
        Field &field,
        const std::vector<MetricPublisherSnapshot> &snapshots)
    {
        double max_difference = 0.0;
        for (const MetricPublisherSnapshot &snapshot : snapshots)
        {
            FieldBlock &block = field.field(snapshot.name, snapshot.block_index);
            if (!block.is_allocated() ||
                block.get_lo().i != snapshot.lo.i ||
                block.get_lo().j != snapshot.lo.j ||
                block.get_lo().k != snapshot.lo.k ||
                block.get_hi().i != snapshot.hi.i ||
                block.get_hi().j != snapshot.hi.j ||
                block.get_hi().k != snapshot.hi.k ||
                block.descriptor().ncomp != snapshot.component_count)
            {
                max_difference = std::numeric_limits<double>::infinity();
                continue;
            }

            std::size_t offset = 0;
            for (int i = snapshot.lo.i; i < snapshot.hi.i; ++i)
                for (int j = snapshot.lo.j; j < snapshot.hi.j; ++j)
                    for (int k = snapshot.lo.k; k < snapshot.hi.k; ++k)
                        for (int m = 0; m < snapshot.component_count; ++m)
                        {
                            const double expected = snapshot.values[offset++];
                            const double actual = block(i, j, k, m);
                            if (!std::isfinite(expected) || !std::isfinite(actual))
                                max_difference = std::numeric_limits<double>::infinity();
                            else
                                max_difference = std::max(
                                    max_difference, std::abs(expected - actual));
                        }
        }
        return Z0_TEST::global_max(max_difference);
    }

    void populate_reconstruction_and_diagonal_hodge(Field &fields, Grid &grid);
}


namespace
{
    double max_difference_from_metric_compatibility_facade(Field &field)
    {
        using Reconstruction =
            bool (*)(Field &, int, int, int, int, int, int, int, double[3]);

        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };
        const auto compare = [&](Reconstruction facade,
                                 Reconstruction implementation,
                                 int fid_xi,
                                 int fid_eta,
                                 int fid_zeta,
                                 int iblock,
                                 int i,
                                 int j,
                                 int k)
        {
            double facade_value[3] = {0.0, 0.0, 0.0};
            double implementation_value[3] = {0.0, 0.0, 0.0};
            const bool facade_ok =
                facade(field, fid_xi, fid_eta, fid_zeta, iblock, i, j, k,
                       facade_value);
            const bool implementation_ok =
                implementation(field, fid_xi, fid_eta, fid_zeta, iblock, i,
                               j, k, implementation_value);
            if (facade_ok != implementation_ok)
            {
                max_difference = std::numeric_limits<double>::infinity();
                return;
            }
            if (facade_ok)
                for (int component = 0; component < 3; ++component)
                    update(facade_value[component],
                           implementation_value[component]);
        };

        const int b_xi = field.field_id("B_xi");
        const int b_eta = field.field_id("B_eta");
        const int b_zeta = field.field_id("B_zeta");
        const int e_xi = field.field_id("E_xi");
        const int e_eta = field.field_id("E_eta");
        const int e_zeta = field.field_id("E_zeta");

        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &cell = field.field("Jac", ib);
            if (!cell.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                continue;
            }
            const Int3 lo = cell.inner_lo();
            const Int3 hi = cell.inner_hi();
            const int i = (lo.i + hi.i - 1) / 2;
            const int j = (lo.j + hi.j - 1) / 2;
            const int k = (lo.k + hi.k - 1) / 2;

            compare(METRIC::reconstruct_face_2form_to_cell,
                    mpcns::compat::reconstruct_legacy_face_2form_to_cell,
                    b_xi, b_eta, b_zeta, ib, i, j, k);
            compare(METRIC::reconstruct_edge_1form_to_cell,
                    mpcns::compat::reconstruct_legacy_edge_1form_to_cell,
                    e_xi, e_eta, e_zeta, ib, i, j, k);
            compare(METRIC::reconstruct_face_2form_to_node,
                    mpcns::compat::reconstruct_legacy_face_2form_to_node,
                    b_xi, b_eta, b_zeta, ib, i, j, k);
            compare(METRIC::reconstruct_edge_1form_to_node,
                    mpcns::compat::reconstruct_legacy_edge_1form_to_node,
                    e_xi, e_eta, e_zeta, ib, i, j, k);
        }

        const auto facade_diagnostics = METRIC::diagnose_metric_fields(field);
        const auto implementation_diagnostics =
            mpcns::compat::diagnose_legacy_metric_fields(field);
        update(facade_diagnostics.jac_nonpositive,
               implementation_diagnostics.jac_nonpositive);
        update(facade_diagnostics.area_nonpositive,
               implementation_diagnostics.area_nonpositive);
        update(facade_diagnostics.dl_nonpositive,
               implementation_diagnostics.dl_nonpositive);
        update(facade_diagnostics.Hodge_star_inverse_2form_to_1form_nonfinite,
               implementation_diagnostics.inverse_star_1_nonfinite);
        update(facade_diagnostics.Hodge_star_2form_to_1form_nonfinite,
               implementation_diagnostics.star_2_nonfinite);
        update(facade_diagnostics.near_axis_singular,
               implementation_diagnostics.near_axis_singular);
        update(facade_diagnostics.near_axis_capped,
               implementation_diagnostics.near_axis_capped);
        return Z0_TEST::global_max(max_difference);
    }

    struct SingularEdgeClosureReference
    {
        double dual_area = 0.0;
        double inverse_hodge = 0.0;
        std::array<double, 3> canonical_edge_vector{{0.0, 0.0, 0.0}};
        bool used_incident_polygon = false;
    };

    SingularEdgeClosureReference legacy_singular_edge_closure_reference(
        const mpcns::geometry::SingularEdgeGeometryInput &input)
    {
        SingularEdgeClosureReference result;
        if (input.tangent_sample_count > 0)
            for (int component = 0; component < 3; ++component)
                result.canonical_edge_vector[component] =
                    input.tangent_sum[component] / input.tangent_sample_count;

        if (input.incident_dual_centers.size() >= 3)
        {
            auto points = input.incident_dual_centers;
            std::array<double, 3> center{{0.0, 0.0, 0.0}};
            for (const auto &point : points)
                for (int component = 0; component < 3; ++component)
                    center[component] += point[component];
            for (int component = 0; component < 3; ++component)
                center[component] /= points.size();

            std::array<double, 3> tangent = input.tangent_sum;
            const auto norm = [](const std::array<double, 3> &vector)
            {
                return std::sqrt(vector[0] * vector[0] +
                                 vector[1] * vector[1] +
                                 vector[2] * vector[2]);
            };
            const double tangent_norm = norm(tangent);
            if (tangent_norm > 0.0)
                for (double &component : tangent)
                    component /= tangent_norm;

            const std::array<double, 3> reference =
                std::abs(tangent[0]) < 0.9
                    ? std::array<double, 3>{{1.0, 0.0, 0.0}}
                    : std::array<double, 3>{{0.0, 1.0, 0.0}};
            std::array<double, 3> u{{
                tangent[1] * reference[2] - tangent[2] * reference[1],
                tangent[2] * reference[0] - tangent[0] * reference[2],
                tangent[0] * reference[1] - tangent[1] * reference[0]}};
            const double u_norm = norm(u);
            if (u_norm > 0.0)
                for (double &component : u)
                    component /= u_norm;
            const std::array<double, 3> v{{
                tangent[1] * u[2] - tangent[2] * u[1],
                tangent[2] * u[0] - tangent[0] * u[2],
                tangent[0] * u[1] - tangent[1] * u[0]}};

            std::sort(points.begin(), points.end(),
                      [&](const std::array<double, 3> &lhs,
                          const std::array<double, 3> &rhs)
                      {
                          const auto angle =
                              [&](const std::array<double, 3> &point)
                          {
                              double du = 0.0;
                              double dv = 0.0;
                              for (int component = 0; component < 3;
                                   ++component)
                              {
                                  du += (point[component] - center[component]) *
                                        u[component];
                                  dv += (point[component] - center[component]) *
                                        v[component];
                              }
                              return std::atan2(dv, du);
                          };
                          return angle(lhs) < angle(rhs);
                      });

            std::array<double, 3> area_vector{{0.0, 0.0, 0.0}};
            for (std::size_t index = 0; index < points.size(); ++index)
            {
                const auto &lhs = points[index];
                const auto &rhs = points[(index + 1) % points.size()];
                area_vector[0] += lhs[1] * rhs[2] - lhs[2] * rhs[1];
                area_vector[1] += lhs[2] * rhs[0] - lhs[0] * rhs[2];
                area_vector[2] += lhs[0] * rhs[1] - lhs[1] * rhs[0];
            }
            result.dual_area = 0.5 * norm(area_vector);
            result.used_incident_polygon = result.dual_area > 0.0;
        }

        if (!(result.dual_area > 0.0))
            result.dual_area = input.fallback_dual_area;
        result.inverse_hodge =
            result.dual_area > 0.0 ? input.primal_length / result.dual_area
                                   : 0.0;
        return result;
    }

    double max_difference_from_singular_edge_geometry_closure()
    {
        using Input = mpcns::geometry::SingularEdgeGeometryInput;
        using Vector = mpcns::geometry::CartesianVector;

        std::vector<Input> cases;
        Input polygon;
        polygon.primal_length = 4.0;
        polygon.tangent_sum = {{0.0, 0.0, 2.0}};
        polygon.tangent_sample_count = 2;
        polygon.incident_dual_centers = {
            Vector{{0.0, 1.0, 0.0}},
            Vector{{1.0, 0.0, 0.0}},
            Vector{{0.0, -1.0, 0.0}},
            Vector{{-1.0, 0.0, 0.0}}};
        polygon.fallback_dual_area = 9.0;
        cases.push_back(polygon);

        Input insufficient;
        insufficient.primal_length = 6.0;
        insufficient.tangent_sum = {{0.0, 0.0, 4.0}};
        insufficient.tangent_sample_count = 2;
        insufficient.incident_dual_centers = {
            Vector{{0.0, 0.0, 0.0}},
            Vector{{1.0, 0.0, 0.0}}};
        insufficient.fallback_dual_area = 3.0;
        cases.push_back(insufficient);

        Input degenerate;
        degenerate.primal_length = 5.0;
        degenerate.tangent_sum = {{0.0, 0.0, 3.0}};
        degenerate.tangent_sample_count = 3;
        degenerate.incident_dual_centers = {
            Vector{{-1.0, 0.0, 0.0}},
            Vector{{0.0, 0.0, 0.0}},
            Vector{{1.0, 0.0, 0.0}}};
        degenerate.fallback_dual_area = 4.0;
        cases.push_back(degenerate);

        double max_difference = 0.0;
        const auto update = [&](double expected, double actual)
        {
            if (!std::isfinite(expected) || !std::isfinite(actual))
                max_difference = std::numeric_limits<double>::infinity();
            else
                max_difference = std::max(max_difference,
                                          std::abs(expected - actual));
        };

        for (const Input &input : cases)
        {
            const SingularEdgeClosureReference expected =
                legacy_singular_edge_closure_reference(input);
            const mpcns::geometry::SingularEdgeGeometry actual =
                mpcns::geometry::close_singular_edge_geometry(input);
            update(expected.dual_area, actual.dual_area);
            update(expected.inverse_hodge, actual.inverse_hodge);
            for (int component = 0; component < 3; ++component)
                update(expected.canonical_edge_vector[component],
                       actual.canonical_edge_vector[component]);
            if (expected.used_incident_polygon != actual.used_incident_polygon)
                max_difference = std::numeric_limits<double>::infinity();
        }
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_legacy_primal_edge_vector_stencil(
        Grid &grid,
        Field &field,
        const mpcns::compat::LegacyGeometryFieldIds &geometry_ids)
    {
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
        const auto compare_direction = [&](const Block &block,
                                           FieldBlock &actual,
                                           const Offset &endpoint_offset)
        {
            if (!actual.is_allocated())
            {
                max_difference = std::numeric_limits<double>::infinity();
                return;
            }

            const Int3 lo = actual.get_lo();
            const Int3 hi = actual.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double expected[3] = {
                            block.x(i + endpoint_offset[0],
                                    j + endpoint_offset[1],
                                    k + endpoint_offset[2]) - block.x(i, j, k),
                            block.y(i + endpoint_offset[0],
                                    j + endpoint_offset[1],
                                    k + endpoint_offset[2]) - block.y(i, j, k),
                            block.z(i + endpoint_offset[0],
                                    j + endpoint_offset[1],
                                    k + endpoint_offset[2]) - block.z(i, j, k)};
                        for (int component = 0; component < 3; ++component)
                            update(expected[component],
                                   actual(i, j, k, component));
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
                field.field(geometry_ids.primal_edge_vector.xi, ib),
                xi_endpoint);
            compare_direction(
                block,
                field.field(geometry_ids.primal_edge_vector.eta, ib),
                eta_endpoint);
            compare_direction(
                block,
                field.field(geometry_ids.primal_edge_vector.zeta, ib),
                zeta_endpoint);
        }
        return Z0_TEST::global_max(max_difference);
    }

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
        const mpcns::geometry::GeometryBuildResult build_result =
            mpcns::compat::LegacyGeometryBuildPipeline{}.finalize(
                mpcns::compat::LegacyGeometryBuildRequest{field, grid});
        const mpcns::geometry::GeometrySnapshot &snapshot =
            build_result.snapshot;
        double max_difference =
            snapshot.is_complete() &&
                    build_result.completed_block_count ==
                        static_cast<std::size_t>(grid.nblock)
                ? 0.0
                : std::numeric_limits<double>::infinity();

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
            compare_directional(geometry.lumped_hodge.star_primal_2_to_dual_1,
                                geometry_ids.lumped_star_primal_2_to_dual_1, 1, ib);
            compare_directional(geometry.lumped_hodge.inverse_star_primal_1,
                                geometry_ids.lumped_inverse_star_primal_1, 1, ib);
            compare_field(geometry.reconstruction.cell_gradient_pseudoinverse,
                          field.field(geometry_ids.cell_reconstruction_gradient_pseudoinverse, ib), 9);
            compare_directional(geometry.reconstruction.edge_gradient_pseudoinverse,
                                geometry_ids.edge_reconstruction_gradient_pseudoinverse, 9, ib);
            compare_directional(geometry.reconstruction.edge_area_pseudoinverse,
                                geometry_ids.edge_reconstruction_area_pseudoinverse, 9, ib);
        }
        return Z0_TEST::global_max(max_difference);
    }

    double max_difference_from_state_layout(const Field &field)
    {
        const mpcns::compat::LegacyStateLayoutBinding binding =
            mpcns::compat::borrow_legacy_state_layout(field);
        const mpcns::state::StateLayout &layout = binding.layout;
        if (layout.size() != static_cast<std::size_t>(field.num_fields()))
            return std::numeric_limits<double>::infinity();

        for (int field_id = 0; field_id < field.num_fields(); ++field_id)
        {
            const FieldDescriptor &legacy = field.descriptor(field_id);
            const auto quantity_id = layout.find(legacy.name);
            if (!quantity_id.is_valid() ||
                binding.legacy_field_id(quantity_id) != field_id)
                return std::numeric_limits<double>::infinity();

            const auto &quantity = layout.descriptor(quantity_id);
            const auto expected_kind = [&]()
            {
                switch (legacy.value_kind)
                {
                case FieldValueKind::Scalar: return mpcns::state::QuantityKind::Scalar;
                case FieldValueKind::CartesianVector: return mpcns::state::QuantityKind::CartesianVector;
                case FieldValueKind::ConservedVector: return mpcns::state::QuantityKind::ConservativeState;
                case FieldValueKind::EdgeCovariant1Form: return mpcns::state::QuantityKind::PrimalEdgeForm;
                case FieldValueKind::FaceContravariant2Form: return mpcns::state::QuantityKind::PrimalFaceForm;
                case FieldValueKind::GeometryMetric: return mpcns::state::QuantityKind::GeometryMetric;
                case FieldValueKind::CouplingBufferOnly:
                case FieldValueKind::Temporary: return mpcns::state::QuantityKind::Scalar;
                }
                return mpcns::state::QuantityKind::Scalar;
            }();
            const auto expected_role = [&]()
            {
                switch (legacy.storage_role)
                {
                case FieldStorageRole::State: break;
                case FieldStorageRole::Derived: return mpcns::state::StorageRole::Derived;
                case FieldStorageRole::Workspace: return mpcns::state::StorageRole::Workspace;
                case FieldStorageRole::CouplingBuffer: return mpcns::state::StorageRole::CouplingBuffer;
                }
                if (legacy.value_kind == FieldValueKind::Temporary)
                    return mpcns::state::StorageRole::Workspace;
                if (legacy.value_kind == FieldValueKind::CouplingBufferOnly)
                    return mpcns::state::StorageRole::CouplingBuffer;
                if (legacy.value_kind == FieldValueKind::GeometryMetric)
                    return mpcns::state::StorageRole::Derived;
                return mpcns::state::StorageRole::State;
            }();
            const auto expected_owner = [&]()
            {
                switch (legacy.sync.owner_sync)
                {
                case OwnerSyncPolicy::None: return mpcns::state::OwnerPolicy::None;
                case OwnerSyncPolicy::NodeOwner: return mpcns::state::OwnerPolicy::Node;
                case OwnerSyncPolicy::EdgeOwner: return mpcns::state::OwnerPolicy::Edge;
                case OwnerSyncPolicy::FaceOwner: return mpcns::state::OwnerPolicy::Face;
                }
                return mpcns::state::OwnerPolicy::None;
            }();
            if (quantity.canonical_name != legacy.name ||
                quantity.external_name != legacy.name ||
                quantity.kind != expected_kind ||
                quantity.storage_role != expected_role ||
                quantity.location != legacy.location ||
                quantity.component_count != legacy.ncomp ||
                quantity.ghost_layers != legacy.nghost ||
                quantity.allocation_domain != legacy.physics ||
                quantity.synchronization.group != legacy.sync.group ||
                quantity.synchronization.needs_coupling != legacy.sync.do_coupling ||
                quantity.synchronization.needs_physical_boundary != legacy.sync.do_physical ||
                quantity.synchronization.needs_halo != legacy.sync.do_halo ||
                quantity.synchronization.halo_level != legacy.sync.halo_level ||
                quantity.synchronization.owner_policy != expected_owner ||
                quantity.synchronization.orientation_aware != legacy.sync.orientation_aware)
                return std::numeric_limits<double>::infinity();
        }
        return 0.0;
    }

    double state_layout_rejection_difference()
    {
        const auto rejects = [](std::vector<mpcns::state::QuantityDescriptor> quantities)
        {
            try
            {
                (void)mpcns::state::StateLayout{std::move(quantities)};
            }
            catch (const std::invalid_argument &)
            {
                return true;
            }
            return false;
        };

        mpcns::state::QuantityDescriptor node;
        node.canonical_name = "node";
        node.external_name = "node";
        node.kind = mpcns::state::QuantityKind::PrimalZeroForm;
        node.location = StaggerLocation::Node;
        const bool duplicate_rejected = rejects({node, node});

        auto invalid_edge = node;
        invalid_edge.canonical_name = "bad_edge";
        invalid_edge.external_name = "bad_edge";
        invalid_edge.kind = mpcns::state::QuantityKind::PrimalEdgeForm;
        const bool incompatible_location_rejected = rejects({invalid_edge});

        const mpcns::state::StateLayout layout{{node}};
        bool invalid_id_rejected = false;
        try
        {
            (void)layout.descriptor(mpcns::state::QuantityId{99});
        }
        catch (const std::out_of_range &)
        {
            invalid_id_rejected = true;
        }
        return duplicate_rejected && incompatible_location_rejected && invalid_id_rejected
                   ? 0.0
                   : std::numeric_limits<double>::infinity();
    }

    double max_difference_from_legacy_state_view(Field &field)
    {
        // Bind a genuinely inactive production-style FieldBlock: its allocation
        // domain deliberately matches no Z0 block.  This is registered only
        // after the catalog/layout equivalence check above, and has no sync
        // contract or downstream consumer.
        constexpr const char *inactive_name = "state_view_test_inactive";
        if (!field.has_field(inactive_name))
        {
            FieldDescriptor inactive{inactive_name, StaggerLocation::Cell, 2, 1};
            inactive.physics = "__state_view_test_inactive_domain__";
            field.register_field(inactive);
        }
        const auto binding = mpcns::compat::borrow_legacy_state_layout(field);
        mpcns::compat::LegacyStateViewAdapter adapter{field, binding};
        static_assert(std::is_same<decltype(std::declval<mpcns::state::ConstStateView>()(0, 0, 0)),
                                   const double &>::value,
                      "ConstStateView must not provide mutable element access");

        const std::array<const char *, 8> names{{
            "U", "phi", "E_xi", "E_eta", "E_zeta", "B_xi", "B_eta", "B_zeta"}};
        double max_difference = 0.0;
        bool inactive_checked = false;
        for (const char *name : names)
        {
            const auto id = binding.layout.find(name);
            if (!id.is_valid())
                return std::numeric_limits<double>::infinity();
            for (int ib = 0; ib < field.num_blocks(); ++ib)
            {
                FieldBlock &legacy = field.field(binding.legacy_field_id(id), ib);
                if (!legacy.is_allocated())
                {
                    try { (void)adapter.view(id, ib); }
                    catch (const std::logic_error &) { inactive_checked = true; continue; }
                    return std::numeric_limits<double>::infinity();
                }
                const auto view = adapter.view(id, ib);
                const auto const_view = adapter.const_view(id, ib);
                const Int3 lo = legacy.get_lo();
                const Int3 hi = legacy.get_hi();
                if (view.extent().lo.i != lo.i || view.extent().lo.j != lo.j ||
                    view.extent().lo.k != lo.k || view.extent().hi.i != hi.i ||
                    view.extent().hi.j != hi.j || view.extent().hi.k != hi.k ||
                    view.component_count() != legacy.descriptor().ncomp)
                    return std::numeric_limits<double>::infinity();
                const Int3 inner_lo = legacy.inner_lo();
                const std::array<Int3, 3> samples{{
                    lo, inner_lo, Int3{hi.i - 1, hi.j - 1, hi.k - 1}}};
                const std::array<int, 2> components{{0, legacy.descriptor().ncomp - 1}};
                for (const Int3 &sample : samples)
                    for (const int component : components)
                    {
                        const double original = legacy(sample.i, sample.j, sample.k, component);
                        max_difference = std::max(max_difference,
                            std::abs(const_view.at(sample.i, sample.j, sample.k, component) - original));
                        view.at(sample.i, sample.j, sample.k, component) = original + 0.375;
                        max_difference = std::max(max_difference,
                            std::abs(legacy(sample.i, sample.j, sample.k, component) - (original + 0.375)));
                        legacy(sample.i, sample.j, sample.k, component) = original - 0.125;
                        max_difference = std::max(max_difference,
                            std::abs(view.at(sample.i, sample.j, sample.k, component) - (original - 0.125)));
                        legacy(sample.i, sample.j, sample.k, component) = original;
                    }
                bool out_of_range_rejected = false;
                try { (void)const_view.at(hi.i, lo.j, lo.k, components.front()); }
                catch (const std::out_of_range &) { out_of_range_rejected = true; }
                if (!out_of_range_rejected)
                    return std::numeric_limits<double>::infinity();
            }
        }
        const auto inactive_id = binding.layout.find(inactive_name);
        if (!inactive_id.is_valid())
            return std::numeric_limits<double>::infinity();
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            if (field.field(binding.legacy_field_id(inactive_id), ib).is_allocated())
                return std::numeric_limits<double>::infinity();
            try { (void)adapter.view(inactive_id, ib); }
            catch (const std::logic_error &) { inactive_checked = true; continue; }
            return std::numeric_limits<double>::infinity();
        }

        bool invalid_id_checked = false;
        try { (void)adapter.view(mpcns::state::QuantityId{99999}, 0); }
        catch (const std::out_of_range &) { invalid_id_checked = true; }
        return inactive_checked && invalid_id_checked
                   ? max_difference
                   : std::numeric_limits<double>::infinity();
    }

    double simulation_state_workspace_difference(Field &field)
    {
        constexpr const char *workspace_name = "state_view_test_workspace";
        if (!field.has_field(workspace_name))
        {
            // Jstar-shaped scratch: it remains a primal face 2-form while
            // the independent role marks its storage lifetime.
            FieldDescriptor workspace{workspace_name, StaggerLocation::FaceXi, 1, 1};
            workspace.value_kind = FieldValueKind::FaceContravariant2Form;
            workspace.storage_role = FieldStorageRole::Workspace;
            field.register_field(workspace);
        }

        const auto binding = mpcns::compat::borrow_legacy_state_layout(field);
        const auto quantity = binding.layout.find(workspace_name);
        if (!quantity.is_valid())
            return std::numeric_limits<double>::infinity();
        if (binding.layout.descriptor(quantity).kind != mpcns::state::QuantityKind::PrimalFaceForm ||
            binding.layout.descriptor(quantity).storage_role != mpcns::state::StorageRole::Workspace ||
            field.num_blocks() < 2)
            return std::numeric_limits<double>::infinity();

        mpcns::state::SimulationState state{binding.layout};
        static_assert(!std::is_copy_constructible<mpcns::state::SimulationState>::value &&
                          !std::is_move_constructible<mpcns::state::SimulationState>::value,
                      "SimulationState identity must not be copied or moved");
        const auto handle = state.handle(quantity);
        mpcns::compat::LegacyStateViewAdapter adapter{field, binding};
        mpcns::state::Workspace workspace{state};
        double max_difference = 0.0;
        std::array<double, 2> original_values{{0.0, 0.0}};

        for (int ib = 0; ib < 2; ++ib)
        {
            FieldBlock &legacy = field.field(binding.legacy_field_id(quantity), ib);
            if (!legacy.is_allocated())
                return std::numeric_limits<double>::infinity();
            const auto borrowed = adapter.view(quantity, ib);
            if (borrowed.local_block_id() != ib ||
                &borrowed(legacy.get_lo().i, legacy.get_lo().j, legacy.get_lo().k, 0) !=
                    &legacy(legacy.get_lo().i, legacy.get_lo().j, legacy.get_lo().k, 0))
                return std::numeric_limits<double>::infinity();
            workspace.bind(handle, ib, borrowed);
            const Int3 lo = legacy.get_lo();
            const double original = legacy(lo.i, lo.j, lo.k, 0);
            original_values[static_cast<std::size_t>(ib)] = original;
            workspace.view(handle, ib)(lo.i, lo.j, lo.k, 0) = original + 0.625 + ib;
            max_difference = std::max(max_difference,
                std::abs(legacy(lo.i, lo.j, lo.k, 0) - (original + 0.625 + ib)));
        }

        // Rebinding block zero replaces only its view and leaves block one.
        workspace.bind(handle, 0, adapter.view(quantity, 0));
        if (!workspace.is_bound(handle, 0) || !workspace.is_bound(handle, 1))
            return std::numeric_limits<double>::infinity();
        bool swapped_block_rejected = false;
        try { workspace.bind(handle, 0, adapter.view(quantity, 1)); }
        catch (const std::invalid_argument &) { swapped_block_rejected = true; }
        for (int ib = 0; ib < 2; ++ib)
        {
            FieldBlock &legacy = field.field(binding.legacy_field_id(quantity), ib);
            const Int3 lo = legacy.get_lo();
            max_difference = std::max(max_difference,
                std::abs(legacy(lo.i, lo.j, lo.k, 0) -
                         (original_values[static_cast<std::size_t>(ib)] + 0.625 + ib)));
        }

        bool non_workspace_rejected = false;
        try { workspace.bind(state.handle(binding.layout.find("B_xi")), 0,
                             adapter.view(binding.layout.find("B_xi"), 0)); }
        catch (const std::invalid_argument &) { non_workspace_rejected = true; }
        bool component_mismatch_rejected = false;
        const auto valid_view = adapter.view(quantity, 0);
        try
        {
            double one_value = 0.0;
            mpcns::state::StateView wrong_components{
                &one_value, valid_view.extent(), valid_view.inner_extent(), 2,
                1, 1, 1, 1, valid_view.location(), quantity, 0};
            workspace.bind(handle, 0, wrong_components);
        }
        catch (const std::invalid_argument &) { component_mismatch_rejected = true; }
        bool foreign_handle_rejected = false;
        mpcns::state::SimulationState foreign_state{binding.layout};
        try { workspace.bind(foreign_state.handle(quantity), 0, valid_view); }
        catch (const std::out_of_range &) { foreign_handle_rejected = true; }

        workspace.release(handle, 0);
        if (workspace.is_bound(handle, 0) || !workspace.is_bound(handle, 1))
            return std::numeric_limits<double>::infinity();
        workspace.clear();
        bool unbound_rejected = false;
        try { (void)workspace.view(handle, 1); }
        catch (const std::logic_error &) { unbound_rejected = true; }
        for (int ib = 0; ib < 2; ++ib)
        {
            FieldBlock &legacy = field.field(binding.legacy_field_id(quantity), ib);
            const Int3 lo = legacy.get_lo();
            legacy(lo.i, lo.j, lo.k, 0) = original_values[static_cast<std::size_t>(ib)];
        }
        return swapped_block_rejected && non_workspace_rejected && component_mismatch_rejected && foreign_handle_rejected &&
                       unbound_rejected && !workspace.is_bound(handle, 0) && !workspace.is_bound(handle, 1)
                   ? max_difference
                   : std::numeric_limits<double>::infinity();
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
        const mpcns::geometry::GeometrySnapshot hodge_snapshot =
            mpcns::compat::LegacyGeometryBuildPipeline{}.finalize(
                mpcns::compat::LegacyGeometryBuildRequest{field, *field.grd}).snapshot;
        const auto state_binding = mpcns::compat::borrow_legacy_state_layout(field);
        mpcns::compat::LegacyStateViewAdapter state_adapter{field, state_binding};
        const auto e_xi_id = state_binding.layout.find("E_xi");
        const auto e_eta_id = state_binding.layout.find("E_eta");
        const auto e_zeta_id = state_binding.layout.find("E_zeta");
        const auto b_xi_id = state_binding.layout.find("B_xi");
        const auto b_eta_id = state_binding.layout.find("B_eta");
        const auto b_zeta_id = state_binding.layout.find("B_zeta");
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

            const mpcns::state::PrimalEdgeFormStateView edge{
                state_adapter.view(e_xi_id, ib), state_adapter.view(e_eta_id, ib),
                state_adapter.view(e_zeta_id, ib)};
            const mpcns::state::PrimalFaceFormStateView face{
                state_adapter.view(b_xi_id, ib), state_adapter.view(b_eta_id, ib),
                state_adapter.view(b_zeta_id, ib)};
            mpcns::operators::FaradayOperator::compute_rhs(edge, face);
            // The legacy CT oracle below computes +d1(E); Faraday RHS is
            // defined as -d1(E), so convert only for this exact oracle check.
            for (auto *block : {&face_xi, &face_eta, &face_zeta})
            {
                const Int3 lo = block->inner_lo();
                const Int3 hi = block->inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            (*block)(i, j, k, 0) = -(*block)(i, j, k, 0);
            }
        }
        const double d1_legacy_difference =
            max_difference_from_legacy_ct_d1(field, e, b);
        passed &= print_result("StateView Faraday RHS matches legacy CT stencil exactly",
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
        fill_edge_form(field);
        const double metric_compatibility_facade_difference =
            max_difference_from_metric_compatibility_facade(field);
        passed &= print_result(
            "METRIC reconstruction and diagnostics facade matches compat module exactly",
            metric_compatibility_facade_difference == 0.0,
            max_detail(metric_compatibility_facade_difference));

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
        const double primal_edge_vector_difference =
            max_difference_from_legacy_primal_edge_vector_stencil(
                *field.grd, field, geometry_ids);
        passed &= print_result(
            "MPCNS primal edge vectors match legacy stencil exactly",
            primal_edge_vector_difference == 0.0,
            max_detail(primal_edge_vector_difference));

        const double singular_edge_closure_difference =
            max_difference_from_singular_edge_geometry_closure();
        passed &= print_result(
            "MPCNS singular-edge closure matches legacy polygon/fallback reference",
            singular_edge_closure_difference == 0.0,
            max_detail(singular_edge_closure_difference));

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
        // Use the pre-Unit-16 algorithms as a test-only oracle. They write
        // the same Field backend, so this intentionally overwrites the relocated
        // outputs before the production coordinator rebuilds them.
        Z0_METRIC_REFERENCE::populate_reconstruction_and_diagonal_hodge(
            field, *field.grd);
        const auto legacy_metric_publisher_snapshot =
            Z0_METRIC_REFERENCE::snapshot_metric_publisher_outputs(field);
        METRIC::compute_metric_fields(field, *field.grd);
        const double metric_publisher_difference =
            Z0_METRIC_REFERENCE::max_difference_from_metric_publisher_snapshot(
                field, legacy_metric_publisher_snapshot);
        passed &= print_result(
            "MPCNS reconstruction cache and diagonal Hodge publishers match baseline exactly",
            metric_publisher_difference == 0.0,
            max_detail(metric_publisher_difference));



        const double geometry_snapshot_difference =
            max_difference_from_geometry_snapshot(*field.grd, field, geometry_ids);
        passed &= print_result("GeometryBuild pipeline reads exact legacy geometry",
                               geometry_snapshot_difference == 0.0,
                               max_detail(geometry_snapshot_difference));

        const double state_layout_difference = max_difference_from_state_layout(field);
        passed &= print_result("StateLayout exactly describes the legacy Field catalog",
                               state_layout_difference == 0.0,
                               max_detail(state_layout_difference));
        const double state_layout_rejection = state_layout_rejection_difference();
        passed &= print_result("StateLayout rejects invalid IDs, duplicate names, and incompatible forms",
                               state_layout_rejection == 0.0,
                               max_detail(state_layout_rejection));
        const double state_view_difference = max_difference_from_legacy_state_view(field);
        passed &= print_result("StateView is zero-copy, stride-correct, and safely bounded",
                               state_view_difference == 0.0,
                               max_detail(state_view_difference));
        const double simulation_state_workspace = simulation_state_workspace_difference(field);
        passed &= print_result("SimulationState Workspace is multi-block, zero-copy, and safely keyed",
                               simulation_state_workspace == 0.0,
                               max_detail(simulation_state_workspace));

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

            const mpcns::discrete::SnapshotLumpedHodge hodge{
                hodge_snapshot.block(static_cast<std::size_t>(ib)).lumped_hodge};
            const mpcns::operators::ElectromagneticOperators electromagnetic{hodge};
            const mpcns::state::PrimalFaceFormView magnetic_flux{
                b_xi, b_eta, b_zeta};
            const mpcns::state::DualEdgeFormView dual_edge_work{
                work_xi, work_eta, work_zeta};
            const mpcns::state::DualFaceFormView dual_current{
                j_xi, j_eta, j_zeta};
            const mpcns::state::PrimalEdgeFormView electric_current{
                j_xi, j_eta, j_zeta};
            const mpcns::operators::MagneticCurrentWorkspace workspace{
                dual_edge_work, dual_current};
            electromagnetic.compute_current(
                magnetic_flux, workspace, electric_current);

        }
        apply_legacy_current_from_magnetic_flux(field, geometry_ids);
        const double hodge_current_difference = std::max(
            max_difference_inner(field, "E_xi", "J_legacy_xi"),
            std::max(max_difference_inner(field, "E_eta", "J_legacy_eta"),
                     max_difference_inner(field, "E_zeta", "J_legacy_zeta")));
        passed &= print_result("Snapshot-backed MPCNS lumped Hodge and D1 transpose match legacy current CT exactly",
                               hodge_current_difference == 0.0,
                               max_detail(hodge_current_difference));

        return passed;
    }
}
