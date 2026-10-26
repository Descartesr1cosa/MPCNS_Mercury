#include "7_metric/Metric.h"

#include "mpcns/compat/legacy_field_reconstruction.hpp"
#include "mpcns/compat/legacy_geometry_field_registration.hpp"
#include "mpcns/compat/legacy_lumped_metric_publisher.hpp"
#include "mpcns/compat/legacy_metric_diagnostics.hpp"
#include "mpcns/compat/legacy_primal_metric_publisher.hpp"
#include "mpcns/compat/legacy_reconstruction_metric_publisher.hpp"

namespace METRIC
{
void register_metric_fields(Field &fields, int geometry_ghost)
{
    mpcns::compat::register_legacy_geometry_fields(fields, geometry_ghost);
}

void compute_metric_fields(Field &fields, Grid &grid)
{
    // Dependency order is intentional; compatibility publishers retain legacy storage.
    // Do not reorder without an exact field-by-field regression gate.
    mpcns::compat::populate_legacy_primal_face_and_cell_metrics(fields, grid);
    mpcns::compat::populate_legacy_lumped_face_and_edge_metrics(fields, grid);
    mpcns::compat::populate_legacy_reconstruction_cache(fields, grid);
    mpcns::compat::populate_legacy_primal_edge_vectors(fields, grid);
    mpcns::compat::populate_legacy_diagonal_hodge_closure(fields, grid);
}

void build_field_geometry(Field &fields, Grid &grid, int geometry_ghost)
{
    register_metric_fields(fields, geometry_ghost);
    compute_metric_fields(fields, grid);
}

MetricDiagnostics diagnose_metric_fields(Field &fields)
{
    const auto legacy = mpcns::compat::diagnose_legacy_metric_fields(fields);
    return MetricDiagnostics{legacy.jac_nonpositive, legacy.area_nonpositive,
                             legacy.dl_nonpositive, legacy.inverse_star_1_nonfinite,
                             legacy.star_2_nonfinite, legacy.near_axis_singular,
                             legacy.near_axis_capped};
}

bool reconstruct_face_2form_to_cell(Field &fields, int fid_xi, int fid_eta, int fid_zeta, int iblock, int i, int j, int k, double out_xyz[3])
{
    return mpcns::compat::reconstruct_legacy_face_2form_to_cell(fields, fid_xi, fid_eta, fid_zeta, iblock, i, j, k, out_xyz);
}

bool reconstruct_edge_1form_to_cell(Field &fields, int fid_xi, int fid_eta, int fid_zeta, int iblock, int i, int j, int k, double out_xyz[3])
{
    return mpcns::compat::reconstruct_legacy_edge_1form_to_cell(fields, fid_xi, fid_eta, fid_zeta, iblock, i, j, k, out_xyz);
}

bool reconstruct_face_2form_to_node(Field &fields, int fid_xi, int fid_eta, int fid_zeta, int iblock, int i, int j, int k, double out_xyz[3])
{
    return mpcns::compat::reconstruct_legacy_face_2form_to_node(fields, fid_xi, fid_eta, fid_zeta, iblock, i, j, k, out_xyz);
}

bool reconstruct_edge_1form_to_node(Field &fields, int fid_xi, int fid_eta, int fid_zeta, int iblock, int i, int j, int k, double out_xyz[3])
{
    return mpcns::compat::reconstruct_legacy_edge_1form_to_node(fields, fid_xi, fid_eta, fid_zeta, iblock, i, j, k, out_xyz);
}
} // namespace METRIC
