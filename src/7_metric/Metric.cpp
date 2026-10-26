#include "7_metric/Metric.h"

#include "mpcns/compat/legacy_field_reconstruction.hpp"
#include "mpcns/compat/legacy_geometry_build_pipeline.hpp"
#include "mpcns/compat/legacy_metric_diagnostics.hpp"

namespace METRIC
{
void register_metric_fields(Field &fields, int geometry_ghost)
{
    mpcns::compat::LegacyGeometryBuildPipeline{}.register_fields(
        fields, geometry_ghost);
}

void compute_metric_fields(Field &fields, Grid &grid)
{
    // The compatibility pipeline preserves the historical publication order.
    // Do not reorder without an exact field-by-field regression gate.
    const mpcns::compat::LegacyGeometryBuildRequest request{fields, grid};
    mpcns::compat::LegacyGeometryBuildPipeline{}.publish(request);
}

void build_field_geometry(Field &fields, Grid &grid, int geometry_ghost)
{
    const mpcns::compat::LegacyGeometryBuildRequest request{fields, grid};
    // Finalization validates the typed GeometryBuilder hand-off. The returned
    // snapshot borrows existing storage and is intentionally not retained here.
    mpcns::compat::LegacyGeometryBuildPipeline{}.build(
        request, geometry_ghost);
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
