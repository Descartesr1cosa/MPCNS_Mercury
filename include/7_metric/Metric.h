#pragma once

#include <cstdint>

class Field;
class Grid;

namespace METRIC
{
    struct MetricDiagnostics
    {
        std::int64_t jac_nonpositive = 0;
        std::int64_t area_nonpositive = 0;
        std::int64_t dl_nonpositive = 0;
        std::int64_t Hodge_star_inverse_2form_to_1form_nonfinite = 0;
        std::int64_t Hodge_star_2form_to_1form_nonfinite = 0;
        std::int64_t near_axis_singular = 0;
        std::int64_t near_axis_capped = 0;
    };

    void register_metric_fields(Field &fields, int geometry_ghost);
    void compute_metric_fields(Field &fields, Grid &grid);

    // Backward-compatible one-shot entry: register metric fields, then compute them.
    void build_field_geometry(Field &fields, Grid &grid, int geometry_ghost);

    MetricDiagnostics diagnose_metric_fields(Field &fields);
}
