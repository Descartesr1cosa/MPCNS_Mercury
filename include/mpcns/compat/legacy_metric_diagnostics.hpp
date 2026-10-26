#pragma once

#include <cstdint>

class Field;

namespace mpcns::compat
{

struct LegacyMetricDiagnostics
{
    std::int64_t jac_nonpositive = 0;
    std::int64_t area_nonpositive = 0;
    std::int64_t dl_nonpositive = 0;
    std::int64_t inverse_star_1_nonfinite = 0;
    std::int64_t star_2_nonfinite = 0;
    std::int64_t near_axis_singular = 0;
    std::int64_t near_axis_capped = 0;
};

LegacyMetricDiagnostics diagnose_legacy_metric_fields(Field &fields);

} // namespace mpcns::compat
