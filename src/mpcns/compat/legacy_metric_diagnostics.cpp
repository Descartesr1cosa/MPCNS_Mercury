#include "mpcns/compat/legacy_metric_diagnostics.hpp"

#include "3_field/Field.h"

#include <cmath>

namespace mpcns::compat
{
LegacyMetricDiagnostics diagnose_legacy_metric_fields(Field &fields)
{
    LegacyMetricDiagnostics diag;

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
                    diag.inverse_star_1_nonfinite);
        scan_finite(fields.field("Hodge_star_inverse_2form_to_1form_edge_eta_lumped", ib),
                    diag.inverse_star_1_nonfinite);
        scan_finite(fields.field("Hodge_star_inverse_2form_to_1form_edge_zeta_lumped", ib),
                    diag.inverse_star_1_nonfinite);

        scan_finite(fields.field("Hodge_star_2form_to_1form_face_xi_lumped", ib),
                    diag.star_2_nonfinite);
        scan_finite(fields.field("Hodge_star_2form_to_1form_face_eta_lumped", ib),
                    diag.star_2_nonfinite);
        scan_finite(fields.field("Hodge_star_2form_to_1form_face_zeta_lumped", ib),
                    diag.star_2_nonfinite);

        scan_zero(fields.field("Astar_xi", ib), diag.near_axis_singular);
        scan_zero(fields.field("Astar_eta", ib), diag.near_axis_singular);
        scan_zero(fields.field("Astar_zeta", ib), diag.near_axis_singular);

        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_xi_lumped", ib), diag.near_axis_capped);
        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_eta_lumped", ib), diag.near_axis_capped);
        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_zeta_lumped", ib), diag.near_axis_capped);
    }

    return diag;
}
} // namespace mpcns::compat
