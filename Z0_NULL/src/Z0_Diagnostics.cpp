#include "Z0_Diagnostics.h"

#include "3_field/Field.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace Z0
{
    FieldSummary SummarizeField(Field &field, const std::string &field_name)
    {
        FieldSummary s;
        s.min_value = std::numeric_limits<double>::infinity();
        s.max_value = -std::numeric_limits<double>::infinity();

        const int fid = field.field_id(field_name);
        const FieldDescriptor &desc = field.descriptor(fid);
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &fb = field.field(fid, ib);
            if (!fb.is_allocated())
                continue;
            const Int3 lo = fb.get_lo();
            const Int3 hi = fb.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                        {
                            const double v = fb(i, j, k, m);
                            if (std::isnan(v))
                            {
                                ++s.nan_count;
                                continue;
                            }
                            if (std::isinf(v))
                            {
                                ++s.inf_count;
                                continue;
                            }
                            if (v < -0.5e300)
                                ++s.sentinel_count;
                            ++s.finite_count;
                            s.min_value = std::min(s.min_value, v);
                            s.max_value = std::max(s.max_value, v);
                        }
        }
        if (s.finite_count == 0)
        {
            s.min_value = 0.0;
            s.max_value = 0.0;
        }
        return s;
    }

    void PrintStepDiagnostics(Field &field,
                              int step,
                              double time,
                              double dt,
                              std::ostream &os)
    {
        os << "[Z0_NULL] step=" << step << " time=" << time << " dt=" << dt << "\n";
        for (const std::string &name : {"null_phi"})
        {
            if (!field.has_field(name))
                continue;
            const FieldSummary s = SummarizeField(field, name);
            os << "  field=" << name
               << " min=" << s.min_value
               << " max=" << s.max_value
               << " nan=" << s.nan_count
               << " inf=" << s.inf_count
               << " sentinel=" << s.sentinel_count
               << "\n";
        }
    }
}
