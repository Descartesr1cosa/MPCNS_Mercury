#pragma once

#include <iosfwd>
#include <string>

class Field;

namespace Z0
{
    struct FieldSummary
    {
        double min_value = 0.0;
        double max_value = 0.0;
        long long nan_count = 0;
        long long inf_count = 0;
        long long sentinel_count = 0;
        long long finite_count = 0;
    };

    FieldSummary SummarizeField(Field &field, const std::string &field_name);
    void PrintStepDiagnostics(Field &field,
                              int step,
                              double time,
                              double dt,
                              std::ostream &os);
}
