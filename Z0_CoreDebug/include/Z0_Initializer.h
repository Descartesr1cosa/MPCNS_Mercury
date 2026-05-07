#pragma once

#include "0_basic/TYPES.h"
#include "3_field/Field_Type.h"

#include <functional>
#include <string>

class Field;
class Grid;

namespace Z0
{
    constexpr double Sentinel = -1.0e300;

    struct InitContext
    {
        int my_rank = 0;
        int dimension = 3;
        double sentinel = Sentinel;
    };

    using InitFunction = std::function<double(
        const FieldDescriptor &desc,
        int fid,
        int iblock,
        int i,
        int j,
        int k,
        int comp,
        const InitContext &ctx)>;

    double analytic_value(const FieldDescriptor &desc,
                          int fid,
                          int iblock,
                          int i,
                          int j,
                          int k,
                          int comp,
                          const InitContext &ctx);

    double analytic_value_for_rank(int fid,
                                   int comp,
                                   int rank,
                                   int iblock,
                                   int i,
                                   int j,
                                   int k);

    Box3 owned_box_for_location(Grid &grid,
                                int iblock,
                                StaggerLocation loc,
                                int dimension);

    void initialize_field(Field &fields,
                          const std::string &field_name,
                          const InitContext &ctx,
                          InitFunction fn = analytic_value);

    void initialize_all_fields(Field &fields,
                               const InitContext &ctx,
                               InitFunction fn = analytic_value);
}
