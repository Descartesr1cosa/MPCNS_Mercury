#pragma once

#include "3_field/Field_Type.h"

#include <functional>
#include <string>

class Field;

namespace Z0_NULL
{
    struct InitContext
    {
        int my_rank = 0;
        int dimension = 3;
        double sentinel = -9.87654321e200;
        bool fill_allocated_with_sentinel = true;
    };

    using InitFunction = std::function<double(
        const FieldDescriptor &desc,
        int fid,
        int iblock,
        int i,
        int j,
        int k,
        int m,
        const InitContext &ctx)>;

    // Default pattern:
    // fid*10000000 + m*1000000 + my_rank*10000 + iblock*1000 + i*100 + j*10 + k
    double default_init_value(const FieldDescriptor &desc,
                              int fid,
                              int iblock,
                              int i,
                              int j,
                              int k,
                              int m,
                              const InitContext &ctx);

    void initialize_field(Field &fields,
                          const std::string &field_name,
                          const InitContext &ctx,
                          InitFunction fn = default_init_value);

    void initialize_registered_fields(Field &fields,
                                      const InitContext &ctx,
                                      InitFunction fn = default_init_value);

    // Recommended template path for a new physics program: initialize only the
    // representative Z0_NULL physics fields, leaving grid/metric/framework
    // fields untouched.
    void initialize_null_physics_fields(Field &fields,
                                        const InitContext &ctx,
                                        InitFunction fn = default_init_value);

    // Strong framework test path: initializes every FieldCatalog entry,
    // including metric/framework fields. Use only for explicit stress tests.
    void initialize_all_fields(Field &fields,
                               const InitContext &ctx,
                               InitFunction fn = default_init_value);

    void initialize_null_fields(Field &fields);

    void initialize_halo_smoke_fields(Field &fields);
}
