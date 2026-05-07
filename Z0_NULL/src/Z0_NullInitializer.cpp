#include "Z0_NullInitializer.h"

#include "0_basic/MPI_WRAPPER.h"
#include "3_field/2_MPCNS_Field.h"

namespace Z0_NULL
{
    double default_init_value(const FieldDescriptor &desc,
                              int fid,
                              int iblock,
                              int i,
                              int j,
                              int k,
                              int m,
                              const InitContext &ctx)
    {
        (void)desc;

        return static_cast<double>(fid) * 10000000.0 +
               static_cast<double>(m) * 1000000.0 +
               static_cast<double>(ctx.my_rank) * 10000.0 +
               static_cast<double>(iblock) * 1000.0 +
               static_cast<double>(i) * 100.0 +
               static_cast<double>(j) * 10.0 +
               static_cast<double>(k);
    }

    void initialize_field(Field &fields,
                          const std::string &field_name,
                          const InitContext &ctx,
                          InitFunction fn)
    {
        const int fid = fields.field_id(field_name);
        const FieldDescriptor &desc = fields.descriptor(fid);

        for (int ib = 0; ib < fields.num_blocks(); ++ib)
        {
            FieldBlock &fb = fields.field(fid, ib);
            if (!fb.is_allocated())
                continue;

            const Int3 lo = fb.get_lo();
            const Int3 hi = fb.get_hi();

            if (ctx.fill_allocated_with_sentinel)
            {
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            for (int m = 0; m < desc.ncomp; ++m)
                                fb(i, j, k, m) = ctx.sentinel;
            }

            // TODO: once FieldBlock exposes owned_box/active_box, keep ghost
            // cells sentinel-filled and initialize only owned interior data.
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                            fb(i, j, k, m) = fn(desc, fid, ib, i, j, k, m, ctx);
        }
    }

    void initialize_registered_fields(Field &fields,
                                      const InitContext &ctx,
                                      InitFunction fn)
    {
        for (const auto &request : fields.halo_requests())
            initialize_field(fields, request.field_name, ctx, fn);
    }

    void initialize_all_fields(Field &fields,
                               const InitContext &ctx,
                               InitFunction fn)
    {
        for (const auto &desc : fields.descriptors())
            initialize_field(fields, desc.name, ctx, fn);
    }

    void initialize_null_fields(Field &fields)
    {
        InitContext ctx;
        PARALLEL::mpi_rank(&ctx.my_rank);
        initialize_all_fields(fields, ctx);
    }

    void initialize_halo_smoke_fields(Field &fields)
    {
        InitContext ctx;
        PARALLEL::mpi_rank(&ctx.my_rank);
        initialize_registered_fields(fields, ctx);
    }
}
