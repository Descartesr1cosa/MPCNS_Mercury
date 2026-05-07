#include "Z0_NullInitializer.h"

#include "0_basic/MPI_WRAPPER.h"
#include "3_field/2_MPCNS_Field.h"

namespace
{
    void fill_field(Field &fld, const std::string &field_name)
    {
        int rank = 0;
        PARALLEL::mpi_rank(&rank);

        const int fid = fld.field_id(field_name);
        const auto &desc = fld.descriptor(fid);

        for (int ib = 0; ib < fld.num_blocks(); ++ib)
        {
            FieldBlock &fb = fld.field(fid, ib);
            if (!fb.is_allocated())
                continue;

            const Int3 lo = fb.get_lo();
            const Int3 hi = fb.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                        {
                            fb(i, j, k, m) =
                                1.0e6 * rank +
                                1.0e4 * fid +
                                1.0e2 * ib +
                                1.0e-1 * i +
                                1.0e-2 * j +
                                1.0e-3 * k +
                                1.0e-4 * m;
                        }
        }
    }
}

namespace Z0_NULL
{
    void initialize_null_fields(Field &fields)
    {
        for (const auto &desc : fields.descriptors())
            fill_field(fields, desc.name);
    }

    void initialize_halo_smoke_fields(Field &fields)
    {
        for (const auto &request : fields.halo_requests())
            fill_field(fields, request.field_name);
    }
}
