#include "Z0_Initializer.h"

#include "0_basic/LayoutTraits.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"

namespace Z0
{
    double analytic_value(const FieldDescriptor &desc,
                          int fid,
                          int iblock,
                          int i,
                          int j,
                          int k,
                          int comp,
                          const InitContext &ctx)
    {
        (void)desc;
        return analytic_value_for_rank(fid, comp, ctx.my_rank, iblock, i, j, k);
    }

    double analytic_value_for_rank(int fid,
                                   int comp,
                                   int rank,
                                   int iblock,
                                   int i,
                                   int j,
                                   int k)
    {
        return static_cast<double>(fid) * 1.0e8 +
               static_cast<double>(comp) * 1.0e7 +
               static_cast<double>(rank) * 1.0e5 +
               static_cast<double>(iblock) * 1.0e4 +
               static_cast<double>(i) * 100.0 +
               static_cast<double>(j) * 10.0 +
               static_cast<double>(k);
    }

    Box3 owned_box_for_location(Grid &grid,
                                int iblock,
                                StaggerLocation loc,
                                int dimension)
    {
        const Block &blk = grid.grids(iblock);
        Int3 ncells{blk.mx, blk.my, blk.mz};
        if (dimension < 3)
            ncells.k = 1;
        return LAYOUT::owned_box_from_cells(ncells, loc);
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
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                            fb(i, j, k, m) = ctx.sentinel;

            const Box3 owned = owned_box_for_location(*fields.grd, ib, desc.location, ctx.dimension);
            for (int i = owned.lo.i; i < owned.hi.i; ++i)
                for (int j = owned.lo.j; j < owned.hi.j; ++j)
                    for (int k = owned.lo.k; k < owned.hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                            fb(i, j, k, m) = fn(desc, fid, ib, i, j, k, m, ctx);
        }
    }

    void initialize_all_fields(Field &fields,
                               const InitContext &ctx,
                               InitFunction fn)
    {
        for (const FieldDescriptor &desc : fields.descriptors())
            initialize_field(fields, desc.name, ctx, fn);
    }
}
