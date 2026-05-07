#include "Z1_Initial.h"

#include "Z1_Boundary.h"
#include "0_basic/LayoutTraits.h"
#include "3_field/2_MPCNS_Field.h"

namespace
{
    bool in_box(const Box3 &b, int i, int j, int k)
    {
        return i >= b.lo.i && i < b.hi.i &&
               j >= b.lo.j && j < b.hi.j &&
               k >= b.lo.k && k < b.hi.k;
    }

    double initial_value(const FieldDescriptor &desc, int comp)
    {
        if (desc.name == "rho")
            return 1.0;
        if (desc.name == "U")
            return comp == 0 ? 1.0 : 0.0;
        if (desc.name == "B_xi")
            return 0.1;
        return 0.0;
    }
}

namespace Z1
{
    void InitializeFields(Field &field, Param &param, Z1_Boundary &boundary)
    {
        (void)param;
        for (const FieldDescriptor &desc : field.descriptors())
        {
            const int fid = field.field_id(desc.name);
            for (int ib = 0; ib < field.num_blocks(); ++ib)
            {
                FieldBlock &fb = field.field(fid, ib);
                if (!fb.is_allocated())
                    continue;

                const Int3 lo = fb.get_lo();
                const Int3 hi = fb.get_hi();
                const Box3 owned = LAYOUT::owned_box_from_cells(
                    {field.grd->grids(ib).mx, field.grd->grids(ib).my, field.grd->grids(ib).mz},
                    desc.location);

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            for (int m = 0; m < desc.ncomp; ++m)
                                fb(i, j, k, m) = in_box(owned, i, j, k)
                                                     ? initial_value(desc, m)
                                                     : 0.0;
            }
        }

        boundary.ApplyAllPhysicalBoundaries();
        boundary.SyncGroup("U");
        boundary.SyncGroup("V");
        boundary.SyncGroup("rho");
        boundary.SyncGroup("Bface");
        boundary.SyncGroup("Eedge");
    }
}
