#include "Z0_Initial.h"

#include "Z0_Boundary.h"
#include "0_basic/LayoutTraits.h"
#include "3_field/Field.h"

#include <cstdlib>
#include <string>

namespace
{
    bool owner_only_enabled()
    {
        const char *value = std::getenv("Z0_OWNER_ONLY");
        if (!value)
            return false;
        const std::string s(value);
        return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
    }

    bool in_box(const Box3 &b, int i, int j, int k)
    {
        return i >= b.lo.i && i < b.hi.i &&
               j >= b.lo.j && j < b.hi.j &&
               k >= b.lo.k && k < b.hi.k;
    }

    double initial_value(const FieldDescriptor &desc, int comp)
    {
        (void)desc;
        (void)comp;
        return 0.0;
    }
}

namespace Z0
{
    void InitializeFields(Field &field, Param &param, Z0_Boundary &boundary)
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
        if (!owner_only_enabled())
            boundary.SyncAllRegistered();
    }
}
