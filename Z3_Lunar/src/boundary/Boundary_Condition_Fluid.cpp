#include "1_Boundary.h"

void LunarBoundary::BC_UH_Farfield_H_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                U(i, j, k, 0) = bc_state_.qinf[0];
                U(i, j, k, 1) = bc_state_.qinf[1];
                U(i, j, k, 2) = bc_state_.qinf[2];
                U(i, j, k, 3) = bc_state_.qinf[3];
                U(i, j, k, 4) = bc_state_.qinf[4];
            }
}


void LunarBoundary::BC_Pole_Cell_(FieldBlock &U, Field *fld,
                                    const BOUND::PhysicalRegion &r, int ngh)
{
    if (!U.is_allocated())
        return;
    const Box3 &inner = r.inner_slab;

    // 这里把 ncomp() 换成你 FieldBlock 实际的分量数接口
    const int ncomp = U.descriptor().ncomp;

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
            for (int n = 0; n < ncomp; ++n)
            {
                double temp_U = 0.0;
                double num_d = 0.0;

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    temp_U += U(i, j, k, n);
                    num_d += 1.0;
                }

                if (num_d > 0.0)
                    temp_U /= num_d;
                else
                    temp_U = 0.0;

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                    U(i, j, k, n) = temp_U;
            }

    bound_.DefaultPhysicalCopy(U, fld, r, ngh);
}
