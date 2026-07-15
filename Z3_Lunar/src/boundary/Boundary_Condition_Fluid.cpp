#include "1_Boundary.h"

#include <algorithm>
#include <array>
#include <cmath>

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

void LunarBoundary::BC_UH_Absorbing_(FieldBlock &U, Field *fld,
                                     const BOUND::PhysicalRegion &r, int ngh)
{
    if (!U.is_allocated() || !fld || ngh <= 0)
        return;

    const Box3 &inner = r.inner_slab;
    const int axis = std::abs(r.direction);
    const int sign_to_ghost = (r.direction > 0) ? 1 : -1;
    const Int3 outward = (axis == 1)   ? Int3{sign_to_ghost, 0, 0}
                         : (axis == 2) ? Int3{0, sign_to_ghost, 0}
                                       : Int3{0, 0, sign_to_ghost};
    const int face_shift = (sign_to_ghost > 0) ? 1 : 0;

    FieldBlock *area = nullptr;
    if (axis == 1)
        area = &fld->field("JDxi", r.this_block);
    else if (axis == 2)
        area = &fld->field("JDet", r.this_block);
    else if (axis == 3)
        area = &fld->field("JDze", r.this_block);
    else
        return;

    if (!area->is_allocated())
        return;

    constexpr double rho_floor = 1.0e-12;
    constexpr double p_floor = 1.0e-12;
    const double gm1 = std::max(bc_state_.gamma - 1.0, 1.0e-12);

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
            {
                const int fi = i + outward.i * face_shift;
                const int fj = j + outward.j * face_shift;
                const int fk = k + outward.k * face_shift;

                std::array<double, 3> n = {
                    (*area)(fi, fj, fk, 0),
                    (*area)(fi, fj, fk, 1),
                    (*area)(fi, fj, fk, 2)};
                const double area_norm = std::sqrt(
                    n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                if (!(area_norm > 1.0e-30) || !std::isfinite(area_norm))
                    continue;

                // n points from the lunar surface into the plasma domain.
                const double into_fluid = (sign_to_ghost > 0) ? -1.0 : 1.0;
                for (double &component : n)
                    component *= into_fluid / area_norm;

                const double rho = std::max(U(i, j, k, 0), rho_floor);
                const std::array<double, 3> velocity = {
                    U(i, j, k, 1) / rho,
                    U(i, j, k, 2) / rho,
                    U(i, j, k, 3) / rho};
                const double speed2 = velocity[0] * velocity[0] +
                                      velocity[1] * velocity[1] +
                                      velocity[2] * velocity[2];
                const double pressure = std::max(
                    gm1 * (U(i, j, k, 4) - 0.5 * rho * speed2), p_floor);
                const double vn = velocity[0] * n[0] +
                                  velocity[1] * n[1] +
                                  velocity[2] * n[2];

                std::array<double, 3> ghost_velocity = velocity;
                if (vn > 0.0)
                {
                    // The active state points from the surface into the
                    // domain. Reflect only that normal component in the ghost
                    // state, so the reconstructed face-normal velocity is
                    // zero while tangential velocity is retained.
                    ghost_velocity[0] -= 2.0 * vn * n[0];
                    ghost_velocity[1] -= 2.0 * vn * n[1];
                    ghost_velocity[2] -= 2.0 * vn * n[2];
                }
                // vn <= 0 points into the surface: copy the state to provide
                // a zero-gradient, freely absorbing outflow.

                const double ghost_speed2 =
                    ghost_velocity[0] * ghost_velocity[0] +
                    ghost_velocity[1] * ghost_velocity[1] +
                    ghost_velocity[2] * ghost_velocity[2];
                const double ghost_energy = pressure / gm1 +
                                            0.5 * rho * ghost_speed2;

                for (int layer = 1; layer <= ngh; ++layer)
                {
                    const int ig = i + layer * outward.i;
                    const int jg = j + layer * outward.j;
                    const int kg = k + layer * outward.k;
                    U(ig, jg, kg, 0) = rho;
                    U(ig, jg, kg, 1) = rho * ghost_velocity[0];
                    U(ig, jg, kg, 2) = rho * ghost_velocity[1];
                    U(ig, jg, kg, 3) = rho * ghost_velocity[2];
                    U(ig, jg, kg, 4) = ghost_energy;
                }
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
