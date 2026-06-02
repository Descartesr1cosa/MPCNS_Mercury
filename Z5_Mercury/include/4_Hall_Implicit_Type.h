#pragma once

#include "00_Mercury_Const.h"
#include "3_field/Field_Array.h"
#include "0_basic/TYPES.h"

struct HallFaceScratchBlock_ // For Rusanov Scheme
{
    Int3 clo{}, chi{};
    int ni = 0, nj = 0, nk = 0;

    Vector Ehc;  // vec_length = 3, 分别存 x/y/z
    Scalar beta; // 标量 beta_hall

    // ----- Whistler-only P0 需要的 frozen / linear scratch -----
    Vector Bflat;      // frozen Bcell
    Scalar alpha_flat; // frozen alpha on cell
    Vector dEhc;       // linearized whistler-only cell Hall EMF

    // Geometry caches used by Hall live in Field:
    //   Bcell_from_Bface_w, Jcell_from_Jedge_w, Pface_xi/eta/zeta, dr_xi/eta/zeta.
};
