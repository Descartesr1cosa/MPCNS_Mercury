#pragma once

#include "3_field/FieldArray.h"
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

    Vector dJcell_w; // vec_length = 36
    // [0..11]   -> Jx 对 12 条边的权重
    // [12..23]  -> Jy 对 12 条边的权重
    // [24..35]  -> Jz 对 12 条边的权重

    Vector dBcell_w; // vec_length = 18
    // [0..5]   -> Bx 对 6 个面的权重
    // [6..11]  -> By 对 6 个面的权重
    // [12..17] -> Bz 对 6 个面的权重

    Vector P_xi;   // face projector on xi-faces, vec_length = 6
    Vector P_eta;  // face projector on eta-faces, vec_length = 6
    Vector P_zeta; // face projector on zeta-faces, vec_length = 6

    Vector dr_xi;   // edge direction vector on xi-edges, vec_length = 3
    Vector dr_eta;  // edge direction vector on eta-edges, vec_length = 3
    Vector dr_zeta; // edge direction vector on zeta-edges, vec_length = 3
};
