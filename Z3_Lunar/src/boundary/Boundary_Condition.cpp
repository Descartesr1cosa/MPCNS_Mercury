#include "1_Boundary.h"

void LunarBoundary::BC_Pole_Eedge_Average_Axis(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            double temp_E = 0.0;
            double num_d = 0.0;

            for (int k = inner.lo.k; k < inner.hi.k - 1; ++k) // avoid recounting: k=0 <=> k=inner.hi.k-1
            {
                temp_E += U(i, j, k, 0);
                num_d += 1.0;
            }
            temp_E /= num_d;

            for (int k = inner.lo.k; k < inner.hi.k; ++k)
                U(i, j, k, 0) = temp_E;
        }

    bound_.DefaultPhysicalCopy(U, fld, r, ngh);
}

void LunarBoundary::BC_Pole_Eedge_Axis(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    // ---------------------------------------------------------------------
    // Purpose:
    //   Pole 上的 axis-edge EMF 不再由 Eface->Eedge 的平均值修正。
    //   直接假设 Pole 第一层 u、B_cell 是唯一 Cartesian 常矢量：
    //
    //       E_cart = - u x B_cell
    //       E_axis = E_cart · dr_axis
    //
    //   当前函数应只用于：
    //       xi-Pole  : E_eta 是 axis edge
    //       eta-Pole : E_xi  是 axis edge
    //
    //   rotate edge(E_zeta) 继续由 BC_Pole_Eedge_Zero_Rotate 处理。
    //   norm edge 之后另写。
    // ---------------------------------------------------------------------

    if (!fld)
    {
        bound_.DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    const int ib = r.this_block;
    const int pole_norm = std::abs(r.direction); // 1: xi-pole, 2: eta-pole
    const std::string fname = U.descriptor().name;

    int edge_axis = 0;
    std::string dr_name;

    if (fname == "E_xi")
    {
        edge_axis = 1;
        dr_name = "dr_xi";
    }
    else if (fname == "E_eta")
    {
        edge_axis = 2;
        dr_name = "dr_eta";
    }
    else
    {
        bound_.DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    // Only axis edge:
    //   eta-pole -> E_xi
    //   xi-pole  -> E_eta
    const bool is_axis_edge =
        (pole_norm == 1 && edge_axis == 2) ||
        (pole_norm == 2 && edge_axis == 1);

    if (!is_axis_edge)
    {
        bound_.DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    auto get_block = [&](const std::string &name) -> FieldBlock *
    {
        if (!fld->has_field(name))
            return nullptr;

        FieldBlock &F = fld->field(name, ib);

        if (!F.is_allocated())
            return nullptr;

        return &F;
    };

    FieldBlock *Uplus_ptr = get_block("U_plus");
    FieldBlock *Bcell_ptr = get_block("B_cell");
    FieldBlock *dr_ptr = get_block(dr_name);

    if (!Uplus_ptr || !Bcell_ptr || !dr_ptr)
    {
        bound_.DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    FieldBlock &Uplus = *Uplus_ptr;
    FieldBlock &Bcell = *Bcell_ptr;
    FieldBlock &dr = *dr_ptr;

    const Box3 &inner = r.inner_slab;

    auto dot3 = [](const double a[3], const double b[3]) -> double
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    int iplus = 0, jplus = 0;
    if (r.direction > 0)
    {
        iplus = (edge_axis == 1) ? 0 : -1;
        jplus = (edge_axis == 2) ? 0 : -1;
    }

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
    {
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
            {
                int ic = i + iplus;
                int jc = j + jplus;
                int kc = k;

                const double ux = Uplus(ic, jc, kc, 0);
                const double uy = Uplus(ic, jc, kc, 1);
                const double uz = Uplus(ic, jc, kc, 2);

                const double Bx = Bcell(ic, jc, kc, 0);
                const double By = Bcell(ic, jc, kc, 1);
                const double Bz = Bcell(ic, jc, kc, 2);

                // E = - u_plus x B_cell
                const double Ecart[3] = {
                    -(uy * Bz - uz * By),
                    -(uz * Bx - ux * Bz),
                    -(ux * By - uy * Bx)};

                const double dl[3] = {
                    dr(i, j, k, 0),
                    dr(i, j, k, 1),
                    dr(i, j, k, 2)};

                U(i, j, k, 0) = dot3(Ecart, dl);
            }
        }
    }

    bound_.DefaultPhysicalCopy(U, fld, r, ngh);
}

void LunarBoundary::BC_Pole_Eedge_Zero_Rotate(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
                U(i, j, k, 0) = 0.0;
        }

    bound_.DefaultPhysicalCopy(U, fld, r, ngh);
}

void LunarBoundary::BC_Pole_Eedge_RegulateK_Norm(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;

    const int ax = std::abs(r.direction);
    const int sgn = (r.direction > 0) ? +1 : -1;

    // 当前 handler 只应被用于：
    //   E_xi  with ax == 1
    //   E_eta with ax == 2
    // 若以后扩展 E_zeta，也可 ax == 3。
    Int3 nrm{0, 0, 0};
    if (ax == 1)
        nrm = Int3{sgn, 0, 0};
    else if (ax == 2)
        nrm = Int3{0, sgn, 0};
    else
        nrm = Int3{0, 0, sgn};

    const Int3 ulo = U.get_lo();
    const Int3 uhi = U.get_hi();

    auto in_range = [&](int i, int j, int k) -> bool
    {
        return (i >= ulo.i && i < uhi.i &&
                j >= ulo.j && j < uhi.j &&
                k >= ulo.k && k < uhi.k);
    };

    // ------------------------------------------------------------
    // 取当前 edge family 对应的物理边向量 dr。
    //
    // E_xi   -> dr_xi
    // E_eta  -> dr_eta
    // E_zeta -> dr_zeta
    // ------------------------------------------------------------

    const int ib = r.this_block;

    FieldBlock *dr_ptr = nullptr;

    if (ax == 1)
        dr_ptr = &((*fld).field("dr_xi")[ib]);
    else if (ax == 2)
        dr_ptr = &((*fld).field("dr_eta")[ib]);
    else
        dr_ptr = &((*fld).field("dr_zeta")[ib]);

    FieldBlock &dr = *dr_ptr;

    auto get_dr = [&](int i, int j, int k) -> std::array<double, 3>
    {
        return {dr(i, j, k, 0), dr(i, j, k, 1), dr(i, j, k, 2)};
    };

    auto dot3 = [](const std::array<double, 3> &a,
                   const std::array<double, 3> &b) -> double
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    auto norm3 = [&](const std::array<double, 3> &a) -> double
    {
        return std::sqrt(dot3(a, a));
    };

    auto solve3x3 = [](double A[3][3], double b[3], double x[3]) -> bool
    {
        double M[3][4] = {
            {A[0][0], A[0][1], A[0][2], b[0]},
            {A[1][0], A[1][1], A[1][2], b[1]},
            {A[2][0], A[2][1], A[2][2], b[2]}};

        for (int c = 0; c < 3; ++c)
        {
            int piv = c;
            double amax = std::abs(M[c][c]);

            for (int r0 = c + 1; r0 < 3; ++r0)
            {
                const double av = std::abs(M[r0][c]);
                if (av > amax)
                {
                    amax = av;
                    piv = r0;
                }
            }

            if (amax < 1.0e-300)
                return false;

            if (piv != c)
            {
                for (int q = c; q < 4; ++q)
                    std::swap(M[c][q], M[piv][q]);
            }

            const double inv_piv = 1.0 / M[c][c];

            for (int q = c; q < 4; ++q)
                M[c][q] *= inv_piv;

            for (int r0 = 0; r0 < 3; ++r0)
            {
                if (r0 == c)
                    continue;

                const double fac = M[r0][c];

                for (int q = c; q < 4; ++q)
                    M[r0][q] -= fac * M[c][q];
            }
        }

        x[0] = M[0][3];
        x[1] = M[1][3];
        x[2] = M[2][3];

        return true;
    };

    constexpr double eps_len = 1.0e-14;
    constexpr double eps_reg = 1.0e-12;

    // ============================================================
    // 1. 对第一层物理边界做：
    //
    //      edge 1-form values -> unique Cartesian E vector -> edge 1-form values
    //
    // 对每个固定的 (i,j)，收集所有 k 上的当前 edge：
    //
    //      e_k = U(i,j,k,0)
    //      dr_k = dr_axis(i,j,k,:)
    //
    // 解：
    //
    //      e_k / |dr_k| ≈ E · (dr_k / |dr_k|)
    //
    // 然后：
    //
    //      U(i,j,k,0) = E · dr_k
    //
    // 这样保证同一个 Pole 点的一圈 edge 由唯一 E_vec 生成。
    // ============================================================

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
    {
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            double A[3][3] = {
                {0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0}};

            double b[3] = {0.0, 0.0, 0.0};

            int cnt = 0;

            for (int k = inner.lo.k; k < inner.hi.k - 1; ++k)
            {
                if (!in_range(i, j, k))
                    continue;

                const auto dl = get_dr(i, j, k);
                const double L = norm3(dl);

                // 对 collapsed edge，不参与矢量重构。
                // 后面会投影成 0。
                if (L < eps_len)
                    continue;

                const double invL = 1.0 / L;

                const double tau[3] = {
                    dl[0] * invL,
                    dl[1] * invL,
                    dl[2] * invL};

                const double y = U(i, j, k, 0) * invL;

                for (int a = 0; a < 3; ++a)
                {
                    b[a] += y * tau[a];

                    for (int c = 0; c < 3; ++c)
                        A[a][c] += tau[a] * tau[c];
                }

                ++cnt;
            }

            if (cnt <= 0)
            {
                // 所有 edge 都 collapsed，则直接置零这一圈。
                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    if (in_range(i, j, k))
                        U(i, j, k, 0) = 0.0;
                }
                continue;
            }

            const double trA = A[0][0] + A[1][1] + A[2][2];
            const double reg = eps_reg * std::max(1.0, trA);

            A[0][0] += reg;
            A[1][1] += reg;
            A[2][2] += reg;

            double Evec[3] = {0.0, 0.0, 0.0};

            const bool ok = solve3x3(A, b, Evec);

            if (!ok)
            {
                // 极端退化时，退回到代数平均，但仍然按 1-form 处理。
                // 这里不建议直接保留原值，因为原值可能包含周向噪声。
                double avg = 0.0;
                int navg = 0;

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    if (!in_range(i, j, k))
                        continue;

                    avg += U(i, j, k, 0);
                    ++navg;
                }

                if (navg > 0)
                    avg /= static_cast<double>(navg);

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    if (in_range(i, j, k))
                        U(i, j, k, 0) = avg;
                }

                continue;
            }

            // 用唯一 Evec 重新投影回每个 edge 的 1-form。
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
            {
                if (!in_range(i, j, k))
                    continue;

                const auto dl = get_dr(i, j, k);
                const double L = norm3(dl);

                if (L < eps_len)
                {
                    U(i, j, k, 0) = 0.0;
                }
                else
                {
                    U(i, j, k, 0) =
                        Evec[0] * dl[0] +
                        Evec[1] * dl[1] +
                        Evec[2] * dl[2];
                }
            }
        }
    }

    // ============================================================
    // 2. ghost 从 regulation 后的第一层 copy。
    //
    // 注意：
    // 这里 copy 的是已经由 Evec 投影回来的 1-form。
    // ============================================================

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
    {
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
            {
                if (!in_range(i, j, k))
                    continue;

                const double val = U(i, j, k, 0);

                for (int g = 1; g <= ngh; ++g)
                {
                    const int ig = i + g * nrm.i;
                    const int jg = j + g * nrm.j;
                    const int kg = k + g * nrm.k;

                    if (!in_range(ig, jg, kg))
                        continue;

                    U(ig, jg, kg, 0) = val;
                }
            }
        }
    }
}

void LunarBoundary::BC_Pole_Jedge_RegulateK_Norm(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;
    const int ax = std::abs(r.direction);

    auto zero_ghost = [&]()
    {
        const Box3 g = BoundaryCore::MakeGhostSlabFromInner(inner, r.direction, ngh);
        for (int i = g.lo.i; i < g.hi.i; ++i)
            for (int j = g.lo.j; j < g.hi.j; ++j)
                for (int k = g.lo.k; k < g.hi.k; ++k)
                    U(i, j, k, 0) = 0.0;
    };

    if (!fld || (ax != 1 && ax != 2))
    {
        zero_ghost();
        return;
    }

    const std::string fname = U.descriptor().name;
    const int ib = r.this_block;

    auto get_block = [&](const std::string &name) -> FieldBlock *
    {
        if (!fld->has_field(name))
            return nullptr;

        FieldBlock &F = fld->field(name, ib);
        if (!F.is_allocated())
            return nullptr;

        return &F;
    };

    constexpr double eps_area = 1.0e-300;

    // Pole norm edge:
    // Use beta*Bface to convert the magnetic 2-form on a primal face to the
    // line integral on the corresponding dual edge, then use the uncapped
    // dl/Astar ratio to get the edge current 1-form. The stored alpha_* is
    // intentionally not used here because it is capped to zero when the dual
    // face touches the Pole.
    if (ax == 1 && fname == "J_xi")
    {
        FieldBlock *Bze = get_block("B_zeta");
        FieldBlock *beta_ze = get_block("beta_zeta");
        FieldBlock *dl_xi = get_block("dl_xi");
        FieldBlock *Ast_xi = get_block("Astar_xi");

        if (Bze && beta_ze && dl_xi && Ast_xi)
        {
            for (int i = inner.lo.i; i < inner.hi.i; ++i)
                for (int j = inner.lo.j; j < inner.hi.j; ++j)
                    for (int k = inner.lo.k; k < inner.hi.k; ++k)
                    {
                        const double circ =
                            (*beta_ze)(i, j, k, 0) * (*Bze)(i, j, k, 0) -
                            (*beta_ze)(i, j - 1, k, 0) * (*Bze)(i, j - 1, k, 0);

                        const double Astar = (*Ast_xi)(i, j, k, 0);
                        const double dl = (*dl_xi)(i, j, k, 0);

                        U(i, j, k, 0) = (std::abs(Astar) <= eps_area) ? 0.0 : (dl / Astar) * circ;
                    }
        }
    }
    else if (ax == 2 && fname == "J_eta")
    {
        FieldBlock *Bze = get_block("B_zeta");
        FieldBlock *beta_ze = get_block("beta_zeta");
        FieldBlock *dl_et = get_block("dl_eta");
        FieldBlock *Ast_et = get_block("Astar_eta");

        if (Bze && beta_ze && dl_et && Ast_et)
        {
            for (int i = inner.lo.i; i < inner.hi.i; ++i)
                for (int j = inner.lo.j; j < inner.hi.j; ++j)
                    for (int k = inner.lo.k; k < inner.hi.k; ++k)
                    {
                        const double circ =
                            (*beta_ze)(i - 1, j, k, 0) * (*Bze)(i - 1, j, k, 0) -
                            (*beta_ze)(i, j, k, 0) * (*Bze)(i, j, k, 0);

                        const double Astar = (*Ast_et)(i, j, k, 0);
                        const double dl = (*dl_et)(i, j, k, 0);

                        U(i, j, k, 0) = (std::abs(Astar) <= eps_area) ? 0.0 : (dl / Astar) * circ;
                    }
        }
    }

    zero_ghost();
}

void LunarBoundary::BC_Farfield_Eedge_set_zerocurl(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    double u_inf[3] = {bc_state_.q_pv_inf[0], bc_state_.q_pv_inf[1], bc_state_.q_pv_inf[2]};
    double B_inf[3] = {bc_state_.B_imf[0], bc_state_.B_imf[1], bc_state_.B_imf[2]};

    double E_inf[3]; // = - u_inf \times B_inf
    E_inf[0] = B_inf[1] * u_inf[2] - B_inf[2] * u_inf[1];
    E_inf[1] = B_inf[2] * u_inf[0] - B_inf[0] * u_inf[2];
    E_inf[2] = B_inf[0] * u_inf[1] - B_inf[1] * u_inf[0];

    int dijk[3] = {0, 0, 0};

    if (U.descriptor().name == "E_xi")
    {
        dijk[0] = 1;
    }
    else if (U.descriptor().name == "E_eta")
    {
        dijk[1] = 1;
    }
    else if (U.descriptor().name == "E_zeta")
    {
        dijk[2] = 1;
    }
    else
    {
        bound_.DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    int iblock = r.this_block;
    double3D &x = fld->grd->grids(iblock).x;
    double3D &y = fld->grd->grids(iblock).y;
    double3D &z = fld->grd->grids(iblock).z;

    const Box3 &inner = r.inner_slab;

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
            {
                const int ir = i + dijk[0];
                const int jr = j + dijk[1];
                const int kr = k + dijk[2];
                U(i, j, k, 0) = E_inf[0] * (x(ir, jr, kr) - x(i, j, k)) +
                                E_inf[1] * (y(ir, jr, kr) - y(i, j, k)) +
                                E_inf[2] * (z(ir, jr, kr) - z(i, j, k));
            }
    bound_.DefaultPhysicalCopy(U, fld, r, ngh);
    return;
}

void LunarBoundary::BC_Farfield_Bface(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write

    if (U.descriptor().name == "B_cell")
    {
        for (int i = g.lo.i; i < g.hi.i; ++i)
            for (int j = g.lo.j; j < g.hi.j; ++j)
                for (int k = g.lo.k; k < g.hi.k; ++k)
                {
                    U(i, j, k, 0) = bc_state_.B_imf[0];
                    U(i, j, k, 1) = bc_state_.B_imf[1];
                    U(i, j, k, 2) = bc_state_.B_imf[2];
                }
    }
    else if (U.descriptor().name == "Bind_cell")
    {
        for (int i = g.lo.i; i < g.hi.i; ++i)
            for (int j = g.lo.j; j < g.hi.j; ++j)
                for (int k = g.lo.k; k < g.hi.k; ++k)
                {
                    U(i, j, k, 0) = 0.0;
                    U(i, j, k, 1) = 0.0;
                    U(i, j, k, 2) = 0.0;
                }
    }
    else
    {
        std::cerr << "BC_Farfield_Bface called for unexpected field: "
                  << U.descriptor().name << std::endl;
        return;
    }
}
