#include "1_Boundary.h"

void MercuryBoundary::BC_Pole_Bcell_Collapse_(FieldBlock &U, Field *fld,
                                              const BOUND::PhysicalRegion &r,
                                              int ngh)
{
    const int dir = std::abs(r.direction); // 1: xi, 2: eta
    const int sgn = (r.direction > 0) ? +1 : -1;

    if (dir != 1 && dir != 2)
    {
        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    const std::string field_name = U.descriptor().name;

    const bool is_total_Bcell = (field_name == "B_cell");
    const bool is_bind_Bcell = (field_name == "Bind_cell");

    if (!is_total_Bcell && !is_bind_Bcell)
    {
        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    const int ib = r.this_block;
    const Box3 &inner = r.inner_slab;

    auto get_required = [&](const std::string &name) -> FieldBlock *
    {
        if (!fld->has_field(name))
            return nullptr;

        FieldBlock &F = fld->field(name, ib);

        if (!F.is_allocated())
            return nullptr;

        return &F;
    };

    auto get_optional = [&](const std::string &name) -> FieldBlock *
    {
        if (!fld->has_field(name))
            return nullptr;

        FieldBlock &F = fld->field(name, ib);

        if (!F.is_allocated())
            return nullptr;

        return &F;
    };

    FieldBlock *Bxi = get_required("B_xi");
    FieldBlock *Bet = get_required("B_eta");

    FieldBlock *Axi = get_required("JDxi");
    FieldBlock *Aet = get_required("JDet");

    FieldBlock *Baddxi = get_optional("Badd_xi");
    FieldBlock *Baddet = get_optional("Badd_eta");

    if (!Bxi || !Bet || !Axi || !Aet)
    {
        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    auto face_value = [&](FieldBlock *F, FieldBlock *Fadd,
                          int i, int j, int k) -> double
    {
        double v = (*F)(i, j, k, 0);

        // B_cell    : total B = induced Bface + background Badd
        // Bind_cell : induced B only
        if (is_total_Bcell && Fadd)
            v += (*Fadd)(i, j, k, 0);

        return v;
    };

    auto push_eq = [](double Sx, double Sy, double Sz, double phi,
                      double M[3][3], double b[3])
    {
        const double sn2 = Sx * Sx + Sy * Sy + Sz * Sz;

        if (sn2 <= 1.0e-300)
            return;

        const double S[3] = {Sx, Sy, Sz};

        for (int a = 0; a < 3; ++a)
        {
            b[a] += phi * S[a];

            for (int c = 0; c < 3; ++c)
                M[a][c] += S[a] * S[c];
        }
    };

    auto solve3 = [](double A[3][3], double b[3], double x[3]) -> bool
    {
        const double tr =
            std::abs(A[0][0]) +
            std::abs(A[1][1]) +
            std::abs(A[2][2]);

        const double reg = 1.0e-12 * std::max(1.0, tr);

        A[0][0] += reg;
        A[1][1] += reg;
        A[2][2] += reg;

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

            const double inv = 1.0 / M[c][c];

            for (int q = c; q < 4; ++q)
                M[c][q] *= inv;

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

    // =====================================================================
    // xi-Pole:
    //
    //   norm = xi / i
    //   axis = eta / j
    //   zeta = k
    //
    // Collapsed cell ring:
    //   (i0, j, k = all)
    //
    // Use only:
    //   1. outer norm face: B_xi(i_norm_face, j, k)
    //   2. lower axis face: B_eta(i0, j,   k)
    //   3. upper axis face: B_eta(i0, j+1, k)
    //
    // Do NOT use B_zeta.
    // =====================================================================
    if (dir == 1)
    {
        for (int i0 = inner.lo.i; i0 < inner.hi.i; ++i0)
        {
            // inner range is boundary range for Cell
            const int i_norm_face = (sgn < 0) ? (i0 + 1) : i0;

            // low-side Pole:
            //   collapsed cell lower xi face is degenerate;
            //   outer xi face is upper face, outward +xi.
            //
            // high-side Pole:
            //   collapsed cell upper xi face is degenerate;
            //   outer xi face is lower face, outward -xi.
            const double s_norm = (sgn < 0) ? +1.0 : -1.0;

            for (int j = inner.lo.j; j < inner.hi.j; ++j)
            {
                double M[3][3] = {
                    {0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0}};

                double b[3] = {0.0, 0.0, 0.0};

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    // ----------------------------------------------------
                    // 1. Outer norm face: B_xi(i_norm_face, j, k)
                    // ----------------------------------------------------
                    {
                        const double Sx = s_norm * (*Axi)(i_norm_face, j, k, 0);
                        const double Sy = s_norm * (*Axi)(i_norm_face, j, k, 1);
                        const double Sz = s_norm * (*Axi)(i_norm_face, j, k, 2);

                        const double phi = s_norm * face_value(Bxi, Baddxi,
                                                               i_norm_face, j, k);

                        push_eq(Sx, Sy, Sz, phi, M, b);
                    }

                    // ----------------------------------------------------
                    // 2. Lower axis face: B_eta(i0, j, k)
                    //    outward direction = -eta
                    // ----------------------------------------------------
                    {
                        const double Sx = -(*Aet)(i0, j, k, 0);
                        const double Sy = -(*Aet)(i0, j, k, 1);
                        const double Sz = -(*Aet)(i0, j, k, 2);

                        const double phi = -face_value(Bet, Baddet,
                                                       i0, j, k);

                        push_eq(Sx, Sy, Sz, phi, M, b);
                    }

                    // ----------------------------------------------------
                    // 3. Upper axis face: B_eta(i0, j + 1, k)
                    //    outward direction = +eta
                    // ----------------------------------------------------
                    {
                        const double Sx = (*Aet)(i0, j + 1, k, 0);
                        const double Sy = (*Aet)(i0, j + 1, k, 1);
                        const double Sz = (*Aet)(i0, j + 1, k, 2);

                        const double phi = face_value(Bet, Baddet,
                                                      i0, j + 1, k);

                        push_eq(Sx, Sy, Sz, phi, M, b);
                    }
                }

                double Bp[3] = {0.0, 0.0, 0.0};

                if (!solve3(M, b, Bp))
                    continue;

                // Write the same collapsed B_pole to all k.
                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    U(i0, j, k, 0) = Bp[0];
                    U(i0, j, k, 1) = Bp[1];
                    U(i0, j, k, 2) = Bp[2];
                }
            }
        }
    }

    // =====================================================================
    // eta-Pole:
    //
    //   norm = eta / j
    //   axis = xi / i
    //   zeta = k
    //
    // Collapsed cell ring:
    //   (i, j0, k = all)
    //
    // Use only:
    //   1. outer norm face: B_eta(i, j_norm_face, k)
    //   2. lower axis face: B_xi(i,   j0, k)
    //   3. upper axis face: B_xi(i+1, j0, k)
    //
    // Do NOT use B_zeta.
    // =====================================================================
    if (dir == 2)
    {
        for (int j0 = inner.lo.j; j0 < inner.hi.j; ++j0)
        {
            // inner range is boundary range for Cell
            const int j_norm_face = (sgn < 0) ? (j0 + 1) : j0;

            // low-side Pole:
            //   outer eta face is upper face, outward +eta.
            //
            // high-side Pole:
            //   outer eta face is lower face, outward -eta.
            const double s_norm = (sgn < 0) ? +1.0 : -1.0;

            for (int i = inner.lo.i; i < inner.hi.i; ++i)
            {
                double M[3][3] = {
                    {0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0}};

                double b[3] = {0.0, 0.0, 0.0};

                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    // ----------------------------------------------------
                    // 1. Outer norm face: B_eta(i, j_norm_face, k)
                    // ----------------------------------------------------
                    {
                        const double Sx = s_norm * (*Aet)(i, j_norm_face, k, 0);
                        const double Sy = s_norm * (*Aet)(i, j_norm_face, k, 1);
                        const double Sz = s_norm * (*Aet)(i, j_norm_face, k, 2);

                        const double phi = s_norm * face_value(Bet, Baddet,
                                                               i, j_norm_face, k);

                        push_eq(Sx, Sy, Sz, phi, M, b);
                    }

                    // ----------------------------------------------------
                    // 2. Lower axis face: B_xi(i, j0, k)
                    //    outward direction = -xi
                    // ----------------------------------------------------
                    {
                        const double Sx = -(*Axi)(i, j0, k, 0);
                        const double Sy = -(*Axi)(i, j0, k, 1);
                        const double Sz = -(*Axi)(i, j0, k, 2);

                        const double phi = -face_value(Bxi, Baddxi,
                                                       i, j0, k);

                        push_eq(Sx, Sy, Sz, phi, M, b);
                    }

                    // ----------------------------------------------------
                    // 3. Upper axis face: B_xi(i + 1, j0, k)
                    //    outward direction = +xi
                    // ----------------------------------------------------
                    {
                        const double Sx = (*Axi)(i + 1, j0, k, 0);
                        const double Sy = (*Axi)(i + 1, j0, k, 1);
                        const double Sz = (*Axi)(i + 1, j0, k, 2);

                        const double phi = face_value(Bxi, Baddxi,
                                                      i + 1, j0, k);

                        push_eq(Sx, Sy, Sz, phi, M, b);
                    }
                }

                double Bp[3] = {0.0, 0.0, 0.0};

                if (!solve3(M, b, Bp))
                    continue;

                // Write the same collapsed B_pole to all k.
                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                {
                    U(i, j0, k, 0) = Bp[0];
                    U(i, j0, k, 1) = Bp[1];
                    U(i, j0, k, 2) = Bp[2];
                }
            }
        }
    }

    // ghost 从 collapse 后的 inner_slab copy。
    BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
}

void MercuryBoundary::BC_Solid_Surface_Jcell(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &inner = r.inner_slab;
    const int ncomp = U.descriptor().ncomp;
    const int ax = std::abs(r.direction);
    const int sgn = (r.direction > 0) ? +1 : -1;
    const int inward = -sgn;
    const Int3 lo = U.inner_lo();
    const Int3 hi = U.inner_hi();

    auto inside = [&](int i, int j, int k) -> bool
    {
        return i >= lo.i && i < hi.i &&
               j >= lo.j && j < hi.j &&
               k >= lo.k && k < hi.k;
    };

    for (int i = inner.lo.i; i < inner.hi.i; ++i)
        for (int j = inner.lo.j; j < inner.hi.j; ++j)
        {
            for (int k = inner.lo.k; k < inner.hi.k; ++k)
            {
                for (int n = 0; n < ncomp; n++)
                    U(i, j, k, n) = 0.0;
                int is = i;
                int js = j;
                int ks = k;

                if (ax == 1)
                    is += inward;
                else if (ax == 2)
                    js += inward;
                else
                    ks += inward;

                for (int n = 0; n < ncomp; n++)
                    U(is, js, ks, n) = 0.0;
                // if (!inside(is, js, ks))
                //     continue;

                // for (int n = 0; n < ncomp; n++)
                //     U(i, j, k, n) = U(is, js, ks, n);
            }
        }
    BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
}

void MercuryBoundary::BC_Pole_Bface_Collapse_(FieldBlock &U, Field *fld,
                                              const BOUND::PhysicalRegion &r,
                                              int ngh)
{
    const int dir = std::abs(r.direction); // 1: xi, 2: eta, 3: zeta

    // 这里只处理旋转轴为 zeta/k 的 Pole：norm = xi 或 eta。
    if (dir != 1 && dir != 2)
    {
        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    const auto loc = U.descriptor().location;

    const bool is_Bxi = (loc == StaggerLocation::FaceXi);
    const bool is_Beta = (loc == StaggerLocation::FaceEt);
    const bool is_Bze = (loc == StaggerLocation::FaceZe);

    const Box3 &inner = r.inner_slab;

    // ------------------------------------------------------------
    // 1. norm face | norm=0 -> 0
    //
    // Pole normal xi  : B_xi  is norm face.
    // Pole normal eta : B_eta is norm face.
    // ------------------------------------------------------------
    const bool is_norm_face =
        (dir == 1 && is_Bxi) ||
        (dir == 2 && is_Beta);

    if (is_norm_face)
    {
        for (int i = inner.lo.i; i < inner.hi.i; ++i)
            for (int j = inner.lo.j; j < inner.hi.j; ++j)
                for (int k = inner.lo.k; k < inner.hi.k; ++k)
                    U(i, j, k, 0) = 0.0;

        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    // ------------------------------------------------------------
    // 2. axis face | norm=0 -> k-ring average
    //
    // Pole normal xi  : axis face = B_eta
    // Pole normal eta : axis face = B_xi
    // ------------------------------------------------------------
    const bool is_axis_face =
        (dir == 1 && is_Beta) ||
        (dir == 2 && is_Bxi);

    // if (is_axis_face)
    // {
    //     const int ib = r.this_block;

    //     FieldBlock *Aaxis = nullptr;

    //     if (fld)
    //     {
    //         if (dir == 1 && is_Beta)
    //         {
    //             // xi-Pole: axis face = B_eta, area vector = JDet
    //             if (fld->has_field("JDet"))
    //                 Aaxis = &fld->field("JDet", ib);
    //         }
    //         else if (dir == 2 && is_Bxi)
    //         {
    //             // eta-Pole: axis face = B_xi, area vector = JDxi
    //             if (fld->has_field("JDxi"))
    //                 Aaxis = &fld->field("JDxi", ib);
    //         }
    //     }

    //     auto fallback_average_ring = [&](int i, int j)
    //     {
    //         double sum = 0.0;
    //         int cnt = 0;

    //         for (int k = inner.lo.k; k < inner.hi.k; ++k)
    //         {
    //             sum += U(i, j, k, 0);
    //             ++cnt;
    //         }

    //         if (cnt <= 0)
    //             return;

    //         const double avg = sum / static_cast<double>(cnt);

    //         for (int k = inner.lo.k; k < inner.hi.k; ++k)
    //             U(i, j, k, 0) = avg;
    //     };

    //     auto solve3x3 = [](double A[3][3], double b[3], double x[3]) -> bool
    //     {
    //         double M[3][4] = {
    //             {A[0][0], A[0][1], A[0][2], b[0]},
    //             {A[1][0], A[1][1], A[1][2], b[1]},
    //             {A[2][0], A[2][1], A[2][2], b[2]}};

    //         for (int col = 0; col < 3; ++col)
    //         {
    //             int piv = col;
    //             double maxv = std::abs(M[col][col]);

    //             for (int r = col + 1; r < 3; ++r)
    //             {
    //                 const double v = std::abs(M[r][col]);
    //                 if (v > maxv)
    //                 {
    //                     maxv = v;
    //                     piv = r;
    //                 }
    //             }

    //             if (maxv < 1.0e-300)
    //                 return false;

    //             if (piv != col)
    //             {
    //                 for (int c = col; c < 4; ++c)
    //                     std::swap(M[col][c], M[piv][c]);
    //             }

    //             const double inv = 1.0 / M[col][col];

    //             for (int c = col; c < 4; ++c)
    //                 M[col][c] *= inv;

    //             for (int r = 0; r < 3; ++r)
    //             {
    //                 if (r == col)
    //                     continue;

    //                 const double fac = M[r][col];

    //                 for (int c = col; c < 4; ++c)
    //                     M[r][c] -= fac * M[col][c];
    //             }
    //         }

    //         x[0] = M[0][3];
    //         x[1] = M[1][3];
    //         x[2] = M[2][3];

    //         return true;
    //     };

    //     auto project_induced_axis_ring = [&](int i, int j)
    //     {
    //         if (!Aaxis)
    //         {
    //             fallback_average_ring(i, j);
    //             return;
    //         }

    //         const int klo = inner.lo.k;
    //         const int khi = inner.hi.k;
    //         const int nk = khi - klo;

    //         if (nk <= 0)
    //             return;

    //         double mean_phi = 0.0;
    //         double mean_A[3] = {0.0, 0.0, 0.0};

    //         for (int k = klo; k < khi; ++k)
    //         {
    //             mean_phi += U(i, j, k, 0);

    //             mean_A[0] += (*Aaxis)(i, j, k, 0);
    //             mean_A[1] += (*Aaxis)(i, j, k, 1);
    //             mean_A[2] += (*Aaxis)(i, j, k, 2);
    //         }

    //         const double inv_nk = 1.0 / static_cast<double>(nk);

    //         mean_phi *= inv_nk;
    //         mean_A[0] *= inv_nk;
    //         mean_A[1] *= inv_nk;
    //         mean_A[2] *= inv_nk;

    //         // 求 b:
    //         // min_b sum_k | (phi_k - mean_phi)
    //         //              - (A_k - mean_A) · b |^2
    //         double M[3][3] = {
    //             {0.0, 0.0, 0.0},
    //             {0.0, 0.0, 0.0},
    //             {0.0, 0.0, 0.0}};

    //         double rhs[3] = {0.0, 0.0, 0.0};

    //         for (int k = klo; k < khi; ++k)
    //         {
    //             const double h[3] = {
    //                 (*Aaxis)(i, j, k, 0) - mean_A[0],
    //                 (*Aaxis)(i, j, k, 1) - mean_A[1],
    //                 (*Aaxis)(i, j, k, 2) - mean_A[2]};

    //             const double rphi = U(i, j, k, 0) - mean_phi;

    //             // flux-norm LS. 如果以后想最小化 Bn 误差，
    //             // 可改成 w = 1.0 / (|A|^2 + eps)。
    //             const double w = 1.0;

    //             for (int a = 0; a < 3; ++a)
    //             {
    //                 rhs[a] += w * h[a] * rphi;

    //                 for (int b = 0; b < 3; ++b)
    //                     M[a][b] += w * h[a] * h[b];
    //             }
    //         }

    //         const double tr = M[0][0] + M[1][1] + M[2][2];

    //         double beta[3] = {0.0, 0.0, 0.0};

    //         if (tr > 1.0e-300)
    //         {
    //             // Tikhonov regularization，避免 A_k 近似平行时奇异。
    //             const double reg = 1.0e-10 * tr;

    //             M[0][0] += reg;
    //             M[1][1] += reg;
    //             M[2][2] += reg;

    //             const bool ok = solve3x3(M, rhs, beta);

    //             if (!ok)
    //             {
    //                 beta[0] = 0.0;
    //                 beta[1] = 0.0;
    //                 beta[2] = 0.0;
    //             }
    //         }

    //         // 生成 target induced flux。
    //         // theta=1: 完全使用几何投影；
    //         // theta=0: 退回原始平均。
    //         const double theta = 1.0;

    //         double sum_target = 0.0;
    //         double sum_old = 0.0;

    //         for (int k = klo; k < khi; ++k)
    //         {
    //             const double h[3] = {
    //                 (*Aaxis)(i, j, k, 0) - mean_A[0],
    //                 (*Aaxis)(i, j, k, 1) - mean_A[1],
    //                 (*Aaxis)(i, j, k, 2) - mean_A[2]};

    //             const double geom =
    //                 h[0] * beta[0] +
    //                 h[1] * beta[1] +
    //                 h[2] * beta[2];

    //             const double phi_target =
    //                 mean_phi + theta * geom;

    //             sum_target += phi_target;
    //             sum_old += U(i, j, k, 0);

    //             U(i, j, k, 0) = phi_target;
    //         }

    //         // roundoff correction：严格保持 sum_k U 不变。
    //         const double corr =
    //             (sum_old - sum_target) / static_cast<double>(nk);

    //         for (int k = klo; k < khi; ++k)
    //             U(i, j, k, 0) += corr;
    //     };

    //     for (int i = inner.lo.i; i < inner.hi.i; ++i)
    //     {
    //         for (int j = inner.lo.j; j < inner.hi.j; ++j)
    //         {
    //             project_induced_axis_ring(i, j);
    //         }
    //     }

    //     BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    //     return;
    // }

    if (is_axis_face)
    {
        if (dir == 1)
        {
            // norm=xi, axis=eta:
            // fixed i,j; average over k
            for (int i = inner.lo.i; i < inner.hi.i; ++i)
            {
                for (int j = inner.lo.j; j < inner.hi.j; ++j)
                {
                    double sum = 0.0;
                    int cnt = 0;

                    for (int k = inner.lo.k; k < inner.hi.k; ++k)
                    {
                        sum += U(i, j, k, 0);
                        ++cnt;
                    }

                    if (cnt > 0)
                    {
                        const double avg = sum / static_cast<double>(cnt);
                        for (int k = inner.lo.k; k < inner.hi.k; ++k)
                            U(i, j, k, 0) = avg;
                    }
                }
            }
        }
        else // dir == 2
        {
            // norm=eta, axis=xi:
            // fixed i,j; average over k
            for (int i = inner.lo.i; i < inner.hi.i; ++i)
            {
                for (int j = inner.lo.j; j < inner.hi.j; ++j)
                {
                    double sum = 0.0;
                    int cnt = 0;

                    for (int k = inner.lo.k; k < inner.hi.k; ++k)
                    {
                        sum += U(i, j, k, 0);
                        ++cnt;
                    }

                    if (cnt > 0)
                    {
                        const double avg = sum / static_cast<double>(cnt);
                        for (int k = inner.lo.k; k < inner.hi.k; ++k)
                            U(i, j, k, 0) = avg;
                    }
                }
            }
        }

        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    // ------------------------------------------------------------
    // 3. B_zeta | norm=0 无效：不做 collapse，保持默认 copy。
    // ------------------------------------------------------------
    if (is_Bze)
    {
        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
        return;
    }

    BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
}

void MercuryBoundary::BC_Solid_Surface_Eface_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{

    // const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write
    const Box3 &g = r.inner_slab;

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                U(i, j, k, 0) = 0.0;
                U(i, j, k, 1) = 0.0;
                U(i, j, k, 2) = 0.0;
            }

    BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
}

void MercuryBoundary::BC_Solid_Surface_Eface_ghots_zero(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    bool is_norm = false;
    if (U.descriptor().name == "Eface_xi")
    {
        if (abs(r.direction) == 1)
            is_norm = true;
    }
    else if (U.descriptor().name == "Eface_eta")
    {
        if (abs(r.direction) == 2)
            is_norm = true;
    }
    else if (U.descriptor().name == "Eface_zeta")
    {
        if (abs(r.direction) == 3)
            is_norm = true;
    }
    else
    {
        std::cerr << "BC_Solid_Surface_Eface_ghots_zero called for unexpected field: "
                  << U.descriptor().name << std::endl;
        return;
    }

    if (is_norm)
    {
        const Box3 &g = r.inner_slab;

        for (int i = g.lo.i; i < g.hi.i; ++i)
            for (int j = g.lo.j; j < g.hi.j; ++j)
                for (int k = g.lo.k; k < g.hi.k; ++k)
                {
                    U(i, j, k, 0) = 0.0;
                    U(i, j, k, 1) = 0.0;
                    U(i, j, k, 2) = 0.0;
                }
    }

    const Box3 &g = BoundaryCore::MakeGhostSlabFromInner(r.inner_slab, r.direction, ngh); // ghost slab to write

    for (int i = g.lo.i; i < g.hi.i; ++i)
        for (int j = g.lo.j; j < g.hi.j; ++j)
            for (int k = g.lo.k; k < g.hi.k; ++k)
            {
                U(i, j, k, 0) = 0.0;
                U(i, j, k, 1) = 0.0;
                U(i, j, k, 2) = 0.0;
            }
}

void MercuryBoundary::BC_Solid_Surface_Eedge_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
{
    const Box3 &wall = r.inner_slab;
    const int ax = std::abs(r.direction);
    const int sgn = (r.direction > 0) ? +1 : -1;

    const Int3 outward_step = (ax == 1)   ? Int3{sgn, 0, 0}
                              : (ax == 2) ? Int3{0, sgn, 0}
                                          : Int3{0, 0, sgn};
    // Only tangential E_edge at 1st layer off the wall will be modified
    if (U.descriptor().name == "E_xi" && ax == 1)
        return;
    if (U.descriptor().name == "E_eta" && ax == 2)
        return;
    if (U.descriptor().name == "E_zeta" && ax == 3)
        return;

    const Int3 ulo = U.get_lo();
    const Int3 uhi = U.get_hi();

    for (int i = wall.lo.i; i < wall.hi.i; ++i)
        for (int j = wall.lo.j; j < wall.hi.j; ++j)
            for (int k = wall.lo.k; k < wall.hi.k; ++k)
            {
                const int ii = i - outward_step.i;
                const int jj = j - outward_step.j;
                const int kk = k - outward_step.k;

                U(ii, jj, kk, 0) = 0.0;
            }

    // const Box3 &wall = r.inner_slab;
    // const int ax = std::abs(r.direction);
    // const int sgn = (r.direction > 0) ? +1 : -1;

    // const Int3 outward_step = (ax == 1)   ? Int3{sgn, 0, 0}
    //                           : (ax == 2) ? Int3{0, sgn, 0}
    //                                       : Int3{0, 0, sgn};

    // const Int3 ulo = U.get_lo();
    // const Int3 uhi = U.get_hi();
    // const int buffer_width = 4; // std::max(1, 4 * ngh);

    // auto in_range = [&](int i, int j, int k)
    // {
    //     return (i >= ulo.i && i < uhi.i &&
    //             j >= ulo.j && j < uhi.j &&
    //             k >= ulo.k && k < uhi.k);
    // };

    // for (int i = wall.lo.i; i < wall.hi.i; ++i)
    //     for (int j = wall.lo.j; j < wall.hi.j; ++j)
    //         for (int k = wall.lo.k; k < wall.hi.k; ++k)
    //             for (int d = 0; d <= buffer_width; ++d)
    //             {
    //                 const int ii = i - d * outward_step.i;
    //                 const int jj = j - d * outward_step.j;
    //                 const int kk = k - d * outward_step.k;

    //                 if (!in_range(ii, jj, kk))
    //                     continue;

    //                 const double t = static_cast<double>(d) / static_cast<double>(buffer_width);
    //                 const double damping = t * t * (3.0 - 2.0 * t);
    //                 U(ii, jj, kk, 0) *= damping;
    //             }
}
