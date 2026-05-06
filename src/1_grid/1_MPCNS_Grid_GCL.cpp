#include "0_basic/DEFINE.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "array"
#include "math.h"

void Block::calc_Dual_Grids()
{
    int ngg = jacobi.Getghostmesh();
    // dual grid的i，j，k都向正向偏移半个单位即为正常网格坐标
    dual_x.SetSize(mx + 2 * ngg + 1, my + 2 * ngg + 1, mz + 2 * ngg + 1, ngg);
    dual_y.SetSize(mx + 2 * ngg + 1, my + 2 * ngg + 1, mz + 2 * ngg + 1, ngg);
    dual_z.SetSize(mx + 2 * ngg + 1, my + 2 * ngg + 1, mz + 2 * ngg + 1, ngg);
    for (int i = -ngg; i <= mx + ngg; i++)
        for (int j = -ngg; j <= my + ngg; j++)
            for (int k = -ngg; k <= mz + ngg; k++)
            {
                dual_x(i, j, k) = 0.125 * (x(i, j, k) + x(i - 1, j, k) + x(i, j - 1, k) + x(i - 1, j - 1, k) +
                                           x(i, j, k - 1) + x(i - 1, j, k - 1) + x(i, j - 1, k - 1) + x(i - 1, j - 1, k - 1));
                dual_y(i, j, k) = 0.125 * (y(i, j, k) + y(i - 1, j, k) + y(i, j - 1, k) + y(i - 1, j - 1, k) +
                                           y(i, j, k - 1) + y(i - 1, j, k - 1) + y(i, j - 1, k - 1) + y(i - 1, j - 1, k - 1));
                dual_z(i, j, k) = 0.125 * (z(i, j, k) + z(i - 1, j, k) + z(i, j - 1, k) + z(i - 1, j - 1, k) +
                                           z(i, j, k - 1) + z(i - 1, j, k - 1) + z(i, j - 1, k - 1) + z(i - 1, j - 1, k - 1));
            }
}

// 本函数生成：每个面的协变切矢 a_α (不要单位化) 以及协变度规 g_αβ = a_α · a_β
// 基于对偶网格计算
void Block::calc_Face_Tangent_Vectors(int32_t &ngg)
{
    auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    // const double eps = 1e-300;

    // --- 目标数组：协变切矢（向量）与协变度规（3x3 张量），全部在相应的面心 ---
    GCL_covar_xi_xi.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_xi_eta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_xi_zeta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_xi.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3, 3);

    GCL_covar_eta_xi.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_eta_eta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_eta_zeta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_eta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3, 3);

    GCL_covar_zeta_xi.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_zeta_eta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_zeta_zeta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_covar_zeta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3, 3);

    // ========== 3D 主路径 ==========
    if (dimension == 3)
    {
        // -------- ξ 面：索引 (i in [-1..mx+1], j in [-1..my+1], k in [-1..mz+1]) --------
        // 针对ξ 面，i=-1表示i=-1/2半点的位置
        for (int i = -2; i <= mx + 1; ++i)
            for (int j = -1; j <= my + 1; ++j)
                for (int k = -1; k <= mz + 1; ++k)
                {
                    // xi方向直接使用网格差分计算
                    std::array<double, 3> r_xi = {
                        x(i + 1, j, k) - x(i, j, k),
                        y(i + 1, j, k) - y(i, j, k),
                        z(i + 1, j, k) - z(i, j, k)};
                    // eta方向使用dual网格差分计算
                    std::array<double, 3> r_eta = {
                        0.5 * (dual_x(i + 1, j + 1, k + 1) + dual_x(i + 1, j + 1, k) - dual_x(i + 1, j, k + 1) - dual_x(i + 1, j, k)),
                        0.5 * (dual_y(i + 1, j + 1, k + 1) + dual_y(i + 1, j + 1, k) - dual_y(i + 1, j, k + 1) - dual_y(i + 1, j, k)),
                        0.5 * (dual_z(i + 1, j + 1, k + 1) + dual_z(i + 1, j + 1, k) - dual_z(i + 1, j, k + 1) - dual_z(i + 1, j, k))};
                    std::array<double, 3> r_zeta = {
                        0.5 * (dual_x(i + 1, j + 1, k + 1) + dual_x(i + 1, j, k + 1) - dual_x(i + 1, j + 1, k) - dual_x(i + 1, j, k)),
                        0.5 * (dual_y(i + 1, j + 1, k + 1) + dual_y(i + 1, j, k + 1) - dual_y(i + 1, j + 1, k) - dual_y(i + 1, j, k)),
                        0.5 * (dual_z(i + 1, j + 1, k + 1) + dual_z(i + 1, j, k + 1) - dual_z(i + 1, j + 1, k) - dual_z(i + 1, j, k))};

                    // 存协变切矢
                    for (int c = 0; c < 3; ++c)
                    {
                        GCL_covar_xi_xi(i, j, k, c) = r_xi[c];
                        GCL_covar_xi_eta(i, j, k, c) = r_eta[c];
                        GCL_covar_xi_zeta(i, j, k, c) = r_zeta[c];
                    }

                    // 协变度规 g_αβ = a_α · a_β （α,β ∈ {ξ,η,ζ}，对称填充）
                    double g_xx = dot(r_xi, r_xi);
                    double g_xe = dot(r_xi, r_eta);
                    double g_xz = dot(r_xi, r_zeta);
                    double g_ee = dot(r_eta, r_eta);
                    double g_ez = dot(r_eta, r_zeta);
                    double g_zz = dot(r_zeta, r_zeta);

                    GCL_covar_xi(i, j, k, 0, 0) = g_xx;
                    GCL_covar_xi(i, j, k, 0, 1) = g_xe;
                    GCL_covar_xi(i, j, k, 0, 2) = g_xz;
                    GCL_covar_xi(i, j, k, 1, 0) = g_xe;
                    GCL_covar_xi(i, j, k, 1, 1) = g_ee;
                    GCL_covar_xi(i, j, k, 1, 2) = g_ez;
                    GCL_covar_xi(i, j, k, 2, 0) = g_xz;
                    GCL_covar_xi(i, j, k, 2, 1) = g_ez;
                    GCL_covar_xi(i, j, k, 2, 2) = g_zz;
                }

        // -------- η 面：索引 (i in [-1..mx+1], j in [-1..my+1], k in [-1..mz+1]) --------
        // 针对η 面，j=-1表示j=-1/2半点的位置
        for (int i = -1; i <= mx + 1; ++i)
            for (int j = -2; j <= my + 1; ++j)
                for (int k = -1; k <= mz + 1; ++k)
                {
                    // xi方向使用dual网格差分计算
                    std::array<double, 3> r_xi = {
                        0.5 * (dual_x(i + 1, j + 1, k + 1) + dual_x(i + 1, j + 1, k) - dual_x(i, j + 1, k + 1) - dual_x(i, j + 1, k)),
                        0.5 * (dual_y(i + 1, j + 1, k + 1) + dual_y(i + 1, j + 1, k) - dual_y(i, j + 1, k + 1) - dual_y(i, j + 1, k)),
                        0.5 * (dual_z(i + 1, j + 1, k + 1) + dual_z(i + 1, j + 1, k) - dual_z(i, j + 1, k + 1) - dual_z(i, j + 1, k))};

                    // eta方向直接使用网格差分计算
                    std::array<double, 3> r_eta = {
                        x(i, j + 1, k) - x(i, j, k),
                        y(i, j + 1, k) - y(i, j, k),
                        z(i, j + 1, k) - z(i, j, k)};

                    // zeta方向使用dual网格差分计算
                    std::array<double, 3> r_zeta = {
                        0.5 * (dual_x(i + 1, j + 1, k + 1) + dual_x(i, j + 1, k + 1) - dual_x(i + 1, j + 1, k) - dual_x(i, j + 1, k)),
                        0.5 * (dual_y(i + 1, j + 1, k + 1) + dual_y(i, j + 1, k + 1) - dual_y(i + 1, j + 1, k) - dual_y(i, j + 1, k)),
                        0.5 * (dual_z(i + 1, j + 1, k + 1) + dual_z(i, j + 1, k + 1) - dual_z(i + 1, j + 1, k) - dual_z(i, j + 1, k))};

                    // 存协变切矢
                    for (int c = 0; c < 3; ++c)
                    {
                        GCL_covar_eta_xi(i, j, k, c) = r_xi[c];
                        GCL_covar_eta_eta(i, j, k, c) = r_eta[c];
                        GCL_covar_eta_zeta(i, j, k, c) = r_zeta[c];
                    }

                    // 协变度规 g_αβ = a_α · a_β （α,β ∈ {ξ,η,ζ}，对称填充）
                    double g_xx = dot(r_xi, r_xi);
                    double g_xe = dot(r_xi, r_eta);
                    double g_xz = dot(r_xi, r_zeta);
                    double g_ee = dot(r_eta, r_eta);
                    double g_ez = dot(r_eta, r_zeta);
                    double g_zz = dot(r_zeta, r_zeta);

                    GCL_covar_eta(i, j, k, 0, 0) = g_xx;
                    GCL_covar_eta(i, j, k, 0, 1) = g_xe;
                    GCL_covar_eta(i, j, k, 0, 2) = g_xz;
                    GCL_covar_eta(i, j, k, 1, 0) = g_xe;
                    GCL_covar_eta(i, j, k, 1, 1) = g_ee;
                    GCL_covar_eta(i, j, k, 1, 2) = g_ez;
                    GCL_covar_eta(i, j, k, 2, 0) = g_xz;
                    GCL_covar_eta(i, j, k, 2, 1) = g_ez;
                    GCL_covar_eta(i, j, k, 2, 2) = g_zz;
                }

        // -------- ζ 面：索引 (i in [-1..mx+1], j in [-1..my+1], k in [-1..mz+1]) --------
        // 针对ζ 面，k=-1表示k=-1/2半点的位置
        for (int i = -1; i <= mx + 1; ++i)
            for (int j = -1; j <= my + 1; ++j)
                for (int k = -2; k <= mz + 1; ++k)
                {
                    // xi方向使用dual网格差分计算
                    std::array<double, 3> r_xi = {
                        0.5 * (dual_x(i + 1, j + 1, k + 1) + dual_x(i + 1, j, k + 1) - dual_x(i, j, k + 1) - dual_x(i, j + 1, k + 1)),
                        0.5 * (dual_y(i + 1, j + 1, k + 1) + dual_y(i + 1, j, k + 1) - dual_y(i, j, k + 1) - dual_y(i, j + 1, k + 1)),
                        0.5 * (dual_z(i + 1, j + 1, k + 1) + dual_z(i + 1, j, k + 1) - dual_z(i, j, k + 1) - dual_z(i, j + 1, k + 1))};

                    // eta方向使用dual网格差分计算
                    std::array<double, 3> r_eta = {
                        0.5 * (dual_x(i + 1, j + 1, k + 1) + dual_x(i, j + 1, k + 1) - dual_x(i + 1, j, k + 1) - dual_x(i, j, k + 1)),
                        0.5 * (dual_y(i + 1, j + 1, k + 1) + dual_y(i, j + 1, k + 1) - dual_y(i + 1, j, k + 1) - dual_y(i, j, k + 1)),
                        0.5 * (dual_z(i + 1, j + 1, k + 1) + dual_z(i, j + 1, k + 1) - dual_z(i + 1, j, k + 1) - dual_z(i, j, k + 1))};

                    // zeta方向直接使用网格差分计算
                    std::array<double, 3> r_zeta = {
                        x(i, j, k + 1) - x(i, j, k),
                        y(i, j, k + 1) - y(i, j, k),
                        z(i, j, k + 1) - z(i, j, k)};

                    // 存协变切矢
                    for (int c = 0; c < 3; ++c)
                    {
                        GCL_covar_zeta_xi(i, j, k, c) = r_xi[c];
                        GCL_covar_zeta_eta(i, j, k, c) = r_eta[c];
                        GCL_covar_zeta_zeta(i, j, k, c) = r_zeta[c];
                    }

                    // 协变度规 g_αβ = a_α · a_β （α,β ∈ {ξ,η,ζ}，对称填充）
                    double g_xx = dot(r_xi, r_xi);
                    double g_xe = dot(r_xi, r_eta);
                    double g_xz = dot(r_xi, r_zeta);
                    double g_ee = dot(r_eta, r_eta);
                    double g_ez = dot(r_eta, r_zeta);
                    double g_zz = dot(r_zeta, r_zeta);

                    GCL_covar_zeta(i, j, k, 0, 0) = g_xx;
                    GCL_covar_zeta(i, j, k, 0, 1) = g_xe;
                    GCL_covar_zeta(i, j, k, 0, 2) = g_xz;
                    GCL_covar_zeta(i, j, k, 1, 0) = g_xe;
                    GCL_covar_zeta(i, j, k, 1, 1) = g_ee;
                    GCL_covar_zeta(i, j, k, 1, 2) = g_ez;
                    GCL_covar_zeta(i, j, k, 2, 0) = g_xz;
                    GCL_covar_zeta(i, j, k, 2, 1) = g_ez;
                    GCL_covar_zeta(i, j, k, 2, 2) = g_zz;
                }
    }
    else // ========== 2D 退化（ζ 无效）：用面内中心差生成协变切矢 ==========
    {
        // To be developed
    }
}

/**
 * @brief 计算Jacob行列式、度量系数
 * @remark 本程序依赖ngg参数,GCL_metric_*，表示 (J ∇ξ), (J ∇η), (J ∇ζ) 的“面量矢量”，全在半点/面心存储。
 */
void Block::calc_Metrics_GCL(int32_t &ngg)
{
    // ---- helpers -----------------------------------------------------------
    auto cross = [](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return std::array<double, 3>{
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
    };

    auto Pdual = [&](int i, int j, int k)
    {
        return std::array<double, 3>{dual_x(i, j, k), dual_y(i, j, k), dual_z(i, j, k)};
    };
    // 四角面面积向量（两三角求和）：对退化天然鲁棒
    auto quad_area = [&](const std::array<double, 3> &P00,
                         const std::array<double, 3> &P10,
                         const std::array<double, 3> &P01,
                         const std::array<double, 3> &P11)
    {
        std::array<double, 3> r1{P10[0] - P00[0], P10[1] - P00[1], P10[2] - P00[2]};
        std::array<double, 3> r2{P01[0] - P00[0], P01[1] - P00[1], P01[2] - P00[2]};
        std::array<double, 3> r3{P11[0] - P10[0], P11[1] - P10[1], P11[2] - P10[2]};
        std::array<double, 3> r4{P01[0] - P10[0], P01[1] - P10[1], P01[2] - P10[2]};
        auto A1 = cross(r1, r2);
        auto A2 = cross(r3, r4);
        return std::array<double, 3>{0.5 * (A1[0] + A2[0]),
                                     0.5 * (A1[1] + A2[1]),
                                     0.5 * (A1[2] + A2[2])};
    };

    // ---- storage -----------------------------------------------------------
    GCL_metric_xi.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_metric_eta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);
    GCL_metric_zeta.SetSize(5 + mx, 5 + my, 5 + mz, 2, 3);

    if (dimension == 3)
    {
        // -------------------- 半点面构造（严格右手） --------------------
        // 目标：确保 cross(r1, r2) 的法向就是 +ξ / +η / +ζ
        // i+1/2 面：用 dual 的  (i+1,*,*)
        auto xi_face_area_half = [&](int i, int j, int k)
        {
            // P00: (j,k), P10: (j+1,k)  → r1 = +η
            // P01: (j,k+1)               → r2 = +ζ
            // P11: (j+1,k+1)
            return quad_area(Pdual(i + 1, j, k),          // P00
                             Pdual(i + 1, j + 1, k),      // P10  (+η)
                             Pdual(i + 1, j, k + 1),      // P01  (+ζ)
                             Pdual(i + 1, j + 1, k + 1)); // P11
        };
        // j+1/2 面
        auto eta_face_area_half = [&](int i, int j, int k)
        {
            // 注意这里与常见 (ξ,ζ) 顺序相反，故意选择 (+ζ)×(+ξ)
            // P00: (i,k), P10: (i,k+1)  → r1 = +ζ
            // P01: (i+1,k)              → r2 = +ξ
            // P11: (i+1,k+1)
            return quad_area(Pdual(i, j + 1, k),          // P00
                             Pdual(i, j + 1, k + 1),      // P10  (+ζ)
                             Pdual(i + 1, j + 1, k),      // P01  (+ξ)
                             Pdual(i + 1, j + 1, k + 1)); // P11
        };
        // k+1/2 面
        auto zeta_face_area_half = [&](int i, int j, int k)
        {
            // P00: (i,j), P10: (i+1,j)  → r1 = +ξ
            // P01: (i,j+1)              → r2 = +η
            // P11: (i+1,j+1)
            return quad_area(Pdual(i, j, k + 1),          // P00
                             Pdual(i + 1, j, k + 1),      // P10  (+ξ)
                             Pdual(i, j + 1, k + 1),      // P01  (+η)
                             Pdual(i + 1, j + 1, k + 1)); // P11
        };

        // ---------- (J ∇ξ) at i+1/2 ----------
        for (int i = -2; i <= mx + 1; ++i)
            for (int j = -1; j <= my + 1; ++j)
                for (int k = -1; k <= mz + 1; ++k)
                {
                    auto A = xi_face_area_half(i, j, k);
                    GCL_metric_xi(i, j, k, 0) = A[0];
                    GCL_metric_xi(i, j, k, 1) = A[1];
                    GCL_metric_xi(i, j, k, 2) = A[2];
                }

        // ---------- (J ∇η) at j+1/2 ----------
        for (int i = -1; i <= mx + 1; ++i)
            for (int j = -2; j <= my + 1; ++j)
                for (int k = -1; k <= mz + 1; ++k)
                {
                    auto A = eta_face_area_half(i, j, k);
                    GCL_metric_eta(i, j, k, 0) = A[0];
                    GCL_metric_eta(i, j, k, 1) = A[1];
                    GCL_metric_eta(i, j, k, 2) = A[2];
                }

        // ---------- (J ∇ζ) at k+1/2 ----------
        for (int i = -1; i <= mx + 1; ++i)
            for (int j = -1; j <= my + 1; ++j)
                for (int k = -2; k <= mz + 1; ++k)
                {
                    auto A = zeta_face_area_half(i, j, k);
                    GCL_metric_zeta(i, j, k, 0) = A[0];
                    GCL_metric_zeta(i, j, k, 1) = A[1];
                    GCL_metric_zeta(i, j, k, 2) = A[2];
                }
    }
    else // dimension == 2
    {
        // 2D: 只有 ξ/η；可把“面”理解为“线段法向长度向量”
        // 如果需要我可以提供与 dual 的 2D 统一实现（本段留空）
    }

    // （可选）恒等式投影：多数情况下可关闭；若要开启，把 for 次数设为 >0
    auto sweep_identity_projection = [&]()
    {
        for (int i = -1; i <= mx; ++i)
            for (int j = -1; j <= my; ++j)
                for (int k = -1; k <= mz; ++k)
                {
                    const double Rx =
                        (GCL_metric_xi(i, j, k, 0) - GCL_metric_xi(i - 1, j, k, 0)) +
                        (GCL_metric_eta(i, j, k, 0) - GCL_metric_eta(i, j - 1, k, 0)) +
                        ((dimension == 3) ? (GCL_metric_zeta(i, j, k, 0) - GCL_metric_zeta(i, j, k - 1, 0)) : 0.0);
                    const double Ry =
                        (GCL_metric_xi(i, j, k, 1) - GCL_metric_xi(i - 1, j, k, 1)) +
                        (GCL_metric_eta(i, j, k, 1) - GCL_metric_eta(i, j - 1, k, 1)) +
                        ((dimension == 3) ? (GCL_metric_zeta(i, j, k, 1) - GCL_metric_zeta(i, j, k - 1, 1)) : 0.0);
                    const double Rz =
                        (GCL_metric_xi(i, j, k, 2) - GCL_metric_xi(i - 1, j, k, 2)) +
                        (GCL_metric_eta(i, j, k, 2) - GCL_metric_eta(i, j - 1, k, 2)) +
                        ((dimension == 3) ? (GCL_metric_zeta(i, j, k, 2) - GCL_metric_zeta(i, j, k - 1, 2)) : 0.0);

                    const double invD = 1.0 / double(dimension);
                    const double fx = Rx * invD, fy = Ry * invD, fz = Rz * invD;

                    GCL_metric_xi(i, j, k, 0) -= fx;
                    GCL_metric_xi(i, j, k, 1) -= fy;
                    GCL_metric_xi(i, j, k, 2) -= fz;

                    GCL_metric_eta(i, j, k, 0) -= fx;
                    GCL_metric_eta(i, j, k, 1) -= fy;
                    GCL_metric_eta(i, j, k, 2) -= fz;

                    if (dimension == 3)
                    {
                        GCL_metric_zeta(i, j, k, 0) -= fx;
                        GCL_metric_zeta(i, j, k, 1) -= fy;
                        GCL_metric_zeta(i, j, k, 2) -= fz;
                    }
                }
    };
    for (int s = 0; s < 0; ++s) // 开启请把 0 改为 1 或更多迭代
        sweep_identity_projection();
}

void Block::calc_modify_Jacobi(int32_t &ngg)
{
    auto face_center = [&](const std::array<double, 3> &P00,
                           const std::array<double, 3> &P10,
                           const std::array<double, 3> &P01,
                           const std::array<double, 3> &P11)
    {
        return std::array<double, 3>{
            0.25 * (P00[0] + P10[0] + P01[0] + P11[0]),
            0.25 * (P00[1] + P10[1] + P01[1] + P11[1]),
            0.25 * (P00[2] + P10[2] + P01[2] + P11[2])};
    };

    auto dot3 = [](const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    auto Pdual = [&](int i, int j, int k)
    {
        return std::array<double, 3>{dual_x(i, j, k), dual_y(i, j, k), dual_z(i, j, k)};
    };

    if (dimension == 3)
    {
        // 计算守恒的jacobi
        for (int i = 0; i <= mx; ++i)
            for (int j = 0; j <= my; ++j)
                for (int k = 0; k <= mz; ++k)
                {
                    double V = 0.0;

                    // 1. xi- 面 (i-1/2)
                    {
                        std::array<double, 3> A = {GCL_metric_xi(i - 1, j, k, 0),
                                                   GCL_metric_xi(i - 1, j, k, 1),
                                                   GCL_metric_xi(i - 1, j, k, 2)};
                        auto P00 = Pdual(i, j, k);
                        auto P10 = Pdual(i, j + 1, k);
                        auto P01 = Pdual(i, j, k + 1);
                        auto P11 = Pdual(i, j + 1, k + 1);
                        auto rc = face_center(P00, P10, P01, P11);
                        V -= dot3(A, rc);
                    }

                    // 2. xi+ 面 (i+1/2)
                    {
                        std::array<double, 3> A = {GCL_metric_xi(i, j, k, 0),
                                                   GCL_metric_xi(i, j, k, 1),
                                                   GCL_metric_xi(i, j, k, 2)};
                        auto P00 = Pdual(i + 1, j, k);
                        auto P10 = Pdual(i + 1, j + 1, k);
                        auto P01 = Pdual(i + 1, j, k + 1);
                        auto P11 = Pdual(i + 1, j + 1, k + 1);
                        auto rc = face_center(P00, P10, P01, P11);
                        V += dot3(A, rc);
                    }

                    // 3. eta- 面 (j-1/2)
                    {
                        std::array<double, 3> A = {GCL_metric_eta(i, j - 1, k, 0),
                                                   GCL_metric_eta(i, j - 1, k, 1),
                                                   GCL_metric_eta(i, j - 1, k, 2)};
                        auto P00 = Pdual(i, j, k);
                        auto P10 = Pdual(i + 1, j, k);
                        auto P01 = Pdual(i, j, k + 1);
                        auto P11 = Pdual(i + 1, j, k + 1);
                        auto rc = face_center(P00, P10, P01, P11);
                        V -= dot3(A, rc);
                    }

                    // 4. eta+ 面 (j+1/2)
                    {
                        std::array<double, 3> A = {GCL_metric_eta(i, j, k, 0),
                                                   GCL_metric_eta(i, j, k, 1),
                                                   GCL_metric_eta(i, j, k, 2)};
                        auto P00 = Pdual(i, j + 1, k);
                        auto P10 = Pdual(i + 1, j + 1, k);
                        auto P01 = Pdual(i, j + 1, k + 1);
                        auto P11 = Pdual(i + 1, j + 1, k + 1);
                        auto rc = face_center(P00, P10, P01, P11);
                        V += dot3(A, rc);
                    }

                    // 5. zeta- 面 (k-1/2)
                    {
                        std::array<double, 3> A = {GCL_metric_zeta(i, j, k - 1, 0),
                                                   GCL_metric_zeta(i, j, k - 1, 1),
                                                   GCL_metric_zeta(i, j, k - 1, 2)};
                        auto P00 = Pdual(i, j, k);
                        auto P10 = Pdual(i + 1, j, k);
                        auto P01 = Pdual(i, j + 1, k);
                        auto P11 = Pdual(i + 1, j + 1, k);
                        auto rc = face_center(P00, P10, P01, P11);
                        V -= dot3(A, rc);
                    }

                    // 6. zeta+ 面 (k+1/2)
                    {
                        std::array<double, 3> A = {GCL_metric_zeta(i, j, k, 0),
                                                   GCL_metric_zeta(i, j, k, 1),
                                                   GCL_metric_zeta(i, j, k, 2)};
                        auto P00 = Pdual(i, j, k + 1);
                        auto P10 = Pdual(i + 1, j, k + 1);
                        auto P01 = Pdual(i, j + 1, k + 1);
                        auto P11 = Pdual(i + 1, j + 1, k + 1);
                        auto rc = face_center(P00, P10, P01, P11);
                        V += dot3(A, rc);
                    }

                    jacobi(i, j, k) = V / 3.0;
                    // // --- 中心几何导数 r_ξ, r_η, r_ζ（节点坐标中心差） ---
                    // // 需要坐标 ghost 至少 ±1
                    // double rxi_x = 0.5 * (x(i + 1, j, k) - x(i - 1, j, k));
                    // double rxi_y = 0.5 * (y(i + 1, j, k) - y(i - 1, j, k));
                    // double rxi_z = 0.5 * (z(i + 1, j, k) - z(i - 1, j, k));

                    // double reta_x = 0.5 * (x(i, j + 1, k) - x(i, j - 1, k));
                    // double reta_y = 0.5 * (y(i, j + 1, k) - y(i, j - 1, k));
                    // double reta_z = 0.5 * (z(i, j + 1, k) - z(i, j - 1, k));

                    // double rzet_x = 0.5 * (x(i, j, k + 1) - x(i, j, k - 1));
                    // double rzet_y = 0.5 * (y(i, j, k + 1) - y(i, j, k - 1));
                    // double rzet_z = 0.5 * (z(i, j, k + 1) - z(i, j, k - 1));

                    // // --- 面度量在单元中心的平均 (Ja·∇ξ)^c, (Ja·∇η)^c, (Ja·∇ζ)^c ---
                    // // 注意索引：xi 面存储在 i-1/2 与 i+1/2 两侧 → 用 (i-1) 与 (i) 两个面平均；eta/zeta 类似
                    // double Ja_xi_x = 0.5 * (GCL_metric_xi(i - 1, j, k, 0) + GCL_metric_xi(i, j, k, 0));
                    // double Ja_xi_y = 0.5 * (GCL_metric_xi(i - 1, j, k, 1) + GCL_metric_xi(i, j, k, 1));
                    // double Ja_xi_z = 0.5 * (GCL_metric_xi(i - 1, j, k, 2) + GCL_metric_xi(i, j, k, 2));

                    // double Ja_et_x = 0.5 * (GCL_metric_eta(i, j - 1, k, 0) + GCL_metric_eta(i, j, k, 0));
                    // double Ja_et_y = 0.5 * (GCL_metric_eta(i, j - 1, k, 1) + GCL_metric_eta(i, j, k, 1));
                    // double Ja_et_z = 0.5 * (GCL_metric_eta(i, j - 1, k, 2) + GCL_metric_eta(i, j, k, 2));

                    // double Ja_zt_x = 0.5 * (GCL_metric_zeta(i, j, k - 1, 0) + GCL_metric_zeta(i, j, k, 0));
                    // double Ja_zt_y = 0.5 * (GCL_metric_zeta(i, j, k - 1, 1) + GCL_metric_zeta(i, j, k, 1));
                    // double Ja_zt_z = 0.5 * (GCL_metric_zeta(i, j, k - 1, 2) + GCL_metric_zeta(i, j, k, 2));

                    // // --- 计算中心 J（同构公式） ---
                    // double Jc =
                    //     (rxi_x * Ja_xi_x + rxi_y * Ja_xi_y + rxi_z * Ja_xi_z + reta_x * Ja_et_x + reta_y * Ja_et_y + reta_z * Ja_et_z + rzet_x * Ja_zt_x + rzet_y * Ja_zt_y + rzet_z * Ja_zt_z) / 3.0;

                    // jacobi(i, j, k) = Jc;
                }
    }
    else
    {
        // 计算守恒的jacobi
        // for (int i = 0; i <= mx; ++i)
        //     for (int j = 0; j <= my; ++j)
        //         for (int k = 0; k <= mz; ++k)
        //         {
        //             // --- 中心几何导数 r_ξ, r_η, r_ζ（节点坐标中心差） ---
        //             // 需要坐标 ghost 至少 ±1
        //             double rxi_x = 0.5 * (x(i + 1, j, k) - x(i - 1, j, k));
        //             double rxi_y = 0.5 * (y(i + 1, j, k) - y(i - 1, j, k));
        //             double rxi_z = 0.5 * (z(i + 1, j, k) - z(i - 1, j, k));

        //             double reta_x = 0.5 * (x(i, j + 1, k) - x(i, j - 1, k));
        //             double reta_y = 0.5 * (y(i, j + 1, k) - y(i, j - 1, k));
        //             double reta_z = 0.5 * (z(i, j + 1, k) - z(i, j - 1, k));

        //             double rzet_x = 0.5 * (x(i, j, k + 1) - x(i, j, k - 1));
        //             double rzet_y = 0.5 * (y(i, j, k + 1) - y(i, j, k - 1));
        //             double rzet_z = 0.5 * (z(i, j, k + 1) - z(i, j, k - 1));

        //             // --- 面度量在单元中心的平均 (Ja·∇ξ)^c, (Ja·∇η)^c, (Ja·∇ζ)^c ---
        //             // 注意索引：xi 面存储在 i-1/2 与 i+1/2 两侧 → 用 (i-1) 与 (i) 两个面平均；eta/zeta 类似
        //             double Ja_xi_x = 0.5 * (GCL_metric_xi(i - 1, j, k, 0) + GCL_metric_xi(i, j, k, 0));
        //             double Ja_xi_y = 0.5 * (GCL_metric_xi(i - 1, j, k, 1) + GCL_metric_xi(i, j, k, 1));
        //             double Ja_xi_z = 0.5 * (GCL_metric_xi(i - 1, j, k, 2) + GCL_metric_xi(i, j, k, 2));

        //             double Ja_et_x = 0.5 * (GCL_metric_eta(i, j - 1, k, 0) + GCL_metric_eta(i, j, k, 0));
        //             double Ja_et_y = 0.5 * (GCL_metric_eta(i, j - 1, k, 1) + GCL_metric_eta(i, j, k, 1));
        //             double Ja_et_z = 0.5 * (GCL_metric_eta(i, j - 1, k, 2) + GCL_metric_eta(i, j, k, 2));

        //             // --- 计算中心 J（同构公式） ---
        //             double Jc =
        //                 ((rxi_x * Ja_xi_x + rxi_y * Ja_xi_y + rxi_z * Ja_xi_z) +
        //                  (reta_x * Ja_et_x + reta_y * Ja_et_y + reta_z * Ja_et_z)) /
        //                 2.0;

        //             jacobi(i, j, k) = Jc;
        //         }
    }
}
/**
 * @brief 循环所有block，给所有的Polar面/轴添加虚网格坐标;
 */
void Grid::MeshTrans_Inner_Pole()
{
    for (int i = 0; i < nblock; i++)
    {
        int num_phy = grids(i).physical_bc.size();
        for (int j = 0; j < num_phy; j++)
        {
            if (grids(i).physical_bc[j].boundary_name == "Pole")
            {
                // 要求Pole旋成方向为维数方向，即2维是\eta, 3维是\zeta, 且不可剖分
                int sub[3], sup[3], cycle[3];
                sub[0] = fmin(fabs(grids(i).physical_bc[j].sub[0]), fabs(grids(i).physical_bc[j].sup[0]));
                sub[1] = fmin(fabs(grids(i).physical_bc[j].sub[1]), fabs(grids(i).physical_bc[j].sup[1]));
                sub[2] = fmin(fabs(grids(i).physical_bc[j].sub[2]), fabs(grids(i).physical_bc[j].sup[2]));
                sup[0] = fmax(fabs(grids(i).physical_bc[j].sub[0]), fabs(grids(i).physical_bc[j].sup[0]));
                sup[1] = fmax(fabs(grids(i).physical_bc[j].sub[1]), fabs(grids(i).physical_bc[j].sup[1]));
                sup[2] = fmax(fabs(grids(i).physical_bc[j].sub[2]), fabs(grids(i).physical_bc[j].sup[2]));
                cycle[0] = grids(i).physical_bc[j].cycle[0];
                cycle[1] = grids(i).physical_bc[j].cycle[1];
                cycle[2] = grids(i).physical_bc[j].cycle[2];
                if (dimension == 2)
                {
                    if (sub[1] != 0 || sup[1] != grids(i).my || fmod(grids(i).my, 2) != 0)
                    {
                        std::cout << "Fatal Error! ! ! Pole requires eta (2D) or zeta (3D) is rotational direction, and could not be split! my (or  mz in 3D) should be even" << std::endl;
                        exit(-1);
                    }
                    int my = grids(i).my;
                    for (int index = 1; index <= ngg + 1; index++)
                        for (int ii = sub[0] + index * cycle[0]; ii <= sup[0] + index * cycle[0]; ii++)
                            for (int jj = sub[1] + index * cycle[1]; jj <= sup[1] + index * cycle[1]; jj++)
                                for (int kk = sub[2] + index * cycle[2]; kk <= sup[2] + index * cycle[2]; kk++)
                                {
                                    grids(i).x(ii, jj, kk) = grids(i).x(ii - cycle[0] * index, jj, kk);
                                    grids(i).y(ii, jj, kk) = grids(i).y(ii - cycle[0] * index, jj, kk);
                                    grids(i).z(ii, jj, kk) = grids(i).z(ii - cycle[0] * index, jj, kk);
                                }
                }
                else if (dimension == 3)
                {
                    if (sub[2] != 0 || sup[2] != grids(i).mz || fmod(grids(i).mz, 2) != 0)
                    {
                        std::cout << "Fatal Error! ! ! Pole requires eta (2D) or zeta (3D) is rotational direction, and could not be split!  my (or  mz in 3D) should be even" << std::endl;
                        exit(-1);
                    }
                    int mz = grids(i).mz;
                    for (int index = 1; index <= ngg + 1; index++)
                        for (int ii = sub[0] + index * cycle[0]; ii <= sup[0] + index * cycle[0]; ii++)
                            for (int jj = sub[1] + index * cycle[1]; jj <= sup[1] + index * cycle[1]; jj++)
                                for (int kk = sub[2] + index * cycle[2]; kk <= sup[2] + index * cycle[2]; kk++)
                                {
                                    grids(i).x(ii, jj, kk) = grids(i).x(ii - cycle[0] * index, jj - cycle[1] * index, kk);
                                    grids(i).y(ii, jj, kk) = grids(i).y(ii - cycle[0] * index, jj - cycle[1] * index, kk);
                                    grids(i).z(ii, jj, kk) = grids(i).z(ii - cycle[0] * index, jj - cycle[1] * index, kk);
                                }
                }
            }
        }
    }
}