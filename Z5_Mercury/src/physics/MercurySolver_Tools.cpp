#include "MercurySolver.h"

void MercurySolver::calc_physical_constant(Param *par)
{
    // -------- constants / reference ----------
    auto cst = par_->GetDou_List("constant");
    auto ref = par_->GetDou_List("REF");

    gamma_ = cst.data["gamma"];
    R_uni = cst.data["R_uni"];

    NA = cst.data["NA"];
    q_e = cst.data["q_e"]; // e (Coulomb)
    k_Boltz = R_uni / NA;
    mu0 = cst.data["mu_mag"];           // μ0
    U_ref = ref.data.at("U");           // m/s
    L_ref = ref.data.at("L_ref");       // m
    B_ref = ref.data.at("B_ref");       // Telsa
    T_ref = ref.data.at("T");           // K
    n_ref = ref.data["n"];              // 1/m^3
    M_ref = ref.data["Molecular_mass"]; // kg/mol (molar mass)
    rho_ref = M_ref * n_ref / NA;       // kg/m^3
    M_H = par->GetDou("mole_mass1");    // kg/mol
    M_Na = par->GetDou("mole_mass2");   // kg/mol
    m_H = M_H / NA;                     // kg/particle
    m_Na = M_Na / NA;                   // kg/particle

    // 无量纲状态方程系数：p = rho * T * coeff
    // coeff = (R_uni * T_ref) / (M * U_ref^2)
    state_coeff_H = (R_uni * T_ref) / (M_H * U_ref * U_ref);
    state_coeff_Na = (R_uni * T_ref) / (M_Na * U_ref * U_ref);

    CFL = par_->GetDou("CFL");

    hall_coef = B_ref * M_ref / (U_ref * q_e * L_ref * mu0 * rho_ref * NA);

    // ambi_coef = M_H * U_ref / (cst.data["q_e"] * ref.data["L_ref"] * ref.data["B_ref"]);

    momentum_induce_coeff = (q_e * L_ref * B_ref * n_ref) / (rho_ref * U_ref); // incude_coeff * n_ns * (u_ns - u_+) \times B  = momentum eqs source, n_ns: m^-3
    momentum_hall_coeff = (B_ref * B_ref) / (mu0 * rho_ref * U_ref * U_ref);   // momentum_hall_coeff * n_ns / n_e * \nabla\times B\times B = momentum eqs source

    inver_MA2 = ref.data["B_ref"] * ref.data["B_ref"] / (U_ref * U_ref * cst.data["mu_mag"] * rho_ref);

    inver_Rem = par->GetDou("eta_max_mercury") / (cst.data["mu_mag"] * U_ref * ref.data["L_ref"]);

    ne_hall_floor = 0.01;                              // nondimensional
    ne_hall_floor_dimensional = ne_hall_floor * n_ref; // dimension
}

// void MercurySolver::Hall_Num_Limiter(double rhoH, double rhoNa, double *num)
// {
//     // const double nH_m = (rho_ref * rhoH) / m_H;    // #/m^3
//     // const double nNa_m = (rho_ref * rhoNa) / m_Na; // #/m^3

//     const double rhoH_pos = std::max(rhoH, 0.0);
//     const double rhoNa_pos = std::max(rhoNa, 0.0);

//     const double nH_ = rhoH_pos / M_H * M_ref;    // non-dimensional
//     const double nNa_ = rhoNa_pos / M_Na * M_ref; // non-dimensional
//     const double ne_true = nH_ + nNa_;

//     // 平滑 floor 后的总电子数密度
//     const double ne_lim = std::sqrt(ne_true * ne_true +
//                                     ne_hall_floor * ne_hall_floor);

//     double nH_lim = 0.0;
//     double nNa_lim = 0.0;

//     if (ne_true > 1e-300)
//     {
//         // 统一相似缩放：保持组分比例不变，同时保证 ne = nH + nNa
//         const double scale = ne_lim / ne_true;
//         nH_lim = scale * nH_;
//         nNa_lim = scale * nNa_;
//     }
//     else
//     {
//         // 真空/极低密区的退化情形：
//         // 这里比例本来就没有定义，给一个兜底组成。
//         // 最简单先全给 H；如果你有更合适的远场组成，也可改成远场比例。
//         nH_lim = ne_lim;
//         nNa_lim = 0.0;
//     }

//     num[0] = nH_lim;
//     num[1] = nNa_lim;
//     num[2] = ne_lim;
// }

NumInfo MercurySolver::Hall_Num_Limiter(double rhoH, double rhoNa)
{
    // const double nH_m = (rho_ref * rhoH) / m_H;    // #/m^3
    // const double nNa_m = (rho_ref * rhoNa) / m_Na; // #/m^3
    // const double rhoH_pos = std::max(rhoH, 0.0);
    // const double rhoNa_pos = std::max(rhoNa, 0.0);
    // const double nH_ = rhoH_pos / M_H * M_ref;    // non-dimensional
    // const double nNa_ = rhoNa_pos / M_Na * M_ref; // non-dimensional
    // const double ne_true = nH_ + nNa_;

    NumInfo nd;

    const double rhoH_pos = std::max(rhoH, 0.0);
    const double rhoNa_pos = std::max(rhoNa, 0.0);

    nd.nH_true = rhoH_pos / M_H * M_ref;
    nd.nNa_true = rhoNa_pos / M_Na * M_ref;
    nd.ne_true = nd.nH_true + nd.nNa_true;

    nd.ne_eff = std::sqrt(nd.ne_true * nd.ne_true + ne_hall_floor * ne_hall_floor);

    // 这里只防止 0/0；不要用 ne_hall_floor 作为 chi 的分母。
    // chi 是组分比例，不是 Hall regularization。
    constexpr double ne_ratio_tiny = 1.0e-300;

    if (nd.ne_true > ne_ratio_tiny)
    {
        nd.chiH = nd.nH_true / nd.ne_true;
        nd.chiNa = nd.nNa_true / nd.ne_true;
    }
    else
    {
        // 真空处比例无定义。这里保持旧逻辑：退化为 H。
        // 由于速度计算里 rho 很小时 uH/uNa 会置 0，所以 u_plus 通常仍为 0。
        nd.chiH = 1.0;
        nd.chiNa = 0.0;
    }

    // 可选：低密度 MHD 源项 taper。
    // 正常区 ne_true >> ne_floor: w ~= chi
    // 稀薄区 ne_true << ne_floor: w -> 0
    nd.wH_mhd = nd.nH_true / nd.ne_eff;
    nd.wNa_mhd = nd.nNa_true / nd.ne_eff;
    nd.mhd_taper = nd.wH_mhd + nd.wNa_mhd;

    return nd;
}

void MercurySolver::calc_PV()
{
    const double rho_floor = 1e-12;
    const double p_floor = 1e-12;

    auto fill_one = [&](int fidU, int fidPV, double coeff)
    {
        const int nb = fld_->num_blocks();
        for (int ib = 0; ib < nb; ++ib)
        {
            FieldBlock &U = fld_->field(fidU, ib);
            FieldBlock &PV = fld_->field(fidPV, ib);
            if (!U.is_allocated() || !PV.is_allocated())
                continue;

            const Int3 lo = PV.get_lo(); // 含 ghost：[-ng, ...]
            const Int3 hi = PV.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double rho0 = U(i, j, k, 0);
                        const double rho = std::max(rho0, rho_floor);

                        const double u = U(i, j, k, 1) / rho;
                        const double v = U(i, j, k, 2) / rho;
                        const double w = U(i, j, k, 3) / rho;

                        const double E = U(i, j, k, 4);
                        const double ke = 0.5 * rho * (u * u + v * v + w * w);

                        // 压力：p = (gamma-1) * (E - ke)
                        double eint = E - ke;
                        if (eint < 0.0)
                            eint = 0.0;
                        double p = (gamma_ - 1.0) * eint;
                        p = std::max(p, p_floor);

                        // 温度：T = p / (rho * coeff)
                        // 对 Na 用 coeff_Na -> 等价于 Fortran 的 T = p/(n kB)，其中 n= rho/m_Na
                        double T = p / (rho * coeff);

                        PV(i, j, k, 0) = u;
                        PV(i, j, k, 1) = v;
                        PV(i, j, k, 2) = w;
                        PV(i, j, k, 3) = p;
                        PV(i, j, k, 4) = T;
                    }
        }
    };

    fill_one(fid_.fid_U_H, fid_.fid_PV_H, state_coeff_H);
    fill_one(fid_.fid_U_Na, fid_.fid_PV_Na, state_coeff_Na);
}

void MercurySolver::calc_Uplus()
{
    const double rho_eps = 1e-20;
    const double inv23 = M_H / M_Na; // 与 Fortran sm2≈23*sm1 对齐

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UN = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &Up = fld_->field(fid_.fid_U_plus, ib); // 新增

        if (!UH.is_allocated() || !UN.is_allocated() || !Up.is_allocated())
            continue;

        const Int3 lo = Up.get_lo();
        const Int3 hi = Up.get_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double rhoH0 = std::max(UH(i, j, k, 0), 0.0);
                    const double rhoNa0 = std::max(UN(i, j, k, 0), 0.0);

                    const double uH = (rhoH0 > rho_eps) ? UH(i, j, k, 1) / rhoH0 : 0.0;
                    const double vH = (rhoH0 > rho_eps) ? UH(i, j, k, 2) / rhoH0 : 0.0;
                    const double wH = (rhoH0 > rho_eps) ? UH(i, j, k, 3) / rhoH0 : 0.0;

                    const double uNa = (rhoNa0 > rho_eps) ? UN(i, j, k, 1) / rhoNa0 : 0.0;
                    const double vNa = (rhoNa0 > rho_eps) ? UN(i, j, k, 2) / rhoNa0 : 0.0;
                    const double wNa = (rhoNa0 > rho_eps) ? UN(i, j, k, 3) / rhoNa0 : 0.0;

                    // double num[3];
                    // Hall_Num_Limiter(rhoH0, rhoNa0, num)
                    // const double nH = num[0];
                    // const double nNa = num[1];
                    // const double ne = num[2];
                    // Up(i, j, k, 0) = (nH * uH + nNa * uNa) / ne;
                    // Up(i, j, k, 1) = (nH * vH + nNa * vNa) / ne;
                    // Up(i, j, k, 2) = (nH * wH + nNa * wNa) / ne;

                    NumInfo num = Hall_Num_Limiter(rhoH0, rhoNa0);
                    const double chi_H = num.chiH;
                    const double chi_Na = num.chiNa;
                    Up(i, j, k, 0) = chi_H * uH + chi_Na * uNa;
                    Up(i, j, k, 1) = chi_H * vH + chi_Na * vNa;
                    Up(i, j, k, 2) = chi_H * wH + chi_Na * wNa;

                    // const double nH = rhoH0;
                    // const double nNa = rhoNa0 * inv23;
                    // const double nt = nH + nNa;

                    // if (nt <= 0.0)
                    // {
                    //     Up(i, j, k, 0) = 0.0;
                    //     Up(i, j, k, 1) = 0.0;
                    //     Up(i, j, k, 2) = 0.0;
                    // }
                    // else
                    // {
                    //     Up(i, j, k, 0) = (nH * uH + nNa * uNa) / nt;
                    //     Up(i, j, k, 1) = (nH * vH + nNa * vNa) / nt;
                    //     Up(i, j, k, 2) = (nH * wH + nNa * wNa) / nt;
                    // }
                }
    }
}

void MercurySolver::calc_Bcell()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &Bcell = fld_->field(fid_.fid_Bcell, ib);
        auto &Bindcell = fld_->field(fid_.fid_Bindcell, ib);

        auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
        auto &Beta = fld_->field(fid_.fid_B.eta, ib);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);

        auto &Baddxi = fld_->field(fid_.fid_Badd.xi, ib);
        auto &Baddeta = fld_->field(fid_.fid_Badd.eta, ib);
        auto &Baddzeta = fld_->field(fid_.fid_Badd.zeta, ib);

        auto &W = hall_face_scratch_[ib].dBcell_w;

        if (!Bcell.is_allocated() || !Bindcell.is_allocated() ||
            !Bxi.is_allocated() || !Beta.is_allocated() || !Bzeta.is_allocated())
            continue;

        Int3 lo = Bcell.inner_lo();
        Int3 hi = Bcell.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    // -------- total B = B + Badd --------
                    const double phi_tot[6] = {
                        -Bxi(i, j, k, 0) - Baddxi(i, j, k, 0),           // xi-
                        Bxi(i + 1, j, k, 0) + Baddxi(i + 1, j, k, 0),    // xi+
                        -Beta(i, j, k, 0) - Baddeta(i, j, k, 0),         // eta-
                        Beta(i, j + 1, k, 0) + Baddeta(i, j + 1, k, 0),  // eta+
                        -Bzeta(i, j, k, 0) - Baddzeta(i, j, k, 0),       // zeta-
                        Bzeta(i, j, k + 1, 0) + Baddzeta(i, j, k + 1, 0) // zeta+
                    };

                    // -------- induced B = B only --------
                    const double phi_ind[6] = {
                        -Bxi(i, j, k, 0),
                        Bxi(i + 1, j, k, 0),
                        -Beta(i, j, k, 0),
                        Beta(i, j + 1, k, 0),
                        -Bzeta(i, j, k, 0),
                        Bzeta(i, j, k + 1, 0)};

                    double Bx_tot = 0.0, By_tot = 0.0, Bz_tot = 0.0;
                    double Bx_ind = 0.0, By_ind = 0.0, Bz_ind = 0.0;

                    for (int n = 0; n < 6; ++n)
                    {
                        const double wt_x = W(i, j, k, n);
                        const double wt_y = W(i, j, k, 6 + n);
                        const double wt_z = W(i, j, k, 12 + n);

                        Bx_tot += wt_x * phi_tot[n];
                        By_tot += wt_y * phi_tot[n];
                        Bz_tot += wt_z * phi_tot[n];

                        Bx_ind += wt_x * phi_ind[n];
                        By_ind += wt_y * phi_ind[n];
                        Bz_ind += wt_z * phi_ind[n];
                    }

                    Bcell(i, j, k, 0) = Bx_tot;
                    Bcell(i, j, k, 1) = By_tot;
                    Bcell(i, j, k, 2) = Bz_tot;

                    Bindcell(i, j, k, 0) = Bx_ind;
                    Bindcell(i, j, k, 1) = By_ind;
                    Bindcell(i, j, k, 2) = Bz_ind;
                }
    }

    mercury_bound_.Sync("B_cell");

    // const int nblock = fld_->num_blocks();

    // const double eps = 1e-300;
    // const double delta = 1e-15; // same spirit as your reference

    // auto dot = [&](const std::array<double, 3> &a, const std::array<double, 3> &b)
    // {
    //     return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    // };
    // auto norm = [&](const std::array<double, 3> &a)
    // {
    //     return std::sqrt(dot(a, a));
    // };

    // for (int ib = 0; ib < nblock; ++ib)
    // {
    //     auto &Bcell = fld_->field(fid_.fid_Bcell, ib);
    //     auto &Bindcell = fld_->field(fid_.fid_Bindcell, ib);
    //     // auto &U = fld_->field(fid_.fid_U, ib);

    //     auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
    //     auto &Beta = fld_->field(fid_.fid_B.eta, ib);
    //     auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);

    //     auto &Baddxi = fld_->field(fid_.fid_Badd.xi, ib);
    //     auto &Baddeta = fld_->field(fid_.fid_Badd.eta, ib);
    //     auto &Baddzeta = fld_->field(fid_.fid_Badd.zeta, ib);

    //     auto &Jac = fld_->field(fid_.fid_Jac, ib);
    //     auto &A_xi = fld_->field(fid_.fid_metric.xi, ib);   // JDxi
    //     auto &A_eta = fld_->field(fid_.fid_metric.eta, ib); // JDet
    //     auto &A_ze = fld_->field(fid_.fid_metric.zeta, ib); // JDze

    //     auto &x = grd_->grids(ib).x;
    //     auto &y = grd_->grids(ib).y;
    //     auto &z = grd_->grids(ib).z;

    //     auto &cx = grd_->grids(ib).dual_x;
    //     auto &cy = grd_->grids(ib).dual_y;
    //     auto &cz = grd_->grids(ib).dual_z;

    //     // face centers (consistent with your Jac construction style)
    //     auto xfc_xi = [&](int i, int j, int k) -> std::array<double, 3>
    //     {
    //         return {
    //             0.25 * (x(i, j, k) + x(i, j + 1, k) + x(i, j, k + 1) + x(i, j + 1, k + 1)),
    //             0.25 * (y(i, j, k) + y(i, j + 1, k) + y(i, j, k + 1) + y(i, j + 1, k + 1)),
    //             0.25 * (z(i, j, k) + z(i, j + 1, k) + z(i, j, k + 1) + z(i, j + 1, k + 1))};
    //     };
    //     auto xfc_eta = [&](int i, int j, int k) -> std::array<double, 3>
    //     {
    //         return {
    //             0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j, k + 1) + x(i + 1, j, k + 1)),
    //             0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j, k + 1) + y(i + 1, j, k + 1)),
    //             0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j, k + 1) + z(i + 1, j, k + 1))};
    //     };
    //     auto xfc_zeta = [&](int i, int j, int k) -> std::array<double, 3>
    //     {
    //         return {
    //             0.25 * (x(i, j, k) + x(i + 1, j, k) + x(i, j + 1, k) + x(i + 1, j + 1, k)),
    //             0.25 * (y(i, j, k) + y(i + 1, j, k) + y(i, j + 1, k) + y(i + 1, j + 1, k)),
    //             0.25 * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k))};
    //     };

    //     Int3 lo = Jac.inner_lo();
    //     Int3 hi = Jac.inner_hi();

    //     for (int i = lo.i; i < hi.i; ++i)
    //         for (int j = lo.j; j < hi.j; ++j)
    //             for (int k = lo.k; k < hi.k; ++k)
    //             {
    //                 // cell center
    //                 std::array<double, 3> Xc = {
    //                     cx(i + 1, j + 1, k + 1),
    //                     cy(i + 1, j + 1, k + 1),
    //                     cz(i + 1, j + 1, k + 1)};

    //                 struct FaceEq
    //                 {
    //                     std::array<double, 3> n; // normalized S
    //                     double phi;              // normalized Phi
    //                     double w;                // weight
    //                 };
    //                 FaceEq eqs[6];
    //                 int K = 0;

    //                 auto push = [&](const std::array<double, 3> &S,
    //                                 double Phi,
    //                                 const std::array<double, 3> &Xf)
    //                 {
    //                     double s_norm = norm(S) + eps;
    //                     std::array<double, 3> nvec = {S[0] / s_norm, S[1] / s_norm, S[2] / s_norm};
    //                     double phi_hat = Phi / s_norm;

    //                     double dx = Xf[0] - Xc[0];
    //                     double dy = Xf[1] - Xc[1];
    //                     double dz = Xf[2] - Xc[2];
    //                     double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

    //                     double w = 1.0 / std::sqrt(dist * dist + delta * delta);

    //                     eqs[K++] = {nvec, phi_hat, w};
    //                 };

    //                 // ---------- xi- face (at i) ----------
    //                 std::array<double, 3> S_xm = {
    //                     -A_xi(i, j, k, 0),
    //                     -A_xi(i, j, k, 1),
    //                     -A_xi(i, j, k, 2)};
    //                 double Phi_xm = -Bxi(i, j, k, 0) - Baddxi(i, j, k, 0);
    //                 push(S_xm, Phi_xm, xfc_xi(i, j, k));

    //                 // ---------- xi+ face (at i+1) ----------
    //                 std::array<double, 3> S_xp = {
    //                     A_xi(i + 1, j, k, 0),
    //                     A_xi(i + 1, j, k, 1),
    //                     A_xi(i + 1, j, k, 2)};
    //                 double Phi_xp = Bxi(i + 1, j, k, 0) + Baddxi(i + 1, j, k, 0);
    //                 push(S_xp, Phi_xp, xfc_xi(i + 1, j, k));

    //                 // ---------- eta- face (at j) ----------
    //                 std::array<double, 3> S_em = {
    //                     -A_eta(i, j, k, 0),
    //                     -A_eta(i, j, k, 1),
    //                     -A_eta(i, j, k, 2)};
    //                 double Phi_em = -Beta(i, j, k, 0) - Baddeta(i, j, k, 0);
    //                 push(S_em, Phi_em, xfc_eta(i, j, k));

    //                 // ---------- eta+ face (at j+1) ----------
    //                 std::array<double, 3> S_ep = {
    //                     A_eta(i, j + 1, k, 0),
    //                     A_eta(i, j + 1, k, 1),
    //                     A_eta(i, j + 1, k, 2)};
    //                 double Phi_ep = Beta(i, j + 1, k, 0) + Baddeta(i, j + 1, k, 0);
    //                 push(S_ep, Phi_ep, xfc_eta(i, j + 1, k));

    //                 // ---------- zeta- face (at k) ----------
    //                 std::array<double, 3> S_zm = {
    //                     -A_ze(i, j, k, 0),
    //                     -A_ze(i, j, k, 1),
    //                     -A_ze(i, j, k, 2)};
    //                 double Phi_zm = -Bzeta(i, j, k, 0) - Baddzeta(i, j, k, 0);
    //                 push(S_zm, Phi_zm, xfc_zeta(i, j, k));

    //                 // ---------- zeta+ face (at k+1) ----------
    //                 std::array<double, 3> S_zp = {
    //                     A_ze(i, j, k + 1, 0),
    //                     A_ze(i, j, k + 1, 1),
    //                     A_ze(i, j, k + 1, 2)};
    //                 double Phi_zp = Bzeta(i, j, k + 1, 0) + Baddzeta(i, j, k + 1, 0);
    //                 push(S_zp, Phi_zp, xfc_zeta(i, j, k + 1));

    //                 // ---------- build normal equations N = A^T W A, r = A^T W phi ----------
    //                 double N00 = 0, N01 = 0, N02 = 0, N11 = 0, N12 = 0, N22 = 0;
    //                 double rx = 0, ry = 0, rz = 0;

    //                 for (int t = 0; t < K; ++t)
    //                 {
    //                     double w = eqs[t].w;
    //                     const auto &n = eqs[t].n;
    //                     double phi = eqs[t].phi;

    //                     N00 += w * n[0] * n[0];
    //                     N01 += w * n[0] * n[1];
    //                     N02 += w * n[0] * n[2];
    //                     N11 += w * n[1] * n[1];
    //                     N12 += w * n[1] * n[2];
    //                     N22 += w * n[2] * n[2];

    //                     rx += w * phi * n[0];
    //                     ry += w * phi * n[1];
    //                     rz += w * phi * n[2];
    //                 }

    //                 auto det3 = [&](double a, double b, double c, double d, double e, double f)
    //                 {
    //                     // | a b c |
    //                     // | b d e |
    //                     // | c e f |
    //                     return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
    //                 };

    //                 double det = det3(N00, N01, N02, N11, N12, N22);
    //                 double reg = 1e-14 * (N00 + N11 + N22);

    //                 if (std::abs(det) < reg)
    //                 {
    //                     N00 += reg;
    //                     N11 += reg;
    //                     N22 += reg;
    //                     det = det3(N00, N01, N02, N11, N12, N22);
    //                 }

    //                 // cofactors of symmetric matrix
    //                 double C00 = (N11 * N22 - N12 * N12);
    //                 double C01 = (N02 * N12 - N01 * N22);
    //                 double C02 = (N01 * N12 - N02 * N11);
    //                 double C11 = (N00 * N22 - N02 * N02);
    //                 double C12 = (N01 * N02 - N00 * N12);
    //                 double C22 = (N00 * N11 - N01 * N01);

    //                 double inv = 1.0 / det;

    //                 double Bx_tot = inv * (C00 * rx + C01 * ry + C02 * rz);
    //                 double By_tot = inv * (C01 * rx + C11 * ry + C12 * rz);
    //                 double Bz_tot = inv * (C02 * rx + C12 * ry + C22 * rz);

    //                 // // energy consistency (keep your existing style)
    //                 // const double Bx_old = Bcell(i, j, k, 0);
    //                 // const double By_old = Bcell(i, j, k, 1);
    //                 // const double Bz_old = Bcell(i, j, k, 2);
    //                 // const double Delta_Eb =
    //                 //     0.5 * inver_MA2 * (Bx_tot * Bx_tot + By_tot * By_tot + Bz_tot * Bz_tot) -
    //                 //     0.5 * inver_MA2 * (Bx_old * Bx_old + By_old * By_old + Bz_old * Bz_old);

    //                 Bcell(i, j, k, 0) = Bx_tot;
    //                 Bcell(i, j, k, 1) = By_tot;
    //                 Bcell(i, j, k, 2) = Bz_tot;
    //             }

    //     for (int i = lo.i; i < hi.i; ++i)
    //         for (int j = lo.j; j < hi.j; ++j)
    //             for (int k = lo.k; k < hi.k; ++k)
    //             {
    //                 // cell center
    //                 std::array<double, 3> Xc = {
    //                     cx(i + 1, j + 1, k + 1),
    //                     cy(i + 1, j + 1, k + 1),
    //                     cz(i + 1, j + 1, k + 1)};

    //                 struct FaceEq
    //                 {
    //                     std::array<double, 3> n; // normalized S
    //                     double phi;              // normalized Phi
    //                     double w;                // weight
    //                 };
    //                 FaceEq eqs[6];
    //                 int K = 0;

    //                 auto push = [&](const std::array<double, 3> &S,
    //                                 double Phi,
    //                                 const std::array<double, 3> &Xf)
    //                 {
    //                     double s_norm = norm(S) + eps;
    //                     std::array<double, 3> nvec = {S[0] / s_norm, S[1] / s_norm, S[2] / s_norm};
    //                     double phi_hat = Phi / s_norm;

    //                     double dx = Xf[0] - Xc[0];
    //                     double dy = Xf[1] - Xc[1];
    //                     double dz = Xf[2] - Xc[2];
    //                     double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

    //                     double w = 1.0 / std::sqrt(dist * dist + delta * delta);

    //                     eqs[K++] = {nvec, phi_hat, w};
    //                 };

    //                 // ---------- xi- face (at i) ----------
    //                 std::array<double, 3> S_xm = {
    //                     -A_xi(i, j, k, 0),
    //                     -A_xi(i, j, k, 1),
    //                     -A_xi(i, j, k, 2)};
    //                 double Phi_xm = -Bxi(i, j, k, 0);
    //                 push(S_xm, Phi_xm, xfc_xi(i, j, k));

    //                 // ---------- xi+ face (at i+1) ----------
    //                 std::array<double, 3> S_xp = {
    //                     A_xi(i + 1, j, k, 0),
    //                     A_xi(i + 1, j, k, 1),
    //                     A_xi(i + 1, j, k, 2)};
    //                 double Phi_xp = Bxi(i + 1, j, k, 0);
    //                 push(S_xp, Phi_xp, xfc_xi(i + 1, j, k));

    //                 // ---------- eta- face (at j) ----------
    //                 std::array<double, 3> S_em = {
    //                     -A_eta(i, j, k, 0),
    //                     -A_eta(i, j, k, 1),
    //                     -A_eta(i, j, k, 2)};
    //                 double Phi_em = -Beta(i, j, k, 0);
    //                 push(S_em, Phi_em, xfc_eta(i, j, k));

    //                 // ---------- eta+ face (at j+1) ----------
    //                 std::array<double, 3> S_ep = {
    //                     A_eta(i, j + 1, k, 0),
    //                     A_eta(i, j + 1, k, 1),
    //                     A_eta(i, j + 1, k, 2)};
    //                 double Phi_ep = Beta(i, j + 1, k, 0);
    //                 push(S_ep, Phi_ep, xfc_eta(i, j + 1, k));

    //                 // ---------- zeta- face (at k) ----------
    //                 std::array<double, 3> S_zm = {
    //                     -A_ze(i, j, k, 0),
    //                     -A_ze(i, j, k, 1),
    //                     -A_ze(i, j, k, 2)};
    //                 double Phi_zm = -Bzeta(i, j, k, 0);
    //                 push(S_zm, Phi_zm, xfc_zeta(i, j, k));

    //                 // ---------- zeta+ face (at k+1) ----------
    //                 std::array<double, 3> S_zp = {
    //                     A_ze(i, j, k + 1, 0),
    //                     A_ze(i, j, k + 1, 1),
    //                     A_ze(i, j, k + 1, 2)};
    //                 double Phi_zp = Bzeta(i, j, k + 1, 0);
    //                 push(S_zp, Phi_zp, xfc_zeta(i, j, k + 1));

    //                 // ---------- build normal equations N = A^T W A, r = A^T W phi ----------
    //                 double N00 = 0, N01 = 0, N02 = 0, N11 = 0, N12 = 0, N22 = 0;
    //                 double rx = 0, ry = 0, rz = 0;

    //                 for (int t = 0; t < K; ++t)
    //                 {
    //                     double w = eqs[t].w;
    //                     const auto &n = eqs[t].n;
    //                     double phi = eqs[t].phi;

    //                     N00 += w * n[0] * n[0];
    //                     N01 += w * n[0] * n[1];
    //                     N02 += w * n[0] * n[2];
    //                     N11 += w * n[1] * n[1];
    //                     N12 += w * n[1] * n[2];
    //                     N22 += w * n[2] * n[2];

    //                     rx += w * phi * n[0];
    //                     ry += w * phi * n[1];
    //                     rz += w * phi * n[2];
    //                 }

    //                 auto det3 = [&](double a, double b, double c, double d, double e, double f)
    //                 {
    //                     // | a b c |
    //                     // | b d e |
    //                     // | c e f |
    //                     return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
    //                 };

    //                 double det = det3(N00, N01, N02, N11, N12, N22);
    //                 double reg = 1e-14 * (N00 + N11 + N22);

    //                 if (std::abs(det) < reg)
    //                 {
    //                     N00 += reg;
    //                     N11 += reg;
    //                     N22 += reg;
    //                     det = det3(N00, N01, N02, N11, N12, N22);
    //                 }

    //                 // cofactors of symmetric matrix
    //                 double C00 = (N11 * N22 - N12 * N12);
    //                 double C01 = (N02 * N12 - N01 * N22);
    //                 double C02 = (N01 * N12 - N02 * N11);
    //                 double C11 = (N00 * N22 - N02 * N02);
    //                 double C12 = (N01 * N02 - N00 * N12);
    //                 double C22 = (N00 * N11 - N01 * N01);

    //                 double inv = 1.0 / det;

    //                 double Bx_tot = inv * (C00 * rx + C01 * ry + C02 * rz);
    //                 double By_tot = inv * (C01 * rx + C11 * ry + C12 * rz);
    //                 double Bz_tot = inv * (C02 * rx + C12 * ry + C22 * rz);

    //                 // // energy consistency (keep your existing style)
    //                 // const double Bx_old = Bcell(i, j, k, 0);
    //                 // const double By_old = Bcell(i, j, k, 1);
    //                 // const double Bz_old = Bcell(i, j, k, 2);
    //                 // const double Delta_Eb =
    //                 //     0.5 * inver_MA2 * (Bx_tot * Bx_tot + By_tot * By_tot + Bz_tot * Bz_tot) -
    //                 //     0.5 * inver_MA2 * (Bx_old * Bx_old + By_old * By_old + Bz_old * Bz_old);

    //                 Bindcell(i, j, k, 0) = Bx_tot;
    //                 Bindcell(i, j, k, 1) = By_tot;
    //                 Bindcell(i, j, k, 2) = Bz_tot;
    //             }
    // }
    // mercury_bound_.Sync("B_cell");
}

void MercurySolver::calc_divB()
{
    const int nblock = fld_->num_blocks();

    auto get_comp = [](const Int3 &a, int ax) -> int
    {
        if (ax == 0)
            return a.i;
        if (ax == 1)
            return a.j;
        return a.k;
    };

    auto set_comp = [](Int3 &a, int ax, int v)
    {
        if (ax == 0)
            a.i = v;
        else if (ax == 1)
            a.j = v;
        else
            a.k = v;
    };

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &divB = fld_->field(fid_.fid_divB, ib);
        auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
        auto &Beta = fld_->field(fid_.fid_B.eta, ib);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);
        auto &Jac = fld_->field(fid_.fid_Jac, ib);

        Int3 lo = divB.inner_lo();
        Int3 hi = divB.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    divB(i, j, k, 0) = (Bxi(i + 1, j, k, 0) - Bxi(i, j, k, 0) +
                                        Beta(i, j + 1, k, 0) - Beta(i, j, k, 0) +
                                        Bzeta(i, j, k + 1, 0) - Bzeta(i, j, k, 0)) /
                                       Jac(i, j, k, 0);
                }

        // Pole collapse intentionally identifies degenerate face fluxes, so the
        // adjacent first off-axis shell is not a meaningful divB diagnostic.
        for (const auto &p : topo_->physical_patches)
        {
            if (p.this_block != ib)
                continue;
            if (p.bc_name != "Pole")
                continue;

            const int ax = std::abs(p.direction) - 1;
            const bool high_side = (p.direction > 0);

            Int3 plo = lo;
            Int3 phi = hi;

            for (int d = 0; d < 3; ++d)
            {
                if (d == ax)
                    continue;

                const int tlo = std::max(get_comp(lo, d), get_comp(p.this_box_node.lo, d));
                const int thi = std::min(get_comp(hi, d), get_comp(p.this_box_node.hi, d) - 1);
                set_comp(plo, d, tlo);
                set_comp(phi, d, thi);
            }

            if (high_side)
            {
                set_comp(plo, ax, get_comp(hi, ax) - 1);
                set_comp(phi, ax, get_comp(hi, ax));
            }
            else
            {
                set_comp(plo, ax, get_comp(lo, ax));
                set_comp(phi, ax, get_comp(lo, ax) + 1);
            }

            if (!(plo.i < phi.i && plo.j < phi.j && plo.k < phi.k))
                continue;

            for (int i = plo.i; i < phi.i; ++i)
                for (int j = plo.j; j < phi.j; ++j)
                    for (int k = plo.k; k < phi.k; ++k)
                        divB(i, j, k, 0) = 0.0;
        }
    }
}

void MercurySolver::calc_Jcell()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &Jcell = fld_->field(fid_.fid_Jcell, ib);

        auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
        auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);

        auto &W = hall_face_scratch_[ib].dJcell_w;

        if (!Jcell.is_allocated() || !Jxi.is_allocated() ||
            !Jeta.is_allocated() || !Jzeta.is_allocated())
            continue;

        Int3 lo = Jcell.inner_lo();
        Int3 hi = Jcell.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double s0 = Jxi(i, j, k, 0);
                    const double s1 = Jxi(i, j + 1, k, 0);
                    const double s2 = Jxi(i, j, k + 1, 0);
                    const double s3 = Jxi(i, j + 1, k + 1, 0);

                    const double s4 = Jeta(i, j, k, 0);
                    const double s5 = Jeta(i + 1, j, k, 0);
                    const double s6 = Jeta(i, j, k + 1, 0);
                    const double s7 = Jeta(i + 1, j, k + 1, 0);

                    const double s8 = Jzeta(i, j, k, 0);
                    const double s9 = Jzeta(i + 1, j, k, 0);
                    const double s10 = Jzeta(i, j + 1, k, 0);
                    const double s11 = Jzeta(i + 1, j + 1, k, 0);

                    const double s[12] = {s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11};

                    double Jx = 0.0, Jy = 0.0, Jz = 0.0;
                    for (int n = 0; n < 12; ++n)
                    {
                        const double sn = s[n];
                        Jx += W(i, j, k, n) * sn;
                        Jy += W(i, j, k, 12 + n) * sn;
                        Jz += W(i, j, k, 24 + n) * sn;
                    }

                    Jcell(i, j, k, 0) = Jx;
                    Jcell(i, j, k, 1) = Jy;
                    Jcell(i, j, k, 2) = Jz;
                }
    }

    mercury_bound_.Sync("J_cell");

    // const int nblock = fld_->num_blocks();

    // constexpr double eps = 1e-25;

    // for (int ib = 0; ib < nblock; ++ib)
    // {
    //     auto &Jcell = fld_->field(fid_.fid_Jcell, ib);

    //     auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
    //     auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
    //     auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);

    //     auto &dl_xi = fld_->field("dl_xi", ib);
    //     auto &dl_eta = fld_->field("dl_eta", ib);
    //     auto &dl_zeta = fld_->field("dl_zeta", ib);

    //     auto &x = grd_->grids(ib).x;
    //     auto &y = grd_->grids(ib).y;
    //     auto &z = grd_->grids(ib).z;

    //     if (!Jcell.is_allocated() || !Jxi.is_allocated() || !Jeta.is_allocated() || !Jzeta.is_allocated())
    //         continue;

    //     auto dot3 = [&](double ax, double ay, double az,
    //                     double bx, double by, double bz) -> double
    //     {
    //         return ax * bx + ay * by + az * bz;
    //     };

    //     auto unit_t_xi = [&](int i, int j, int k,
    //                          double &tx, double &ty, double &tz)
    //     {
    //         const double L = std::max(dl_xi(i, j, k, 0), eps);
    //         tx = (x(i + 1, j, k) - x(i, j, k)) / L;
    //         ty = (y(i + 1, j, k) - y(i, j, k)) / L;
    //         tz = (z(i + 1, j, k) - z(i, j, k)) / L;
    //     };

    //     auto unit_t_eta = [&](int i, int j, int k,
    //                           double &tx, double &ty, double &tz)
    //     {
    //         const double L = std::max(dl_eta(i, j, k, 0), eps);
    //         tx = (x(i, j + 1, k) - x(i, j, k)) / L;
    //         ty = (y(i, j + 1, k) - y(i, j, k)) / L;
    //         tz = (z(i, j + 1, k) - z(i, j, k)) / L;
    //     };

    //     auto unit_t_zeta = [&](int i, int j, int k,
    //                            double &tx, double &ty, double &tz)
    //     {
    //         const double L = std::max(dl_zeta(i, j, k, 0), eps);
    //         tx = (x(i, j, k + 1) - x(i, j, k)) / L;
    //         ty = (y(i, j, k + 1) - y(i, j, k)) / L;
    //         tz = (z(i, j, k + 1) - z(i, j, k)) / L;
    //     };

    //     struct Eq
    //     {
    //         double tx, ty, tz; // unit tangent
    //         double rhs;        // J_edge / |dl|
    //         double w;          // weight
    //     };

    //     Int3 lo = Jcell.inner_lo();
    //     Int3 hi = Jcell.inner_hi();

    //     for (int i = lo.i; i < hi.i; ++i)
    //         for (int j = lo.j; j < hi.j; ++j)
    //             for (int k = lo.k; k < hi.k; ++k)
    //             {
    //                 Eq eqs[12];
    //                 int K = 0;

    //                 auto push = [&](double tx, double ty, double tz,
    //                                 double Jint, double L, double w = 1.0)
    //                 {
    //                     L = std::max(L, eps);
    //                     eqs[K++] = {tx, ty, tz, Jint / L, w};
    //                 };

    //                 double tx, ty, tz;

    //                 // =====================================================
    //                 // 4 xi-edges around cell(i,j,k)
    //                 // =====================================================
    //                 unit_t_xi(i, j, k, tx, ty, tz);
    //                 push(tx, ty, tz, Jxi(i, j, k, 0), dl_xi(i, j, k, 0));

    //                 unit_t_xi(i, j + 1, k, tx, ty, tz);
    //                 push(tx, ty, tz, Jxi(i, j + 1, k, 0), dl_xi(i, j + 1, k, 0));

    //                 unit_t_xi(i, j, k + 1, tx, ty, tz);
    //                 push(tx, ty, tz, Jxi(i, j, k + 1, 0), dl_xi(i, j, k + 1, 0));

    //                 unit_t_xi(i, j + 1, k + 1, tx, ty, tz);
    //                 push(tx, ty, tz, Jxi(i, j + 1, k + 1, 0), dl_xi(i, j + 1, k + 1, 0));

    //                 // =====================================================
    //                 // 4 eta-edges
    //                 // =====================================================
    //                 unit_t_eta(i, j, k, tx, ty, tz);
    //                 push(tx, ty, tz, Jeta(i, j, k, 0), dl_eta(i, j, k, 0));

    //                 unit_t_eta(i + 1, j, k, tx, ty, tz);
    //                 push(tx, ty, tz, Jeta(i + 1, j, k, 0), dl_eta(i + 1, j, k, 0));

    //                 unit_t_eta(i, j, k + 1, tx, ty, tz);
    //                 push(tx, ty, tz, Jeta(i, j, k + 1, 0), dl_eta(i, j, k + 1, 0));

    //                 unit_t_eta(i + 1, j, k + 1, tx, ty, tz);
    //                 push(tx, ty, tz, Jeta(i + 1, j, k + 1, 0), dl_eta(i + 1, j, k + 1, 0));

    //                 // =====================================================
    //                 // 4 zeta-edges
    //                 // =====================================================
    //                 unit_t_zeta(i, j, k, tx, ty, tz);
    //                 push(tx, ty, tz, Jzeta(i, j, k, 0), dl_zeta(i, j, k, 0));

    //                 unit_t_zeta(i + 1, j, k, tx, ty, tz);
    //                 push(tx, ty, tz, Jzeta(i + 1, j, k, 0), dl_zeta(i + 1, j, k, 0));

    //                 unit_t_zeta(i, j + 1, k, tx, ty, tz);
    //                 push(tx, ty, tz, Jzeta(i, j + 1, k, 0), dl_zeta(i, j + 1, k, 0));

    //                 unit_t_zeta(i + 1, j + 1, k, tx, ty, tz);
    //                 push(tx, ty, tz, Jzeta(i + 1, j + 1, k, 0), dl_zeta(i + 1, j + 1, k, 0));

    //                 // =====================================================
    //                 // Weighted least squares:
    //                 //   minimize sum w | t·Jcell - J_edge/|dl| |^2
    //                 // =====================================================
    //                 double N00 = 0.0, N01 = 0.0, N02 = 0.0;
    //                 double N11 = 0.0, N12 = 0.0, N22 = 0.0;
    //                 double r0 = 0.0, r1 = 0.0, r2 = 0.0;

    //                 for (int n = 0; n < K; ++n)
    //                 {
    //                     const double w = eqs[n].w;
    //                     const double tx = eqs[n].tx;
    //                     const double ty = eqs[n].ty;
    //                     const double tz = eqs[n].tz;
    //                     const double b = eqs[n].rhs;

    //                     N00 += w * tx * tx;
    //                     N01 += w * tx * ty;
    //                     N02 += w * tx * tz;
    //                     N11 += w * ty * ty;
    //                     N12 += w * ty * tz;
    //                     N22 += w * tz * tz;

    //                     r0 += w * tx * b;
    //                     r1 += w * ty * b;
    //                     r2 += w * tz * b;
    //                 }

    //                 auto det3 = [&](double a, double b, double c,
    //                                 double d, double e, double f) -> double
    //                 {
    //                     // | a b c |
    //                     // | b d e |
    //                     // | c e f |
    //                     return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
    //                 };

    //                 double det = det3(N00, N01, N02, N11, N12, N22);
    //                 const double reg = 1e-14 * (N00 + N11 + N22 + 1.0);

    //                 if (std::abs(det) < reg)
    //                 {
    //                     N00 += reg;
    //                     N11 += reg;
    //                     N22 += reg;
    //                     det = det3(N00, N01, N02, N11, N12, N22);
    //                 }

    //                 // inverse of symmetric 3x3 normal matrix
    //                 const double C00 = (N11 * N22 - N12 * N12);
    //                 const double C01 = (N02 * N12 - N01 * N22);
    //                 const double C02 = (N01 * N12 - N02 * N11);
    //                 const double C11 = (N00 * N22 - N02 * N02);
    //                 const double C12 = (N01 * N02 - N00 * N12);
    //                 const double C22 = (N00 * N11 - N01 * N01);

    //                 const double invdet = 1.0 / det;

    //                 const double Jx =
    //                     invdet * (C00 * r0 + C01 * r1 + C02 * r2);
    //                 const double Jy =
    //                     invdet * (C01 * r0 + C11 * r1 + C12 * r2);
    //                 const double Jz =
    //                     invdet * (C02 * r0 + C12 * r1 + C22 * r2);

    //                 Jcell(i, j, k, 0) = Jx;
    //                 Jcell(i, j, k, 1) = Jy;
    //                 Jcell(i, j, k, 2) = Jz;
    //             }
    // }
    // mercury_bound_.Sync("J_cell");
}

void MercurySolver::calc_Jcell_from_Bcell_metric_()
{
    const int nb = fld_->num_blocks();

    constexpr double tiny = 1.0e-300;

    auto dd = [](double a, double b, double c,
                 double fp_i, double fm_i,
                 double fp_j, double fm_j,
                 double fp_k, double fm_k) -> double
    {
        return 0.5 * (a * (fp_i - fm_i) +
                      b * (fp_j - fm_j) +
                      c * (fp_k - fm_k));
    };

    auto dd_axis = [](int axis,
                      double a, double b, double c,
                      double fp_i, double fm_i,
                      double fp_j, double fm_j,
                      double fp_k, double fm_k) -> double
    {
        if (axis == 0)
            return 0.5 * a * (fp_i - fm_i);
        if (axis == 1)
            return 0.5 * b * (fp_j - fm_j);
        return 0.5 * c * (fp_k - fm_k);
    };

    auto get_comp = [](const Int3 &a, int ax) -> int
    {
        if (ax == 0)
            return a.i;
        if (ax == 1)
            return a.j;
        return a.k;
    };

    auto set_comp = [](Int3 &a, int ax, int v)
    {
        if (ax == 0)
            a.i = v;
        else if (ax == 1)
            a.j = v;
        else
            a.k = v;
    };

    auto derv_at = [&](FieldBlock &Jac,
                       FieldBlock &Axi,
                       FieldBlock &Aet,
                       FieldBlock &Aze,
                       int i, int j, int k,
                       double &ax, double &ay, double &az,
                       double &bx, double &by, double &bz,
                       double &cx, double &cy, double &cz)
    {
        const double V = std::abs(Jac(i, j, k, 0));

        if (V <= tiny)
        {
            ax = ay = az = 0.0;
            bx = by = bz = 0.0;
            cx = cy = cz = 0.0;
            return;
        }

        ax = 0.5 * (Axi(i + 1, j, k, 0) + Axi(i, j, k, 0)) / V;
        ay = 0.5 * (Axi(i + 1, j, k, 1) + Axi(i, j, k, 1)) / V;
        az = 0.5 * (Axi(i + 1, j, k, 2) + Axi(i, j, k, 2)) / V;

        bx = 0.5 * (Aet(i, j + 1, k, 0) + Aet(i, j, k, 0)) / V;
        by = 0.5 * (Aet(i, j + 1, k, 1) + Aet(i, j, k, 1)) / V;
        bz = 0.5 * (Aet(i, j + 1, k, 2) + Aet(i, j, k, 2)) / V;

        cx = 0.5 * (Aze(i, j, k + 1, 0) + Aze(i, j, k, 0)) / V;
        cy = 0.5 * (Aze(i, j, k + 1, 1) + Aze(i, j, k, 1)) / V;
        cz = 0.5 * (Aze(i, j, k + 1, 2) + Aze(i, j, k, 2)) / V;
    };

    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Bcell = fld_->field(fid_.fid_Bindcell, ib);
        FieldBlock &Jcell = fld_->field(fid_.fid_Jcell, ib);

        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);
        FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);
        FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib);

        if (!Bcell.is_allocated() || !Jcell.is_allocated())
            continue;

        if (!Jac.is_allocated() || !Axi.is_allocated() ||
            !Aet.is_allocated() || !Aze.is_allocated())
            continue;

        Int3 lo = Jcell.inner_lo();
        Int3 hi = Jcell.inner_hi();

        // 直接计算整个 inner 区域。
        // i±1, j±1, k±1 全部依赖已经处理好的 ghost / halo。
        for (int i = lo.i; i < hi.i; ++i)
        {
            for (int j = lo.j; j < hi.j; ++j)
            {
                for (int k = lo.k; k < hi.k; ++k)
                {
                    double ax, ay, az;
                    double bx, by, bz;
                    double cx, cy, cz;

                    derv_at(Jac, Axi, Aet, Aze,
                            i, j, k,
                            ax, ay, az,
                            bx, by, bz,
                            cx, cy, cz);

                    const double dBx_dx = dd(ax, bx, cx,
                                             Bcell(i + 1, j, k, 0), Bcell(i - 1, j, k, 0),
                                             Bcell(i, j + 1, k, 0), Bcell(i, j - 1, k, 0),
                                             Bcell(i, j, k + 1, 0), Bcell(i, j, k - 1, 0));

                    const double dBx_dy = dd(ay, by, cy,
                                             Bcell(i + 1, j, k, 0), Bcell(i - 1, j, k, 0),
                                             Bcell(i, j + 1, k, 0), Bcell(i, j - 1, k, 0),
                                             Bcell(i, j, k + 1, 0), Bcell(i, j, k - 1, 0));

                    const double dBx_dz = dd(az, bz, cz,
                                             Bcell(i + 1, j, k, 0), Bcell(i - 1, j, k, 0),
                                             Bcell(i, j + 1, k, 0), Bcell(i, j - 1, k, 0),
                                             Bcell(i, j, k + 1, 0), Bcell(i, j, k - 1, 0));

                    const double dBy_dx = dd(ax, bx, cx,
                                             Bcell(i + 1, j, k, 1), Bcell(i - 1, j, k, 1),
                                             Bcell(i, j + 1, k, 1), Bcell(i, j - 1, k, 1),
                                             Bcell(i, j, k + 1, 1), Bcell(i, j, k - 1, 1));

                    const double dBy_dy = dd(ay, by, cy,
                                             Bcell(i + 1, j, k, 1), Bcell(i - 1, j, k, 1),
                                             Bcell(i, j + 1, k, 1), Bcell(i, j - 1, k, 1),
                                             Bcell(i, j, k + 1, 1), Bcell(i, j, k - 1, 1));

                    const double dBy_dz = dd(az, bz, cz,
                                             Bcell(i + 1, j, k, 1), Bcell(i - 1, j, k, 1),
                                             Bcell(i, j + 1, k, 1), Bcell(i, j - 1, k, 1),
                                             Bcell(i, j, k + 1, 1), Bcell(i, j, k - 1, 1));

                    const double dBz_dx = dd(ax, bx, cx,
                                             Bcell(i + 1, j, k, 2), Bcell(i - 1, j, k, 2),
                                             Bcell(i, j + 1, k, 2), Bcell(i, j - 1, k, 2),
                                             Bcell(i, j, k + 1, 2), Bcell(i, j, k - 1, 2));

                    const double dBz_dy = dd(ay, by, cy,
                                             Bcell(i + 1, j, k, 2), Bcell(i - 1, j, k, 2),
                                             Bcell(i, j + 1, k, 2), Bcell(i, j - 1, k, 2),
                                             Bcell(i, j, k + 1, 2), Bcell(i, j, k - 1, 2));

                    const double dBz_dz = dd(az, bz, cz,
                                             Bcell(i + 1, j, k, 2), Bcell(i - 1, j, k, 2),
                                             Bcell(i, j + 1, k, 2), Bcell(i, j - 1, k, 2),
                                             Bcell(i, j, k + 1, 2), Bcell(i, j, k - 1, 2));

                    Jcell(i, j, k, 0) = dBz_dy - dBy_dz;
                    Jcell(i, j, k, 1) = dBx_dz - dBz_dx;
                    Jcell(i, j, k, 2) = dBy_dx - dBx_dy;
                }
            }
        }

        // Pole cells are degenerate in the normal direction and collapsed
        // around zeta/k. Recompute those cells with only the axial
        // computational derivative: eta for xi-normal Poles, xi for eta-normal Poles.
        for (const auto &p : topo_->physical_patches)
        {
            if (p.this_block != ib)
                continue;
            if (p.bc_name != "Pole")
                continue;

            const int dir = std::abs(p.direction);
            if (dir != 1 && dir != 2)
                continue;

            const int norm_axis = dir - 1;
            const int axial_axis = (dir == 1) ? 1 : 0;
            const bool high_side = (p.direction > 0);

            Int3 plo = lo;
            Int3 phi = hi;

            for (int d = 0; d < 3; ++d)
            {
                if (d == norm_axis)
                    continue;

                const int tlo = std::max(get_comp(lo, d), get_comp(p.this_box_node.lo, d));
                const int thi = std::min(get_comp(hi, d), get_comp(p.this_box_node.hi, d) - 1);
                set_comp(plo, d, tlo);
                set_comp(phi, d, thi);
            }

            if (high_side)
            {
                set_comp(plo, norm_axis, get_comp(hi, norm_axis) - 1);
                set_comp(phi, norm_axis, get_comp(hi, norm_axis));
            }
            else
            {
                set_comp(plo, norm_axis, get_comp(lo, norm_axis));
                set_comp(phi, norm_axis, get_comp(lo, norm_axis) + 1);
            }

            if (!(plo.i < phi.i && plo.j < phi.j && plo.k < phi.k))
                continue;

            auto calc_curl_at = [&](int i, int j, int k,
                                    double ax, double ay, double az,
                                    double bx, double by, double bz,
                                    double cx, double cy, double cz,
                                    double &Jx, double &Jy, double &Jz)
            {
                const double dBx_dx = dd_axis(axial_axis, ax, bx, cx,
                                              Bcell(i + 1, j, k, 0), Bcell(i - 1, j, k, 0),
                                              Bcell(i, j + 1, k, 0), Bcell(i, j - 1, k, 0),
                                              Bcell(i, j, k + 1, 0), Bcell(i, j, k - 1, 0));

                const double dBx_dy = dd_axis(axial_axis, ay, by, cy,
                                              Bcell(i + 1, j, k, 0), Bcell(i - 1, j, k, 0),
                                              Bcell(i, j + 1, k, 0), Bcell(i, j - 1, k, 0),
                                              Bcell(i, j, k + 1, 0), Bcell(i, j, k - 1, 0));

                const double dBx_dz = dd_axis(axial_axis, az, bz, cz,
                                              Bcell(i + 1, j, k, 0), Bcell(i - 1, j, k, 0),
                                              Bcell(i, j + 1, k, 0), Bcell(i, j - 1, k, 0),
                                              Bcell(i, j, k + 1, 0), Bcell(i, j, k - 1, 0));

                const double dBy_dx = dd_axis(axial_axis, ax, bx, cx,
                                              Bcell(i + 1, j, k, 1), Bcell(i - 1, j, k, 1),
                                              Bcell(i, j + 1, k, 1), Bcell(i, j - 1, k, 1),
                                              Bcell(i, j, k + 1, 1), Bcell(i, j, k - 1, 1));

                const double dBy_dy = dd_axis(axial_axis, ay, by, cy,
                                              Bcell(i + 1, j, k, 1), Bcell(i - 1, j, k, 1),
                                              Bcell(i, j + 1, k, 1), Bcell(i, j - 1, k, 1),
                                              Bcell(i, j, k + 1, 1), Bcell(i, j, k - 1, 1));

                const double dBy_dz = dd_axis(axial_axis, az, bz, cz,
                                              Bcell(i + 1, j, k, 1), Bcell(i - 1, j, k, 1),
                                              Bcell(i, j + 1, k, 1), Bcell(i, j - 1, k, 1),
                                              Bcell(i, j, k + 1, 1), Bcell(i, j, k - 1, 1));

                const double dBz_dx = dd_axis(axial_axis, ax, bx, cx,
                                              Bcell(i + 1, j, k, 2), Bcell(i - 1, j, k, 2),
                                              Bcell(i, j + 1, k, 2), Bcell(i, j - 1, k, 2),
                                              Bcell(i, j, k + 1, 2), Bcell(i, j, k - 1, 2));

                const double dBz_dy = dd_axis(axial_axis, ay, by, cy,
                                              Bcell(i + 1, j, k, 2), Bcell(i - 1, j, k, 2),
                                              Bcell(i, j + 1, k, 2), Bcell(i, j - 1, k, 2),
                                              Bcell(i, j, k + 1, 2), Bcell(i, j, k - 1, 2));

                const double dBz_dz = dd_axis(axial_axis, az, bz, cz,
                                              Bcell(i + 1, j, k, 2), Bcell(i - 1, j, k, 2),
                                              Bcell(i, j + 1, k, 2), Bcell(i, j - 1, k, 2),
                                              Bcell(i, j, k + 1, 2), Bcell(i, j, k - 1, 2));

                Jx = dBz_dy - dBy_dz;
                Jy = dBx_dz - dBz_dx;
                Jz = dBy_dx - dBx_dy;
            };

            const int k_count = phi.k - plo.k;
            const double inv_k_count = 1.0 / static_cast<double>(k_count);

            for (int i = plo.i; i < phi.i; ++i)
            {
                for (int j = plo.j; j < phi.j; ++j)
                {
                    double ax_avg = 0.0, ay_avg = 0.0, az_avg = 0.0;
                    double bx_avg = 0.0, by_avg = 0.0, bz_avg = 0.0;
                    double cx_avg = 0.0, cy_avg = 0.0, cz_avg = 0.0;

                    for (int k = plo.k; k < phi.k; ++k)
                    {
                        double ax, ay, az;
                        double bx, by, bz;
                        double cx, cy, cz;

                        derv_at(Jac, Axi, Aet, Aze,
                                i, j, k,
                                ax, ay, az,
                                bx, by, bz,
                                cx, cy, cz);

                        if (axial_axis == 0)
                        {
                            ax_avg += ax;
                            ay_avg += ay;
                            az_avg += az;
                        }
                        else
                        {
                            bx_avg += bx;
                            by_avg += by;
                            bz_avg += bz;
                        }
                    }

                    if (axial_axis == 0)
                    {
                        ax_avg *= inv_k_count;
                        ay_avg *= inv_k_count;
                        az_avg *= inv_k_count;
                    }
                    else
                    {
                        bx_avg *= inv_k_count;
                        by_avg *= inv_k_count;
                        bz_avg *= inv_k_count;
                    }

                    double Jx, Jy, Jz;
                    calc_curl_at(i, j, plo.k,
                                 ax_avg, ay_avg, az_avg,
                                 bx_avg, by_avg, bz_avg,
                                 cx_avg, cy_avg, cz_avg,
                                 Jx, Jy, Jz);

                    for (int k = plo.k; k < phi.k; ++k)
                    {
                        Jcell(i, j, k, 0) = Jx;
                        Jcell(i, j, k, 1) = Jy;
                        Jcell(i, j, k, 2) = Jz;
                    }
                }
            }
        }

        // Coupled-Solid wall cells are zeroed by the boundary handler below.
        // For the first fluid layer next to that wall layer, avoid the wall
        // Bcell in the normal derivative to suppress numerical wall-current
        // contamination.
        for (const auto &p : topo_->physical_patches)
        {
            if (p.this_block != ib)
                continue;
            if (p.bc_name != "Coupled-Solid")
                continue;

            const int dir = std::abs(p.direction);
            if (dir < 1 || dir > 3)
                continue;

            const int norm_axis = dir - 1;
            const bool high_side = (p.direction > 0);
            const int inward = high_side ? -1 : +1;
            const int nlo = get_comp(lo, norm_axis);
            const int nhi = get_comp(hi, norm_axis);

            if (nhi - nlo < 3)
                continue;

            Int3 plo = lo;
            Int3 phi = hi;

            for (int d = 0; d < 3; ++d)
            {
                if (d == norm_axis)
                    continue;

                const int tlo = std::max(get_comp(lo, d), get_comp(p.this_box_node.lo, d));
                const int thi = std::min(get_comp(hi, d), get_comp(p.this_box_node.hi, d) - 1);
                set_comp(plo, d, tlo);
                set_comp(phi, d, thi);
            }

            if (high_side)
            {
                set_comp(plo, norm_axis, nhi - 2);
                set_comp(phi, norm_axis, nhi - 1);
            }
            else
            {
                set_comp(plo, norm_axis, nlo + 1);
                set_comp(phi, norm_axis, nlo + 2);
            }

            if (!(plo.i < phi.i && plo.j < phi.j && plo.k < phi.k))
                continue;

            auto B_at_axis_offset = [&](int i, int j, int k, int axis, int offset, int comp) -> double
            {
                if (axis == 0)
                    return Bcell(i + offset, j, k, comp);
                if (axis == 1)
                    return Bcell(i, j + offset, k, comp);
                return Bcell(i, j, k + offset, comp);
            };

            auto dcomp_axis = [&](int i, int j, int k, int axis, int comp) -> double
            {
                if (axis == norm_axis)
                {
                    if (inward > 0)
                        return B_at_axis_offset(i, j, k, axis, +1, comp) -
                               B_at_axis_offset(i, j, k, axis, 0, comp);

                    return B_at_axis_offset(i, j, k, axis, 0, comp) -
                           B_at_axis_offset(i, j, k, axis, -1, comp);
                }

                return 0.5 * (B_at_axis_offset(i, j, k, axis, +1, comp) -
                              B_at_axis_offset(i, j, k, axis, -1, comp));
            };

            auto dphys = [&](int i, int j, int k, int comp,
                             double qxi, double qet, double qze) -> double
            {
                return qxi * dcomp_axis(i, j, k, 0, comp) +
                       qet * dcomp_axis(i, j, k, 1, comp) +
                       qze * dcomp_axis(i, j, k, 2, comp);
            };

            for (int i = plo.i; i < phi.i; ++i)
            {
                for (int j = plo.j; j < phi.j; ++j)
                {
                    for (int k = plo.k; k < phi.k; ++k)
                    {
                        double ax, ay, az;
                        double bx, by, bz;
                        double cx, cy, cz;

                        derv_at(Jac, Axi, Aet, Aze,
                                i, j, k,
                                ax, ay, az,
                                bx, by, bz,
                                cx, cy, cz);

                        const double dBx_dx = dphys(i, j, k, 0, ax, bx, cx);
                        const double dBx_dy = dphys(i, j, k, 0, ay, by, cy);
                        const double dBx_dz = dphys(i, j, k, 0, az, bz, cz);

                        const double dBy_dx = dphys(i, j, k, 1, ax, bx, cx);
                        const double dBy_dy = dphys(i, j, k, 1, ay, by, cy);
                        const double dBy_dz = dphys(i, j, k, 1, az, bz, cz);

                        const double dBz_dx = dphys(i, j, k, 2, ax, bx, cx);
                        const double dBz_dy = dphys(i, j, k, 2, ay, by, cy);
                        const double dBz_dz = dphys(i, j, k, 2, az, bz, cz);

                        Jcell(i, j, k, 0) = dBz_dy - dBy_dz;
                        Jcell(i, j, k, 1) = dBx_dz - dBz_dx;
                        Jcell(i, j, k, 2) = dBy_dx - dBx_dy;
                    }
                }
            }
        }

        // Intersection of a Pole collapse and the first fluid layer next to a
        // Coupled-Solid wall. This must run after both treatments above: keep
        // the Pole k-collapse, but use a one-sided axial derivative away from
        // the wall so the wall-layer Bcell is not sampled.
        for (const auto &ppole : topo_->physical_patches)
        {
            if (ppole.this_block != ib)
                continue;
            if (ppole.bc_name != "Pole")
                continue;

            const int pole_dir = std::abs(ppole.direction);
            if (pole_dir != 1 && pole_dir != 2)
                continue;

            const int pole_norm_axis = pole_dir - 1;
            const int axial_axis = (pole_dir == 1) ? 1 : 0;
            const bool pole_high_side = (ppole.direction > 0);

            for (const auto &pwall : topo_->physical_patches)
            {
                if (pwall.this_block != ib)
                    continue;
                if (pwall.bc_name != "Coupled-Solid")
                    continue;

                const int wall_dir = std::abs(pwall.direction);
                if (wall_dir - 1 != axial_axis)
                    continue;

                const bool wall_high_side = (pwall.direction > 0);
                const int wall_inward = wall_high_side ? -1 : +1;
                const int axial_lo = get_comp(lo, axial_axis);
                const int axial_hi = get_comp(hi, axial_axis);

                if (axial_hi - axial_lo < 3)
                    continue;

                Int3 plo = lo;
                Int3 phi = hi;

                if (pole_high_side)
                {
                    set_comp(plo, pole_norm_axis, get_comp(hi, pole_norm_axis) - 1);
                    set_comp(phi, pole_norm_axis, get_comp(hi, pole_norm_axis));
                }
                else
                {
                    set_comp(plo, pole_norm_axis, get_comp(lo, pole_norm_axis));
                    set_comp(phi, pole_norm_axis, get_comp(lo, pole_norm_axis) + 1);
                }

                if (wall_high_side)
                {
                    set_comp(plo, axial_axis, axial_hi - 2);
                    set_comp(phi, axial_axis, axial_hi - 1);
                }
                else
                {
                    set_comp(plo, axial_axis, axial_lo + 1);
                    set_comp(phi, axial_axis, axial_lo + 2);
                }

                const int klo = std::max({lo.k, ppole.this_box_node.lo.k, pwall.this_box_node.lo.k});
                const int khi = std::min({hi.k, ppole.this_box_node.hi.k - 1, pwall.this_box_node.hi.k - 1});
                plo.k = klo;
                phi.k = khi;

                if (!(plo.i < phi.i && plo.j < phi.j && plo.k < phi.k))
                    continue;

                auto B_axis = [&](int i, int j, int k, int offset, int comp) -> double
                {
                    if (axial_axis == 0)
                        return Bcell(i + offset, j, k, comp);
                    return Bcell(i, j + offset, k, comp);
                };

                const int k_count = phi.k - plo.k;
                const double inv_k_count = 1.0 / static_cast<double>(k_count);

                for (int i = plo.i; i < phi.i; ++i)
                {
                    for (int j = plo.j; j < phi.j; ++j)
                    {
                        double ax_avg = 0.0, ay_avg = 0.0, az_avg = 0.0;
                        double bx_avg = 0.0, by_avg = 0.0, bz_avg = 0.0;

                        for (int k = plo.k; k < phi.k; ++k)
                        {
                            double ax, ay, az;
                            double bx, by, bz;
                            double cx, cy, cz;

                            derv_at(Jac, Axi, Aet, Aze,
                                    i, j, k,
                                    ax, ay, az,
                                    bx, by, bz,
                                    cx, cy, cz);

                            if (axial_axis == 0)
                            {
                                ax_avg += ax;
                                ay_avg += ay;
                                az_avg += az;
                            }
                            else
                            {
                                bx_avg += bx;
                                by_avg += by;
                                bz_avg += bz;
                            }
                        }

                        ax_avg *= inv_k_count;
                        ay_avg *= inv_k_count;
                        az_avg *= inv_k_count;
                        bx_avg *= inv_k_count;
                        by_avg *= inv_k_count;
                        bz_avg *= inv_k_count;

                        const int k0 = plo.k;
                        double dB_daxis[3] = {0.0, 0.0, 0.0};
                        for (int comp = 0; comp < 3; ++comp)
                        {
                            if (wall_inward > 0)
                                dB_daxis[comp] = B_axis(i, j, k0, +1, comp) -
                                                 B_axis(i, j, k0, 0, comp);
                            else
                                dB_daxis[comp] = B_axis(i, j, k0, 0, comp) -
                                                 B_axis(i, j, k0, -1, comp);
                        }

                        const double qx = (axial_axis == 0) ? ax_avg : bx_avg;
                        const double qy = (axial_axis == 0) ? ay_avg : by_avg;
                        const double qz = (axial_axis == 0) ? az_avg : bz_avg;

                        const double dBx_dx = qx * dB_daxis[0];
                        const double dBx_dy = qy * dB_daxis[0];
                        const double dBx_dz = qz * dB_daxis[0];

                        const double dBy_dx = qx * dB_daxis[1];
                        const double dBy_dy = qy * dB_daxis[1];
                        const double dBy_dz = qz * dB_daxis[1];

                        const double dBz_dx = qx * dB_daxis[2];
                        const double dBz_dy = qy * dB_daxis[2];
                        const double dBz_dz = qz * dB_daxis[2];

                        const double Jx = dBz_dy - dBy_dz;
                        const double Jy = dBx_dz - dBz_dx;
                        const double Jz = dBy_dx - dBx_dy;

                        for (int k = plo.k; k < phi.k; ++k)
                        {
                            Jcell(i, j, k, 0) = Jx;
                            Jcell(i, j, k, 1) = Jy;
                            Jcell(i, j, k, 2) = Jz;
                        }
                    }
                }
            }
        }
    }

    // 后续所有 Pole / wall / coupling / halo 处理都交给边界系统。
    mercury_bound_.Sync("J_cell");
}
