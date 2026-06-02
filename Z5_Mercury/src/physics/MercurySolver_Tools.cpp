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

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UN = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &PVH = fld_->field(fid_.fid_PV_H, ib);
        FieldBlock &PVN = fld_->field(fid_.fid_PV_Na, ib);
        FieldBlock &Up = fld_->field(fid_.fid_U_plus, ib);

        if (!Up.is_allocated())
            continue;

        const Int3 lo = Up.get_lo();
        const Int3 hi = Up.get_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    Up(i, j, k, 0) = 0.0;
                    Up(i, j, k, 1) = 0.0;
                    Up(i, j, k, 2) = 0.0;
                }

        if (!UH.is_allocated() || !UN.is_allocated() ||
            !PVH.is_allocated() || !PVN.is_allocated())
            continue;

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double rhoH0 = std::max(UH(i, j, k, 0), 0.0);
                    const double rhoNa0 = std::max(UN(i, j, k, 0), 0.0);

                    const double uH = (rhoH0 > rho_eps) ? PVH(i, j, k, 0) : 0.0;
                    const double vH = (rhoH0 > rho_eps) ? PVH(i, j, k, 1) : 0.0;
                    const double wH = (rhoH0 > rho_eps) ? PVH(i, j, k, 2) : 0.0;

                    const double uNa = (rhoNa0 > rho_eps) ? PVN(i, j, k, 0) : 0.0;
                    const double vNa = (rhoNa0 > rho_eps) ? PVN(i, j, k, 1) : 0.0;
                    const double wNa = (rhoNa0 > rho_eps) ? PVN(i, j, k, 2) : 0.0;

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

    mercury_bound_.Sync("Uplus");
}

void MercurySolver::UpdateFluidDerivedFields_()
{
    calc_PV();
    calc_Uplus();
}

void MercurySolver::UpdateMagneticDerivedFields_()
{
    // B_face is the CT state; J_cell is reconstructed from mimetic J_edge.
    calc_Bcell();
    Calc_J_Edge();
    calc_Jcell();
}

void MercurySolver::UpdateDerivedFields_()
{
    UpdateMagneticDerivedFields_();
    UpdateFluidDerivedFields_();
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

        auto &W = fld_->field(fid_.fid_Bcell_from_Bface_w, ib);

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

    if (topo_)
    {
        for (const auto &p : topo_->physical_patches)
        {
            if (p.bc_name != "Pole")
                continue;

            const int dir = std::abs(p.direction);
            if (dir != 1 && dir != 2)
                continue;

            const int ib = p.this_block;
            if (ib < 0 || ib >= nblock)
                continue;

            auto &Bindcell = fld_->field(fid_.fid_Bindcell, ib);
            auto &Bzeta = fld_->field(fid_.fid_B.zeta, ib);
            // auto &Baddzeta = fld_->field(fid_.fid_Badd.zeta, ib);
            auto &Aze = fld_->field(fid_.fid_metric.zeta, ib);

            if (!Bzeta.is_allocated() || !Aze.is_allocated())
                continue;

            const int sgn = (p.direction > 0) ? +1 : -1;

            auto write_zeta_flux = [&](int i, int j, int k)
            {
                const int kc = k;
                const double phi =
                    Bindcell(i, j, kc, 0) * Aze(i, j, k, 0) +
                    Bindcell(i, j, kc, 1) * Aze(i, j, k, 1) +
                    Bindcell(i, j, kc, 2) * Aze(i, j, k, 2);

                Bzeta(i, j, k, 0) = phi; //- Baddzeta(i, j, k, 0);
            };

            const Box3 &node = p.this_box_node;

            if (dir == 1)
            {
                const int icell = (sgn < 0) ? node.lo.i : (node.lo.i - 1);

                const int jface_lo = node.lo.j;
                const int jface_hi = node.hi.j - 1;
                const int kface_lo = node.lo.k;
                const int kface_hi = node.hi.k;

                for (int j = jface_lo; j < jface_hi; ++j)
                    for (int k = kface_lo; k < kface_hi; ++k)
                        write_zeta_flux(icell, j, k);
            }
            else
            {
                const int jcell = (sgn < 0) ? node.lo.j : (node.lo.j - 1);

                const int iface_lo = node.lo.i;
                const int iface_hi = node.hi.i - 1;
                const int kface_lo = node.lo.k;
                const int kface_hi = node.hi.k;

                for (int i = iface_lo; i < iface_hi; ++i)
                    for (int k = kface_lo; k < kface_hi; ++k)
                        write_zeta_flux(i, jcell, k);
            }
        }
    }

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

void MercurySolver::Calc_J_Edge()
{
    //  ComputeJ_AtEdges_Inner_();
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Bxi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_J.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_J.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, iblk);

        auto &beta_xi = fld_->field(fid_.Face_beta.xi, iblk); // Hodge *: 2-form face -> 1-form face
        auto &beta_eta = fld_->field(fid_.Face_beta.eta, iblk);
        auto &beta_zeta = fld_->field(fid_.Face_beta.zeta, iblk);

        auto &alpha_xi = fld_->field(fid_.Edge_alpha.xi, iblk); // Hodge *: 2-form edge -> 1-form edge
        auto &alpha_eta = fld_->field(fid_.Edge_alpha.eta, iblk);
        auto &alpha_zeta = fld_->field(fid_.Edge_alpha.zeta, iblk);

        // compute J (edge 1-form) from face B (2-form)
        // multiper 用 +1.0 J =curl B。
        CTOperators::CurlAdjFaceToEdge(iblk,
                                       Bxi, Beta, Bzeta,
                                       beta_xi, beta_eta, beta_zeta,
                                       Jxi, Jeta, Jzeta,
                                       /*multiper=*/1.0);

        //  J_edge^(1-form) = (alpha_edge / mu0) * Jcirc_edge
        //     alpha = |e|/|S*|  (⋆1^{-1})
        {
            // Edge_xi
            Int3 lo = Jxi.inner_lo(), hi = Jxi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jxi(i, j, k, 0) *= (alpha_xi(i, j, k, 0));
        }
        {
            // Edge_eta
            Int3 lo = Jeta.inner_lo(), hi = Jeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jeta(i, j, k, 0) *= (alpha_eta(i, j, k, 0));
        }
        {
            // Edge_zeta
            Int3 lo = Jzeta.inner_lo(), hi = Jzeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jzeta(i, j, k, 0) *= (alpha_zeta(i, j, k, 0));
        }
    }

    mercury_bound_.Sync("Jedge"); // ApplyBC_EdgeJ_();
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

        auto &W = fld_->field(fid_.fid_Jcell_from_Jedge_w, ib);

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

    if (topo_)
    {
        auto solve3x3 = [](double A[3][3], double b[3], double x[3]) -> bool
        {
            const double trA = A[0][0] + A[1][1] + A[2][2];
            const double reg = 1.0e-12 * std::max(1.0, trA);

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

        constexpr double eps_len = 1.0e-14;

        for (const auto &p : topo_->physical_patches)
        {
            if (p.bc_name != "Pole")
                continue;

            const int dir = std::abs(p.direction);
            if (dir != 1 && dir != 2)
                continue;

            const int ib = p.this_block;
            if (ib < 0 || ib >= nblock)
                continue;

            auto &Jcell = fld_->field(fid_.fid_Jcell, ib);
            auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
            auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
            auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);
            auto &dr_xi = fld_->field(fid_.Edge_dr.xi, ib);
            auto &dr_eta = fld_->field(fid_.Edge_dr.eta, ib);
            auto &dr_zeta = fld_->field(fid_.Edge_dr.zeta, ib);

            if (!Jcell.is_allocated() || !Jxi.is_allocated() ||
                !Jeta.is_allocated() || !Jzeta.is_allocated() ||
                !dr_xi.is_allocated() || !dr_eta.is_allocated() ||
                !dr_zeta.is_allocated())
                continue;

            const bool high_side = (p.direction > 0);
            const Box3 &node = p.this_box_node;

            auto push_edge = [&](FieldBlock &Je, FieldBlock &dr,
                                 int i, int j, int k,
                                 double A[3][3], double b[3], int &cnt)
            {
                const double dx = dr(i, j, k, 0);
                const double dy = dr(i, j, k, 1);
                const double dz = dr(i, j, k, 2);
                const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

                if (L < eps_len)
                    return;

                const double tau[3] = {dx / L, dy / L, dz / L};
                const double y = Je(i, j, k, 0) / L;

                for (int a = 0; a < 3; ++a)
                {
                    b[a] += y * tau[a];
                    for (int c = 0; c < 3; ++c)
                        A[a][c] += tau[a] * tau[c];
                }

                ++cnt;
            };

            auto write_ring = [&](int ic, int jc,
                                  double A[3][3], double b[3], int cnt,
                                  int klo, int khi)
            {
                if (cnt <= 0)
                    return;

                double Jp[3] = {0.0, 0.0, 0.0};
                if (!solve3x3(A, b, Jp))
                    return;

                for (int k = klo; k < khi; ++k)
                {
                    Jcell(ic, jc, k, 0) = Jp[0];
                    Jcell(ic, jc, k, 1) = Jp[1];
                    Jcell(ic, jc, k, 2) = Jp[2];
                }
            };

            if (dir == 1)
            {
                const int ic = high_side ? (node.lo.i - 1) : node.lo.i;
                const int it = high_side ? ic : (ic + 1);

                for (int j = node.lo.j; j < node.hi.j - 1; ++j)
                {
                    double A[3][3] = {
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
                    double b[3] = {0.0, 0.0, 0.0};
                    int cnt = 0;

                    for (int k = node.lo.k; k < node.hi.k - 1; ++k)
                    {
                        // norm edges on the Pole layer
                        push_edge(Jxi, dr_xi, ic, j, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j + 1, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j, k + 1, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j + 1, k + 1, A, b, cnt);

                        // axis/rotate edges one layer away from the Pole
                        push_edge(Jeta, dr_eta, it, j, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, it, j, k + 1, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, it, j, k, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, it, j + 1, k, A, b, cnt);
                    }

                    write_ring(ic, j, A, b, cnt, node.lo.k, node.hi.k - 1);
                }
            }
            else
            {
                const int jc = high_side ? (node.lo.j - 1) : node.lo.j;
                const int jt = high_side ? jc : (jc + 1);

                for (int i = node.lo.i; i < node.hi.i - 1; ++i)
                {
                    double A[3][3] = {
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
                    double b[3] = {0.0, 0.0, 0.0};
                    int cnt = 0;

                    for (int k = node.lo.k; k < node.hi.k - 1; ++k)
                    {
                        // norm edges on the Pole layer
                        push_edge(Jeta, dr_eta, i, jc, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, i + 1, jc, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, i, jc, k + 1, A, b, cnt);
                        push_edge(Jeta, dr_eta, i + 1, jc, k + 1, A, b, cnt);

                        // axis/rotate edges one layer away from the Pole
                        push_edge(Jxi, dr_xi, i, jt, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, i, jt, k + 1, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, i, jt, k, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, i + 1, jt, k, A, b, cnt);
                    }

                    write_ring(i, jc, A, b, cnt, node.lo.k, node.hi.k - 1);
                }
            }
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
