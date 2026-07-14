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

double MercurySolver::MercuryResistivityShape_(double radius) const
{
    const double r0 = resist_control.radial_inner;
    const double r1 = resist_control.radial_outer;
    const double width = resist_control.radial_width;

    // Do not truncate at r0/r1: those are the tanh transition centers, where
    // the window is approximately 1/2.  The previous hard cutoff therefore
    // introduced two O(eta_max) coefficient jumps and circular current sheets.
    const double shape = 0.5 *
        (std::tanh((radius - r0) / width) -
         std::tanh((radius - r1) / width));

    // Keep the implicit unknown set compact, but truncate only in the
    // numerically negligible tail rather than at the transition center.
    constexpr double tail_floor = 1.0e-10;
    if (!std::isfinite(shape) || shape <= tail_floor)
        return 0.0;
    return std::min(1.0, shape);
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

    // Uplus is reconstructed over the complete allocated range, including
    // ghosts, from U_H/U_Na and PV after Ucell synchronization.  It exists
    // only on Fluid blocks, so a topology-wide component-copy halo would
    // incorrectly try to pack inactive Solid-block storage.
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
