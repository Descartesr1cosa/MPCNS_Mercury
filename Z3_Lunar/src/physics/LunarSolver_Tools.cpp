#include "LunarSolver.h"

void LunarSolver::calc_physical_constant(Param *par)
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
    m_H = M_H / NA;                     // kg/particle

    // 无量纲状态方程系数：p = rho * T * coeff
    // coeff = (R_uni * T_ref) / (M * U_ref^2)
    state_coeff_H = (R_uni * T_ref) / (M_H * U_ref * U_ref);

    CFL = par_->GetDou("CFL");

    hall_coef = B_ref * M_ref / (U_ref * q_e * L_ref * mu0 * rho_ref * NA);

    // ambi_coef = M_H * U_ref / (cst.data["q_e"] * ref.data["L_ref"] * ref.data["B_ref"]);

    momentum_induce_coeff = (q_e * L_ref * B_ref * n_ref) / (rho_ref * U_ref); // incude_coeff * n_ns * (u_ns - u_+) \times B  = momentum eqs source, n_ns: m^-3
    momentum_hall_coeff = (B_ref * B_ref) / (mu0 * rho_ref * U_ref * U_ref);   // momentum_hall_coeff * n_ns / n_e * \nabla\times B\times B = momentum eqs source

    inver_MA2 = ref.data["B_ref"] * ref.data["B_ref"] / (U_ref * U_ref * cst.data["mu_mag"] * rho_ref);

    inver_Rem = par->GetDou("eta_max_lunar") / (cst.data["mu_mag"] * U_ref * ref.data["L_ref"]);

    ne_hall_floor = 0.01;                              // nondimensional
    ne_hall_floor_dimensional = ne_hall_floor * n_ref; // dimension
}

NumInfo LunarSolver::Hall_Num_Limiter(double rhoH)
{
    NumInfo nd;
    const double rhoH_pos = std::max(rhoH, 0.0);
    nd.nH_true = rhoH_pos / M_H * M_ref;
    nd.ne_true = nd.nH_true;
    nd.ne_eff = std::sqrt(nd.ne_true * nd.ne_true + ne_hall_floor * ne_hall_floor);
    nd.chiH = 1.0;
    nd.wH_mhd = nd.nH_true / nd.ne_eff;
    nd.mhd_taper = nd.wH_mhd;
    return nd;
}

void LunarSolver::calc_PV()
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
}

void LunarSolver::calc_Uplus()
{
    const double rho_eps = 1e-20;

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &PVH = fld_->field(fid_.fid_PV_H, ib);
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

        if (!UH.is_allocated() || !PVH.is_allocated())
            continue;

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double rhoH0 = std::max(UH(i, j, k, 0), 0.0);
                    const double uH = (rhoH0 > rho_eps) ? PVH(i, j, k, 0) : 0.0;
                    const double vH = (rhoH0 > rho_eps) ? PVH(i, j, k, 1) : 0.0;
                    const double wH = (rhoH0 > rho_eps) ? PVH(i, j, k, 2) : 0.0;

                    Up(i, j, k, 0) = uH;
                    Up(i, j, k, 1) = vH;
                    Up(i, j, k, 2) = wH;
                }
    }

    // Uplus is reconstructed over the complete allocated range, including
    // ghosts, from U_H and PV after Ucell synchronization.  It exists
    // only on Fluid blocks, so a topology-wide component-copy halo would
    // incorrectly try to pack inactive non-fluid block storage.
}

void LunarSolver::UpdateFluidDerivedFields_()
{
    calc_PV();
    calc_Uplus();
}

void LunarSolver::UpdateMagneticDerivedFields_()
{
    // B_face is the CT state; J_cell is reconstructed from mimetic J_edge.
    calc_Bcell();
    Calc_J_Edge();
    calc_Jcell();
}

void LunarSolver::UpdateDerivedFields_()
{
    UpdateMagneticDerivedFields_();
    UpdateFluidDerivedFields_();
}
