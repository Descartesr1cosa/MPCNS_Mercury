
#include "MercurySolver.h"

void MercurySolver::Time_Advance()
{
    double Emag0, Emag1;
    if (control_.if_outres)
        Emag0 = ComputeMagEnergy_Cell_();

    // 1) dq/db set to ZERO
    ZeroRHS_();

    // 2) RHS for Flow
    AssembleRHS_Fluid_(); // = Scheme_U_ + AddSourceToRHS_Fluid()

    // 3) RHS for Magnetic fields：CT (Edge EMF -> Curl -> Face RHS)
    AssembleRHS_Induction_CT_();

    // 4) Euler 1st order Advance
    ApplyUpdate_Euler_(); // U += dt*RHS, B_face += dt*RHS_Bface

    // 5) 低密度/负压修复（尽量按 Fortran：邻域平均 + 重建 E）
    // RepairNonPhysical_();

    // 2) Hall 子步：只更新 Bface
    // {
    //     mercury_bound_.Sync("Ucell");
    //     mercury_bound_.Sync("Bface");

    //     calc_PV();
    //     calc_Uplus();
    //     calc_Bcell();
    //     // Calc_J_Edge();
    //     // calc_Jcell();
    // }

    // {
    //     mercury_bound_.Sync("Ucell");
    //     mercury_bound_.Sync("Bface");
    //     calc_PV();
    //     calc_Uplus();
    //     calc_Bcell();
    //     Calc_J_Edge();
    //     calc_Jcell();
    //     DebugFindExtremaInner({}, {"J_cell", "B_cell"});
    //     for (int i : {16})
    //         for (int k : {10})
    //             for (int j : {-1, 0, 1})
    //             {
    //                 DebugDumpPointFields(
    //                     1, 0, i, j, k,
    //                     {},
    //                     {"J_zeta", "J_xi"});

    //                 DebugDumpPointPartners(
    //                     1, 0, i, j, k,
    //                     {},
    //                     {"J_zeta", "J_xi"},
    //                     1,    // ngh
    //                     true, // include_topo: 只看精确耦合/接口对应点
    //                     true, // include_halo: 先不要混 halo
    //                     false // include_physical
    //                 );
    //             }
    // }
#if HALL_IMPLICIT == 1

    hall_implicit_.SolveOneStep(dt, control_.if_outres);
    calc_Bcell();

    if (control_.if_outres)
    {
        Emag1 = ComputeMagEnergy_Cell_();
        if (par_->GetInt("myid") == 0)
        {
            double dE = Emag1 - Emag0;
            double rel = dE / std::max(std::abs(Emag0), 1e-300);
            std::cout << "[HallOnlyEnergy] dt=" << dt
                      << " Emag0=" << Emag0
                      << " Emag1=" << Emag1
                      << " dE=" << dE
                      << " rel=" << rel
                      << std::endl
                      << std::endl
                      << std::endl;
        }
    }
#elif HALL_IMPLICIT == 0
    // Calc_J_Edge();
    // calc_Jcell();
    // const int nsub_max = par_->GetInt("hall_subcycle_max_steps");
    // double dt_hall_min = dt / (double)(nsub_max + 0.0);
    // for (int nsub = 0; nsub < nsub_max; nsub++)
    // {
    //     AddHallEdgeEMF_();

    //     AddHyperResistiveEdgeEMF_();

    //     AddSecondResistiveEdgeEMF_();

    //     // 2) E_edge 也要进边界/halo（否则 curl 更新 B_face 时边界附近会乱）
    //     mercury_bound_.Sync("Ehall"); // 你需要加一个 group：fields={E_xi,E_eta,E_zeta}

    //     // 3) curl(E_edge) -> RHS_Bface，然后 Bface += dt*RHS
    //     for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    //     {
    //         auto &Exi = fld_->field(fid_.fid_Ehall.xi, ib);
    //         auto &Eeta = fld_->field(fid_.fid_Ehall.eta, ib);
    //         auto &Eze = fld_->field(fid_.fid_Ehall.zeta, ib);

    //         auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
    //         auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
    //         auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);
    //         if (!Exi.is_allocated())
    //             continue;

    //         // multiper 的符号你要和你的取向一致：
    //         // 一般 CT 是 B^{n+1} = B^n - dt * curl(E)
    //         CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze, RHSBxi, RHSBeta, RHSBze, /*multiper=*/-1.0);
    //     }

    //     ApplyUpdate_Euler_BfaceOnly_(dt_hall_min);

    //     mercury_bound_.Sync("Bface");
    //     calc_Bcell();
    //     Calc_J_Edge();
    //     calc_Jcell();
    // }
#endif
}

void MercurySolver::ZeroRHS_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSN = fld_->field(fid_.fid_RHS_Na, ib);
        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_et = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_ze = fld_->field(fid_.fid_RHS_b.zeta, ib);

        Int3 lo; // = Jac.inner_lo();
        Int3 hi; // = Jac.inner_hi();

        if (RHSH.is_allocated())
        {
            lo = RHSH.inner_lo();
            hi = RHSH.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            RHSH(i, j, k, m) = 0.0;
        }

        if (RHSN.is_allocated())
        {
            lo = RHSN.inner_lo();
            hi = RHSN.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            RHSN(i, j, k, m) = 0.0;
        }

        if (RHSB_xi.is_allocated())
        {
            lo = RHSB_xi.inner_lo();
            hi = RHSB_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_xi(i, j, k, 0) = 0.0;
        }

        if (RHSB_et.is_allocated())
        {
            lo = RHSB_et.inner_lo();
            hi = RHSB_et.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_et(i, j, k, 0) = 0.0;
        }

        if (RHSB_ze.is_allocated())
        {
            lo = RHSB_ze.inner_lo();
            hi = RHSB_ze.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_ze(i, j, k, 0) = 0.0;
        }
    }
}

void MercurySolver::ApplyUpdate_Euler_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);

        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UN = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &Ub_xi = fld_->field(fid_.fid_B.xi, ib);
        FieldBlock &Ub_eta = fld_->field(fid_.fid_B.eta, ib);
        FieldBlock &Ub_zeta = fld_->field(fid_.fid_B.zeta, ib);

        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSN = fld_->field(fid_.fid_RHS_Na, ib);
        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_eta = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_zeta = fld_->field(fid_.fid_RHS_b.zeta, ib);

        Int3 lo; //= Jac.inner_lo();
        Int3 hi; //= Jac.inner_hi();

        if (RHSH.is_allocated())
        {
            lo = UH.inner_lo();
            hi = UH.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            UH(i, j, k, m) += dt * RHSH(i, j, k, m);
        }

        if (RHSN.is_allocated())
        {
            lo = UN.inner_lo();
            hi = UN.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            UN(i, j, k, m) += dt * RHSN(i, j, k, m);
        }

        if (RHSB_xi.is_allocated())
        {
            lo = Ub_xi.inner_lo();
            hi = Ub_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_xi(i, j, k, 0) += dt * RHSB_xi(i, j, k, 0);
        }

        if (RHSB_eta.is_allocated())
        {
            lo = Ub_eta.inner_lo();
            hi = Ub_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_eta(i, j, k, 0) += dt * RHSB_eta(i, j, k, 0);
        }

        if (RHSB_zeta.is_allocated())
        {
            lo = Ub_zeta.inner_lo();
            hi = Ub_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_zeta(i, j, k, 0) += dt * RHSB_zeta(i, j, k, 0);
        }
    }
}
