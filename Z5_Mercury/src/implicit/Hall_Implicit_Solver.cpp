#include "4_Hall_Implicit.h"

#if HALL_IMPLICIT == 1

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iomanip>

#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"
#include "1_Boundary.h"
#include "operators/Vector.h"

namespace
{
    void setup_like_face_field(Scalar &buf, FieldBlock &F)
    {
        if (!F.is_allocated())
            return;

        const Int3 lo = F.get_lo();
        const Int3 hi = F.get_hi();

        const int ghost = -lo.i;
        if (lo.i != -ghost || lo.j != -ghost || lo.k != -ghost)
        {
            throw std::runtime_error(
                "ImplicitHallSolver: face field is not in standard ghost form.");
        }

        buf.SetSize(hi.i - lo.i, hi.j - lo.j, hi.k - lo.k, ghost);
    }
} // namespace

ImplicitHallSolver::~ImplicitHallSolver()
{
    FinalizePetsc();
}

void ImplicitHallSolver::Setup(Grid *grd,
                               TOPO::Topology *topo,
                               Field *fld,
                               Halo *halo,
                               Param *par,
                               MercuryBoundary *bound,
                               const SolverFields &fid,
                               const TOPO::TopologyEquiv &equiv,
                               const HALO_OWNER::EdgeOwnerSyncPattern &owner_pat,
                               std::vector<HallFaceScratchBlock_> *hall_face_scratch)
{
    grd_ = grd;
    topo_ = topo;
    fld_ = fld;
    halo_ = halo;
    par_ = par;
    bound_ = bound;

    fid_ = fid;
    equiv_ = equiv;
    owner_pat_ = owner_pat;

    if (!fld_ || !bound_)
        throw std::runtime_error("ImplicitHallSolver::Setup: null core pointer.");

    x_local_.resize(static_cast<size_t>(equiv_.n_local_edge_owner), 0.0);
    eh_pred_local_.resize(static_cast<size_t>(equiv_.n_local_edge_owner), 0.0);

    HALO_OWNER::gather_local_owner_edges_sorted(equiv_, owner_edges_sorted_);
    if (static_cast<int>(owner_edges_sorted_.size()) != equiv_.n_local_edge_owner)
    {
        throw std::runtime_error(
            "ImplicitHallSolver::Setup: owner_edges_sorted_ size mismatch.");
    }

    hall_face_scratch_ = hall_face_scratch;
    SetupBfaceSnapshotStorage_();
}

void ImplicitHallSolver::CheckReady_() const
{
    if (!fld_ || !bound_)
        throw std::runtime_error("ImplicitHallSolver not setup.");
    if (!hall_face_scratch_)
        throw std::runtime_error("ImplicitHallSolver hall scratch is not setup.");
    if (!cb_.sync_Bface || !cb_.sync_Ehalledge || !cb_.calc_PV ||
        !cb_.calc_Uplus || !cb_.build_Ehall_from_current_B)
        throw std::runtime_error("ImplicitHallSolver callbacks are not fully bound.");
    if (use_shell_pc_)
    {
        if (!cb_.calc_Bcell_from_current_Bface ||
            !cb_.FillFrozenBflatFromCurrentBcell_ ||
            !cb_.FillFrozenAlphaFlatCell_ ||
            !cb_.sync_dEedge ||
            !cb_.sync_dBface)
        {
            throw std::runtime_error("Shell PC callbacks are not fully bound.");
        }
    }
}

void ImplicitHallSolver::SetupBfaceSnapshotStorage_()
{
    const int nb = fld_->num_blocks();
    Bstar_xi_.resize(nb);
    Bstar_eta_.resize(nb);
    Bstar_ze_.resize(nb);

    for (int ib = 0; ib < nb; ++ib)
    {
        setup_like_face_field(Bstar_xi_[ib], fld_->field(fid_.fid_B.xi, ib));
        setup_like_face_field(Bstar_eta_[ib], fld_->field(fid_.fid_B.eta, ib));
        setup_like_face_field(Bstar_ze_[ib], fld_->field(fid_.fid_B.zeta, ib));
    }
}

void ImplicitHallSolver::InitializePetsc()
{
    PetscBool is_init = PETSC_FALSE;
    PetscInitialized(&is_init);
    if (!is_init)
    {
        throw std::runtime_error("PETSc not initialized before ImplicitHallSolver::InitializePetsc()");
    }

    if (petsc_ready_)
        return;
    CheckReady_();
    CreatePetscObjects_();
    petsc_ready_ = true;
}

void ImplicitHallSolver::FinalizePetsc()
{
    DestroyPetscObjects_();
    petsc_ready_ = false;
}

void ImplicitHallSolver::CreatePetscObjects_()
{
    const PetscInt nloc = static_cast<PetscInt>(equiv_.n_local_edge_owner);
    const PetscInt nglb = static_cast<PetscInt>(equiv_.n_global_edge_owner);

    VecCreateMPI(PETSC_COMM_WORLD, nloc, nglb, &X_);
    VecDuplicate(X_, &F_);

    // shell PC 自己的工作向量：分布必须和 unknown 一致
    if (use_shell_pc_)
    {
        VecDuplicate(X_, &pc_q_);
        VecDuplicate(X_, &pc_res_);
    }

    SNESCreate(PETSC_COMM_WORLD, &snes_);
    SNESSetFunction(snes_, F_, &ImplicitHallSolver::FormFunction_, this);

    MatCreateShell(PETSC_COMM_WORLD, nloc, nloc, nglb, nglb, this, &Jshell_);
    MatShellSetOperation(Jshell_, MATOP_MULT,
                         (void (*)(void))&ImplicitHallSolver::MatMult_WhistlerShell_);
    SNESSetJacobian(snes_, Jshell_, Jshell_, &ImplicitHallSolver::FormJacobian_, this);

    SNESSetType(snes_, SNESNEWTONLS);

    SNESGetKSP(snes_, &ksp_);
    KSPSetType(ksp_, KSPGMRES);

    KSPGetPC(ksp_, &pc_);
    if (use_shell_pc_)
    {
        PCSetType(pc_, PCSHELL);
        PCShellSetContext(pc_, this);
        PCShellSetApply(pc_, &ImplicitHallSolver::PCApplyWhistlerP0_Shell_);
        PCShellSetName(pc_, "WhistlerP0");
    }
    else
    {
        PCSetType(pc_, PCNONE);
    }

    SNESSetFromOptions(snes_);
}

void ImplicitHallSolver::DestroyPetscObjects_()
{
    p0_frozen_ready_ = false;

    if (pc_res_)
    {
        VecDestroy(&pc_res_);
        pc_res_ = nullptr;
    }
    if (pc_q_)
    {
        VecDestroy(&pc_q_);
        pc_q_ = nullptr;
    }

    // 只是 alias，不能 destroy
    pc_ = nullptr;
    ksp_ = nullptr;

    if (Jmf_)
    {
        MatDestroy(&Jmf_);
        Jmf_ = nullptr;
    }
    if (F_)
    {
        VecDestroy(&F_);
        F_ = nullptr;
    }
    if (X_)
    {
        VecDestroy(&X_);
        X_ = nullptr;
    }
    if (snes_)
    {
        SNESDestroy(&snes_);
        snes_ = nullptr;
    }
    if (Jshell_)
    {
        MatDestroy(&Jshell_);
        Jshell_ = nullptr;
    }
}

void ImplicitHallSolver::SolveOneStep(double dt, bool if_outres)
{
    CheckReady_();
    if (!petsc_ready_)
        InitializePetsc();

    dt_ = dt;

    // 1) snapshot B*
    SnapshotCurrentBface_();
    ClearEdgeTriplet_(fid_.fid_Ehall);

    p0_frozen_ready_ = false;
    if (use_shell_pc_)
        PrepareWhistlerP0FrozenState_();

    // 2) Initial guess: current explicit/previous Ehall.
    LoadCurrentEhallIntoSolution_();

    // 3) solve
    SNESSolve(snes_, nullptr, X_);

    // 4) Write solved Ehall and apply the full-step CT update.
    ApplySolvedEhallToBface_();

    if (if_outres)
        PrintSolveDiagnostics_();

    cb_.sync_Bface();
}

void ImplicitHallSolver::LoadCurrentEhallIntoSolution_()
{
    HALO_OWNER::pack_owner_edge_1form_local(
        *fld_, fid_.fid_Ehall, equiv_, owner_edges_sorted_, x_local_);

    PetscScalar *xarr = nullptr;
    VecGetArray(X_, &xarr);
    for (PetscInt i = 0; i < static_cast<PetscInt>(x_local_.size()); ++i)
        xarr[i] = x_local_[static_cast<size_t>(i)];
    VecRestoreArray(X_, &xarr);
}

void ImplicitHallSolver::CurlEhallToRhsB_()
{
    cb_.sync_Ehalledge();
    ClearFaceTriplet_(fid_.fid_RHS_b);

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_Ehall.xi, ib);
        auto &Eet = fld_->field(fid_.fid_Ehall.eta, ib);
        auto &Eze = fld_->field(fid_.fid_Ehall.zeta, ib);
        auto &Rxi = fld_->field(fid_.fid_RHS_b.xi, ib);
        auto &Ret = fld_->field(fid_.fid_RHS_b.eta, ib);
        auto &Rze = fld_->field(fid_.fid_RHS_b.zeta, ib);

        if (!Exi.is_allocated())
            continue;

        CTOperators::CurlEdgeToFace(ib, Exi, Eet, Eze, Rxi, Ret, Rze, -1.0);
    }
}

void ImplicitHallSolver::ApplySolvedEhallToBface_()
{
    UnpackVecToEhallField_(X_);
    RestoreCurrentBfaceFromSnapshot_();
    CurlEhallToRhsB_();
    AddFaceInnerFromRHS_(fid_.fid_B.xi, fid_.fid_RHS_b.xi, dt_);
    AddFaceInnerFromRHS_(fid_.fid_B.eta, fid_.fid_RHS_b.eta, dt_);
    AddFaceInnerFromRHS_(fid_.fid_B.zeta, fid_.fid_RHS_b.zeta, dt_);
}

void ImplicitHallSolver::PrintSolveDiagnostics_()
{
    SNESConvergedReason reason;
    PetscInt its;
    PetscReal fnorm = 0.0;

    SNESGetConvergedReason(snes_, &reason);
    SNESGetIterationNumber(snes_, &its);
    SNESGetFunctionNorm(snes_, &fnorm);

    const double maxEhall = MaxAbsTriplet_(fid_.fid_Ehall);
    const double maxRHSB = MaxAbsTriplet_(fid_.fid_RHS_b);
    const double maxdBhall_est = dt_ * maxRHSB;

    if (par_->GetInt("myid") == 0)
    {
        std::cout << std::scientific << std::setprecision(6)
                  << "[HallImplicit] reason=" << reason
                  << "  its=" << its
                  << "  |F|_2=" << fnorm
                  << "  max|Ehall|=" << maxEhall
                  << "  dt*max|RHS_B|=" << maxdBhall_est
                  << std::endl;
    }
}

double ImplicitHallSolver::MaxAbsTriplet_(const IdTriplet &fid_triplet)
{
    double local_max = 0.0;

    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        for (int fid : {fid_triplet.xi, fid_triplet.eta, fid_triplet.zeta})
        {
            auto &F = fld_->field(fid, ib);
            if (!F.is_allocated())
                continue;

            Int3 lo = F.inner_lo(), hi = F.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        local_max = std::max(local_max, std::abs(F(i, j, k, 0)));
        }
    }

    double global_max = 0.0;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_max;
}

void ImplicitHallSolver::UnpackVecToTempDEdgeField_(Vec X)
{
    ClearEdgeTriplet_(fid_.fid_dE);

    const PetscScalar *xarr = nullptr;
    VecGetArrayRead(X, &xarr);

    for (PetscInt lid = 0; lid < static_cast<PetscInt>(owner_edges_sorted_.size()); ++lid)
    {
        const auto &e = owner_edges_sorted_[static_cast<size_t>(lid)];
        const double val = static_cast<double>(xarr[lid]);

        int fid_edge = -1;
        switch (e.dir)
        {
        case 1:
            fid_edge = fid_.fid_dE.xi;
            break;
        case 2:
            fid_edge = fid_.fid_dE.eta;
            break;
        case 3:
            fid_edge = fid_.fid_dE.zeta;
            break;
        default:
            VecRestoreArrayRead(X, &xarr);
            throw std::runtime_error("UnpackVecToTempDEdgeField_: invalid edge axis.");
        }

        auto &F = fld_->field(fid_edge, e.gblock);
        F(e.i, e.j, e.k, 0) = val;
    }

    VecRestoreArrayRead(X, &xarr);

    cb_.sync_dEedge();
}

PetscErrorCode ImplicitHallSolver::PCApplyWhistlerP0_Shell_(PC pc, Vec in, Vec out)
{
    void *ctx = nullptr;
    PetscCall(PCShellGetContext(pc, &ctx));

    auto *S = static_cast<ImplicitHallSolver *>(ctx);
    PetscCheck(S, PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL,
               "Null ImplicitHallSolver context in PC shell");

    PetscCall(S->ApplyWhistlerP0ApproxInverse_(in, out));
    return 0;
}

void ImplicitHallSolver::PrepareWhistlerP0FrozenState_()
{
    RestoreCurrentBfaceFromSnapshot_();
    cb_.sync_Bface();

    cb_.calc_Bcell_from_current_Bface();
    cb_.FillFrozenBflatFromCurrentBcell_();
    cb_.FillFrozenAlphaFlatCell_();

    p0_frozen_ready_ = true;
}

PetscErrorCode ImplicitHallSolver::ApplyWhistlerP0ApproxInverse_(Vec in, Vec out)
{
    // PetscCall(VecCopy(in, out));
    // return 0;

    // ------------------------------------------------------------
    // 目标：给定 residual r = in，近似求解
    //
    //      P0 z = r
    //
    // 其中
    //      P0 = I - (whistler-only dominant term)
    //
    // 先用最简单的一步/两步 Richardson:
    //
    //      z^{m+1} = z^m + omega * (r - P0 z^m)
    //
    // 初始猜测先取 z^0 = r
    // ------------------------------------------------------------

    // 0) frozen 状态只准备一次（后面可升级成每个 Newton 步更新）
    if (!p0_frozen_ready_)
    {
        PrepareWhistlerP0FrozenState_();
    }

    // 1) 初值：z = r
    PetscCall(VecCopy(in, out));

    // 2) Richardson sweeps
    for (int sweep = 0; sweep < p0_nsweeps_; ++sweep)
    {
        // q = P0 z
        PetscCall(ApplyWhistlerP0Operator_(out, pc_q_));

        // pc_res_ = r - q
        PetscCall(VecWAXPY(pc_res_, -1.0, pc_q_, in)); // pc_res_ = in - pc_q_

        // z <- z + omega * (r - P0 z)
        PetscCall(VecAXPY(out, p0_omega_, pc_res_));
    }

    return 0;
}

PetscErrorCode ImplicitHallSolver::ApplyWhistlerP0Operator_(Vec in, Vec out)
{
    // ------------------------------------------------------------
    // 输入:  in  = 某个 edge-space 向量 z
    // 输出:  out = P0 z
    //
    // P0 z = z - dEpre
    //
    // 其中 dEpre 由下面的线性 whistler 链产生：
    //
    //   z(edge)
    //     -> dB(face)
    //     -> dJedge(edge)
    //     -> dJcell(cell)
    //     -> dEhc(cell)
    //     -> dEface(face)
    //     -> dEpre(edge)
    // ------------------------------------------------------------

    // A) z -> dB
    BuildDeltaBFromDeltaE_(in);

    // B) dB -> dJedge
    Calc_DeltaJ_Edge_FromDeltaB_();

    // C) dJedge -> dJcell
    Calc_DeltaJcell_FromDeltaJedge_Frozen_();

    // D) dJcell + frozen(Bflat, alpha_flat) -> dEhc(cell)
    BuildLinearHallCellEMF_();

    // E) dEhc(cell) -> dEface(face)
    BuildLinearHallFaceEMF_();

    // F) dEface(face) -> dEpre(edge)
    AssembleLinearHallEdgeEMF_();

    // G) out = in - dEpre
    WriteP0Action_(in, out);

    return 0;
}

void ImplicitHallSolver::WriteP0Action_(Vec in, Vec out)
{
    // out = in
    PetscCallAbort(PETSC_COMM_WORLD, VecCopy(in, out));

    // out -= Pack(dEpre)
    SubtractPackedTempDEpreFromVec_(out);
}

void ImplicitHallSolver::BuildDeltaBFromDeltaE_(Vec in)
{
    // 1) 把 Krylov 向量 unpack 到临时 edge 场，并做 owner/non-owner + halo 同步
    UnpackVecToTempDEdgeField_(in);

    // 2) 清零临时 dB face
    ClearFaceTriplet_(fid_.fid_dB);

    // 3) dB = -theta * dt * curl(dE)
    const double factor = -theta_ * dt_;
    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_dE.xi, ib);
        auto &Eet = fld_->field(fid_.fid_dE.eta, ib);
        auto &Eze = fld_->field(fid_.fid_dE.zeta, ib);

        auto &Bxi = fld_->field(fid_.fid_dB.xi, ib);
        auto &Bet = fld_->field(fid_.fid_dB.eta, ib);
        auto &Bze = fld_->field(fid_.fid_dB.zeta, ib);

        if (!Exi.is_allocated())
            continue;

        CTOperators::CurlEdgeToFace(ib, Exi, Eet, Eze, Bxi, Bet, Bze, factor);
    }

    // 4) 补齐 dB halo/interface，供下一步 dB -> dJ 使用
    cb_.sync_dBface();
}

void ImplicitHallSolver::Calc_DeltaJ_Edge_FromDeltaB_()
{
    ClearEdgeTriplet_(fid_.fid_dJ);

    //  ComputeJ_AtEdges_Inner_();
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Bxi = fld_->field(fid_.fid_dB.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_dB.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_dB.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_dJ.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_dJ.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_dJ.zeta, iblk);

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

    cb_.sync_dJedge();
}

void ImplicitHallSolver::Calc_DeltaJcell_FromDeltaJedge_Frozen_()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &Jcell = fld_->field(fid_.fid_dJcell, ib);

        auto &Jxi = fld_->field(fid_.fid_dJ.xi, ib);
        auto &Jeta = fld_->field(fid_.fid_dJ.eta, ib);
        auto &Jzeta = fld_->field(fid_.fid_dJ.zeta, ib);

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

    cb_.sync_dJcell();
}

void ImplicitHallSolver::BuildLinearHallCellEMF_()
{
    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &dJc = fld_->field(fid_.fid_dJcell, ib);

        auto &Bflat = (*hall_face_scratch_)[ib].Bflat;      // 3-comp cell
        auto &alpha = (*hall_face_scratch_)[ib].alpha_flat; // scalar cell
        auto &dEhc = (*hall_face_scratch_)[ib].dEhc;        // 3-comp cell

        if (!dJc.is_allocated())
            continue;

        Int3 lo = dJc.inner_lo();
        Int3 hi = dJc.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double a = alpha(i, j, k);

                    const double dJx = dJc(i, j, k, 0);
                    const double dJy = dJc(i, j, k, 1);
                    const double dJz = dJc(i, j, k, 2);

                    const double Bx = Bflat(i, j, k, 0);
                    const double By = Bflat(i, j, k, 1);
                    const double Bz = Bflat(i, j, k, 2);

                    dEhc(i, j, k, 0) = a * (dJy * Bz - dJz * By);
                    dEhc(i, j, k, 1) = a * (dJz * Bx - dJx * Bz);
                    dEhc(i, j, k, 2) = a * (dJx * By - dJy * Bx);
                }
    }
}

void ImplicitHallSolver::BuildLinearHallFaceEMF_()
{
    ClearFaceTriplet_(fid_.fid_Eface);

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &scratch = (*hall_face_scratch_)[ib];

        auto &dEhc = scratch.dEhc; // 3-comp cell
        auto &Pxi = fld_->field(fid_.Face_projector.xi, ib);
        auto &Pet = fld_->field(fid_.Face_projector.eta, ib);
        auto &Pze = fld_->field(fid_.Face_projector.zeta, ib);

        auto &Efxi = fld_->field(fid_.fid_Eface.xi, ib); // 3-comp face
        auto &Efet = fld_->field(fid_.fid_Eface.eta, ib);
        auto &Efze = fld_->field(fid_.fid_Eface.zeta, ib);

        if (!Efxi.is_allocated() || !Efet.is_allocated() || !Efze.is_allocated())
            continue;

        // ============================================================
        // xi-face : average adjacent cells, then apply cached projector
        // ============================================================
        {
            Int3 lo = Efxi.inner_lo();
            Int3 hi = Efxi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int iL = i - 1;
                        const int iR = i;

                        const double Ex0 = 0.5 * (dEhc(iL, j, k, 0) + dEhc(iR, j, k, 0));
                        const double Ey0 = 0.5 * (dEhc(iL, j, k, 1) + dEhc(iR, j, k, 1));
                        const double Ez0 = 0.5 * (dEhc(iL, j, k, 2) + dEhc(iR, j, k, 2));

                        const double Pxx = Pxi(i, j, k, 0);
                        const double Pxy = Pxi(i, j, k, 1);
                        const double Pxz = Pxi(i, j, k, 2);
                        const double Pyy = Pxi(i, j, k, 3);
                        const double Pyz = Pxi(i, j, k, 4);
                        const double Pzz = Pxi(i, j, k, 5);

                        Efxi(i, j, k, 0) = Pxx * Ex0 + Pxy * Ey0 + Pxz * Ez0;
                        Efxi(i, j, k, 1) = Pxy * Ex0 + Pyy * Ey0 + Pyz * Ez0;
                        Efxi(i, j, k, 2) = Pxz * Ex0 + Pyz * Ey0 + Pzz * Ez0;
                    }
        }

        // ============================================================
        // eta-face
        // ============================================================
        {
            Int3 lo = Efet.inner_lo();
            Int3 hi = Efet.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int jL = j - 1;
                        const int jR = j;

                        const double Ex0 = 0.5 * (dEhc(i, jL, k, 0) + dEhc(i, jR, k, 0));
                        const double Ey0 = 0.5 * (dEhc(i, jL, k, 1) + dEhc(i, jR, k, 1));
                        const double Ez0 = 0.5 * (dEhc(i, jL, k, 2) + dEhc(i, jR, k, 2));

                        const double Pxx = Pet(i, j, k, 0);
                        const double Pxy = Pet(i, j, k, 1);
                        const double Pxz = Pet(i, j, k, 2);
                        const double Pyy = Pet(i, j, k, 3);
                        const double Pyz = Pet(i, j, k, 4);
                        const double Pzz = Pet(i, j, k, 5);

                        Efet(i, j, k, 0) = Pxx * Ex0 + Pxy * Ey0 + Pxz * Ez0;
                        Efet(i, j, k, 1) = Pxy * Ex0 + Pyy * Ey0 + Pyz * Ez0;
                        Efet(i, j, k, 2) = Pxz * Ex0 + Pyz * Ey0 + Pzz * Ez0;
                    }
        }

        // ============================================================
        // zeta-face
        // ============================================================
        {
            Int3 lo = Efze.inner_lo();
            Int3 hi = Efze.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int kL = k - 1;
                        const int kR = k;

                        const double Ex0 = 0.5 * (dEhc(i, j, kL, 0) + dEhc(i, j, kR, 0));
                        const double Ey0 = 0.5 * (dEhc(i, j, kL, 1) + dEhc(i, j, kR, 1));
                        const double Ez0 = 0.5 * (dEhc(i, j, kL, 2) + dEhc(i, j, kR, 2));

                        const double Pxx = Pze(i, j, k, 0);
                        const double Pxy = Pze(i, j, k, 1);
                        const double Pxz = Pze(i, j, k, 2);
                        const double Pyy = Pze(i, j, k, 3);
                        const double Pyz = Pze(i, j, k, 4);
                        const double Pzz = Pze(i, j, k, 5);

                        Efze(i, j, k, 0) = Pxx * Ex0 + Pxy * Ey0 + Pxz * Ez0;
                        Efze(i, j, k, 1) = Pxy * Ex0 + Pyy * Ey0 + Pyz * Ez0;
                        Efze(i, j, k, 2) = Pxz * Ex0 + Pyz * Ey0 + Pzz * Ez0;
                    }
        }
    }

    cb_.sync_dEface();
}

void ImplicitHallSolver::AssembleLinearHallEdgeEMF_()
{
    Int3 sub, sup;
    ClearEdgeTriplet_(fid_.fid_dEpre);

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        double3D &x = grd_->grids(iblk).x;
        double3D &y = grd_->grids(iblk).y;
        double3D &z = grd_->grids(iblk).z;

        Vec3 E, dr;

        // -----------------------------------------
        // dEpre_xi(edge) from face-eta and face-zeta
        // -----------------------------------------
        {
            auto &Exi = fld_->field(fid_.fid_dEpre.xi, iblk);
            auto &Ef_et = fld_->field(fid_.fid_Eface.eta, iblk);
            auto &Ef_ze = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Exi.inner_lo();
            sup = Exi.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_et(i, j, k, 0) + Ef_et(i, j, k - 1, 0) + Ef_ze(i, j, k, 0) + Ef_ze(i, j - 1, k, 0));
                        E.vec[1] = 0.25 * (Ef_et(i, j, k, 1) + Ef_et(i, j, k - 1, 1) + Ef_ze(i, j, k, 1) + Ef_ze(i, j - 1, k, 1));
                        E.vec[2] = 0.25 * (Ef_et(i, j, k, 2) + Ef_et(i, j, k - 1, 2) + Ef_ze(i, j, k, 2) + Ef_ze(i, j - 1, k, 2));

                        dr = {
                            x(i + 1, j, k) - x(i, j, k),
                            y(i + 1, j, k) - y(i, j, k),
                            z(i + 1, j, k) - z(i, j, k)};

                        Exi(i, j, k, 0) = E * dr;
                    }
        }

        // -----------------------------------------
        // dEpre_eta(edge) from face-xi and face-zeta
        // -----------------------------------------
        {
            auto &Eet = fld_->field(fid_.fid_dEpre.eta, iblk);
            auto &Ef_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &Ef_ze = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Eet.inner_lo();
            sup = Eet.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_xi(i, j, k, 0) + Ef_xi(i, j, k - 1, 0) + Ef_ze(i, j, k, 0) + Ef_ze(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (Ef_xi(i, j, k, 1) + Ef_xi(i, j, k - 1, 1) + Ef_ze(i, j, k, 1) + Ef_ze(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (Ef_xi(i, j, k, 2) + Ef_xi(i, j, k - 1, 2) + Ef_ze(i, j, k, 2) + Ef_ze(i - 1, j, k, 2));

                        dr = {
                            x(i, j + 1, k) - x(i, j, k),
                            y(i, j + 1, k) - y(i, j, k),
                            z(i, j + 1, k) - z(i, j, k)};

                        Eet(i, j, k, 0) = E * dr;
                    }
        }

        // -----------------------------------------
        // dEpre_zeta(edge) from face-xi and face-eta
        // -----------------------------------------
        {
            auto &Eze = fld_->field(fid_.fid_dEpre.zeta, iblk);
            auto &Ef_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &Ef_et = fld_->field(fid_.fid_Eface.eta, iblk);

            sub = Eze.inner_lo();
            sup = Eze.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_xi(i, j, k, 0) + Ef_xi(i, j - 1, k, 0) + Ef_et(i, j, k, 0) + Ef_et(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (Ef_xi(i, j, k, 1) + Ef_xi(i, j - 1, k, 1) + Ef_et(i, j, k, 1) + Ef_et(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (Ef_xi(i, j, k, 2) + Ef_xi(i, j - 1, k, 2) + Ef_et(i, j, k, 2) + Ef_et(i - 1, j, k, 2));

                        dr = {
                            x(i, j, k + 1) - x(i, j, k),
                            y(i, j, k + 1) - y(i, j, k),
                            z(i, j, k + 1) - z(i, j, k)};

                        Eze(i, j, k, 0) = E * dr;
                    }
        }
    }
}

void ImplicitHallSolver::SubtractPackedTempDEpreFromVec_(Vec out)
{
    HALO_OWNER::pack_owner_edge_1form_local(
        *fld_, fid_.fid_dEpre, equiv_, owner_edges_sorted_, x_local_);

    PetscScalar *outarr = nullptr;
    VecGetArray(out, &outarr);

    for (PetscInt i = 0; i < static_cast<PetscInt>(x_local_.size()); ++i)
        outarr[i] -= static_cast<PetscScalar>(x_local_[static_cast<size_t>(i)]);

    VecRestoreArray(out, &outarr);
}

PetscErrorCode ImplicitHallSolver::FormJacobian_(SNES, Vec X, Mat, Mat, void *ctx)
{
    auto *S = static_cast<ImplicitHallSolver *>(ctx);

    // Current first pass: keep one frozen whistler state per time step.
    if (!S->p0_frozen_ready_)
        S->PrepareWhistlerP0FrozenState_();
    return 0;
}

PetscErrorCode ImplicitHallSolver::MatMult_WhistlerShell_(Mat A, Vec in, Vec out)
{
    void *ctx = nullptr;
    PetscCall(MatShellGetContext(A, &ctx));

    auto *S = static_cast<ImplicitHallSolver *>(ctx);
    PetscCheck(S, PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL,
               "Null ImplicitHallSolver context in shell MatMult");

    PetscCall(S->ApplyWhistlerJv_(in, out));
    return 0;
}

PetscErrorCode ImplicitHallSolver::ApplyWhistlerJv_(Vec in, Vec out)
{
    // First-pass Jacobian action: reuse the current P0 operator.
    return ApplyWhistlerP0Operator_(in, out);
}
#endif
