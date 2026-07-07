#include "4_Hall_Implicit.h"

#if HALL_IMPLICIT == 1

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iomanip>

#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "1_Boundary.h"
#include "operators/Vector.h"

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
                               const TOPO::Topology &equiv,
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

    x_local_.resize(static_cast<size_t>(equiv_.edges.n_local_owner), 0.0);
    eh_pred_local_.resize(static_cast<size_t>(equiv_.edges.n_local_owner), 0.0);

    HALO_OWNER::gather_local_owner_edges_sorted(equiv_, owner_edges_sorted_);
    if (static_cast<int>(owner_edges_sorted_.size()) != equiv_.edges.n_local_owner)
    {
        throw std::runtime_error(
            "ImplicitHallSolver::Setup: owner_edges_sorted_ size mismatch.");
    }

    // Setup B* snapshot
    {
        const int nb = fld_->num_blocks();
        Bstar_xi_.resize(nb);
        Bstar_eta_.resize(nb);
        Bstar_ze_.resize(nb);

        auto setup_like_face = [](Scalar &buf, auto &F)
        {
            if (!F.is_allocated())
                return false;

            const Int3 lo = F.get_lo();
            const Int3 hi = F.get_hi();

            const int ghost = -lo.i;
            if (lo.i != -ghost || lo.j != -ghost || lo.k != -ghost)
            {
                throw std::runtime_error(
                    "SetupBfaceSnapshotScratch_: face field not in standard ghost form.");
            }

            const int dim1 = hi.i - lo.i;
            const int dim2 = hi.j - lo.j;
            const int dim3 = hi.k - lo.k;

            buf.SetSize(dim1, dim2, dim3, ghost);
            return true;
        };

        for (int ib = 0; ib < nb; ++ib)
        {
            auto &bxi = fld_->field(fid_.fid_B.xi, ib);
            auto &bet = fld_->field(fid_.fid_B.eta, ib);
            auto &bze = fld_->field(fid_.fid_B.zeta, ib);

            auto &xi = Bstar_xi_[ib];
            auto &eta = Bstar_eta_[ib];
            auto &zeta = Bstar_ze_[ib];
            setup_like_face(xi, bxi);
            setup_like_face(eta, bet);
            setup_like_face(zeta, bze);
        }
    }

    hall_face_scratch_ = hall_face_scratch;
}

void ImplicitHallSolver::CheckReady_() const
{
    if (!fld_ || !bound_)
        throw std::runtime_error("ImplicitHallSolver not setup.");
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
    const PetscInt nloc = static_cast<PetscInt>(equiv_.edges.n_local_owner);
    const PetscInt nglb = static_cast<PetscInt>(equiv_.edges.n_global_owner);

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

    // MatCreateSNESMF(snes_, &Jmf_);
    // SNESSetJacobian(snes_, Jmf_, Jmf_, MatMFFDComputeJacobian, nullptr);

    // const PetscInt nloc = static_cast<PetscInt>(equiv_.edges.n_local_owner);
    // const PetscInt nglb = static_cast<PetscInt>(equiv_.edges.n_global_owner);
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

    // const PetscInt nloc = static_cast<PetscInt>(equiv_.edges.n_local_owner);
    // const PetscInt nglb = static_cast<PetscInt>(equiv_.edges.n_global_owner);

    // VecCreateMPI(PETSC_COMM_WORLD, nloc, nglb, &X_);
    // VecDuplicate(X_, &F_);

    // SNESCreate(PETSC_COMM_WORLD, &snes_);
    // SNESSetFunction(snes_, F_, &ImplicitHallSolver::FormFunction_, this);

    // // matrix-free Jacobian
    // MatCreateSNESMF(snes_, &Jmf_);
    // SNESSetJacobian(snes_, Jmf_, Jmf_, MatMFFDComputeJacobian, nullptr);

    // SNESSetType(snes_, SNESNEWTONLS);

    // KSP ksp = nullptr;
    // SNESGetKSP(snes_, &ksp);
    // KSPSetType(ksp, KSPGMRES);

    // PC pc = nullptr;
    // KSPGetPC(ksp, &pc);
    // // PCSetType(pc, PCJACOBI);
    // PCSetType(pc, PCNONE);

    // SNESSetFromOptions(snes_);
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

    // 这一句非常关键：每个时间步都要重新冻结
    p0_frozen_ready_ = false;

    // 可选：提前准备，而不是等第一次 PCApply 时才做
    if (use_shell_pc_)
        PrepareWhistlerP0FrozenState_();

    // 2) 初值：拿当前显式/上一步的 Ehall 当 guess
    HALO_OWNER::pack_owner_edge_1form_local(
        *fld_, fid_.fid_Ehall, equiv_, owner_edges_sorted_, x_local_);

    PetscScalar *xarr = nullptr;
    VecGetArray(X_, &xarr);
    for (PetscInt i = 0; i < static_cast<PetscInt>(x_local_.size()); ++i)
        xarr[i] = x_local_[static_cast<size_t>(i)];
    VecRestoreArray(X_, &xarr);

    // 3) solve
    SNESSolve(snes_, nullptr, X_);

    // 4) 把最终解写回 Ehall
    UnpackVecToEhallField_(X_);

    // 5) 用最终 Ehall 做整步 CT 更新：
    //    B^{n+1} = B* + dt * RHS(Ehall), RHS=-curl(Ehall)
    RestoreCurrentBfaceFromSnapshot_();

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

    if (if_outres)
    {
        SNESConvergedReason reason;
        PetscInt its;
        PetscReal fnorm = 0.0;

        SNESGetConvergedReason(snes_, &reason);
        SNESGetIterationNumber(snes_, &its);
        SNESGetFunctionNorm(snes_, &fnorm);

        auto max_abs_face_delta = [&](int fid,
                                      std::vector<Scalar> &snap) -> double
        {
            double local_max = 0.0;

            for (int ib = 0; ib < fld_->num_blocks(); ++ib)
            {
                auto &F = fld_->field(fid, ib);
                if (!F.is_allocated())
                    continue;

                Int3 lo = F.inner_lo(), hi = F.inner_hi();
                auto &buf = snap[static_cast<size_t>(ib)];

                size_t t = 0;
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k, ++t)
                            local_max = std::max(local_max, std::abs(F(i, j, k, 0) - buf(i, j, k)));
            }

            double global_max = 0.0;
            MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            return global_max;
        };

        auto max_dB_hall = [&]() -> double
        {
            const double dxi = max_abs_face_delta(fid_.fid_B.xi, Bstar_xi_);
            const double det = max_abs_face_delta(fid_.fid_B.eta, Bstar_eta_);
            const double dze = max_abs_face_delta(fid_.fid_B.zeta, Bstar_ze_);
            return std::max(dxi, std::max(det, dze));
        };

        const double maxEhall = MaxAbsTriplet_(fid_.fid_Ehall);
        const double maxB = MaxAbsTriplet_(fid_.fid_B);
        const double maxRHSB = MaxAbsTriplet_(fid_.fid_RHS_b);
        const double maxdBhall_est = dt_ * maxRHSB;
        int rank = par_->GetInt("myid");
        if (rank == 0)
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

    AddFaceInnerFromRHS_(fid_.fid_B.xi, fid_.fid_RHS_b.xi, dt_);
    AddFaceInnerFromRHS_(fid_.fid_B.eta, fid_.fid_RHS_b.eta, dt_);
    AddFaceInnerFromRHS_(fid_.fid_B.zeta, fid_.fid_RHS_b.zeta, dt_);

    cb_.sync_Bface();
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

        auto &Hodge_star_2form_to_1form_face_xi = fld_->field(fid_.Hodge_star_2form_to_1form_face.xi, iblk);
        auto &Hodge_star_2form_to_1form_face_eta = fld_->field(fid_.Hodge_star_2form_to_1form_face.eta, iblk);
        auto &Hodge_star_2form_to_1form_face_zeta = fld_->field(fid_.Hodge_star_2form_to_1form_face.zeta, iblk);

        auto &Hodge_star_inverse_2form_to_1form_edge_xi = fld_->field(fid_.Hodge_star_inverse_2form_to_1form_edge.xi, iblk);
        auto &Hodge_star_inverse_2form_to_1form_edge_eta = fld_->field(fid_.Hodge_star_inverse_2form_to_1form_edge.eta, iblk);
        auto &Hodge_star_inverse_2form_to_1form_edge_zeta = fld_->field(fid_.Hodge_star_inverse_2form_to_1form_edge.zeta, iblk);

        // compute J (edge 1-form) from face B (2-form)
        // multiper 用 +1.0 J =curl B。
        CTOperators::CurlAdjFaceToEdge(iblk,
                                       Bxi, Beta, Bzeta,
                                       Hodge_star_2form_to_1form_face_xi, Hodge_star_2form_to_1form_face_eta, Hodge_star_2form_to_1form_face_zeta,
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
                        Jxi(i, j, k, 0) *= (Hodge_star_inverse_2form_to_1form_edge_xi(i, j, k, 0));
        }
        {
            // Edge_eta
            Int3 lo = Jeta.inner_lo(), hi = Jeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jeta(i, j, k, 0) *= (Hodge_star_inverse_2form_to_1form_edge_eta(i, j, k, 0));
        }
        {
            // Edge_zeta
            Int3 lo = Jzeta.inner_lo(), hi = Jzeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jzeta(i, j, k, 0) *= (Hodge_star_inverse_2form_to_1form_edge_zeta(i, j, k, 0));
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

        auto &W = (*hall_face_scratch_)[ib].dJcell_w;

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

// void ImplicitHallSolver::Calc_DeltaJcell_FromDeltaJedge_Frozen_()
// {
//     const int nblock = fld_->num_blocks();

//     constexpr double eps = 1e-25;

//     for (int ib = 0; ib < nblock; ++ib)
//     {
//         auto &Jcell = fld_->field(fid_.fid_dJcell, ib);

//         auto &Jxi = fld_->field(fid_.fid_dJ.xi, ib);
//         auto &Jeta = fld_->field(fid_.fid_dJ.eta, ib);
//         auto &Jzeta = fld_->field(fid_.fid_dJ.zeta, ib);

//         auto &dl_xi = fld_->field("dl_xi", ib);
//         auto &dl_eta = fld_->field("dl_eta", ib);
//         auto &dl_zeta = fld_->field("dl_zeta", ib);

//         auto &x = grd_->grids(ib).x;
//         auto &y = grd_->grids(ib).y;
//         auto &z = grd_->grids(ib).z;

//         if (!Jcell.is_allocated() || !Jxi.is_allocated() || !Jeta.is_allocated() || !Jzeta.is_allocated())
//             continue;

//         auto dot3 = [&](double ax, double ay, double az,
//                         double bx, double by, double bz) -> double
//         {
//             return ax * bx + ay * by + az * bz;
//         };

//         auto unit_t_xi = [&](int i, int j, int k,
//                              double &tx, double &ty, double &tz)
//         {
//             const double L = std::max(dl_xi(i, j, k, 0), eps);
//             tx = (x(i + 1, j, k) - x(i, j, k)) / L;
//             ty = (y(i + 1, j, k) - y(i, j, k)) / L;
//             tz = (z(i + 1, j, k) - z(i, j, k)) / L;
//         };

//         auto unit_t_eta = [&](int i, int j, int k,
//                               double &tx, double &ty, double &tz)
//         {
//             const double L = std::max(dl_eta(i, j, k, 0), eps);
//             tx = (x(i, j + 1, k) - x(i, j, k)) / L;
//             ty = (y(i, j + 1, k) - y(i, j, k)) / L;
//             tz = (z(i, j + 1, k) - z(i, j, k)) / L;
//         };

//         auto unit_t_zeta = [&](int i, int j, int k,
//                                double &tx, double &ty, double &tz)
//         {
//             const double L = std::max(dl_zeta(i, j, k, 0), eps);
//             tx = (x(i, j, k + 1) - x(i, j, k)) / L;
//             ty = (y(i, j, k + 1) - y(i, j, k)) / L;
//             tz = (z(i, j, k + 1) - z(i, j, k)) / L;
//         };

//         struct Eq
//         {
//             double tx, ty, tz; // unit tangent
//             double rhs;        // J_edge / |dl|
//             double w;          // weight
//         };

//         Int3 lo = Jcell.inner_lo();
//         Int3 hi = Jcell.inner_hi();

//         for (int i = lo.i; i < hi.i; ++i)
//             for (int j = lo.j; j < hi.j; ++j)
//                 for (int k = lo.k; k < hi.k; ++k)
//                 {
//                     Eq eqs[12];
//                     int K = 0;

//                     auto push = [&](double tx, double ty, double tz,
//                                     double Jint, double L, double w = 1.0)
//                     {
//                         L = std::max(L, eps);
//                         eqs[K++] = {tx, ty, tz, Jint / L, w};
//                     };

//                     double tx, ty, tz;

//                     // =====================================================
//                     // 4 xi-edges around cell(i,j,k)
//                     // =====================================================
//                     unit_t_xi(i, j, k, tx, ty, tz);
//                     push(tx, ty, tz, Jxi(i, j, k, 0), dl_xi(i, j, k, 0));

//                     unit_t_xi(i, j + 1, k, tx, ty, tz);
//                     push(tx, ty, tz, Jxi(i, j + 1, k, 0), dl_xi(i, j + 1, k, 0));

//                     unit_t_xi(i, j, k + 1, tx, ty, tz);
//                     push(tx, ty, tz, Jxi(i, j, k + 1, 0), dl_xi(i, j, k + 1, 0));

//                     unit_t_xi(i, j + 1, k + 1, tx, ty, tz);
//                     push(tx, ty, tz, Jxi(i, j + 1, k + 1, 0), dl_xi(i, j + 1, k + 1, 0));

//                     // =====================================================
//                     // 4 eta-edges
//                     // =====================================================
//                     unit_t_eta(i, j, k, tx, ty, tz);
//                     push(tx, ty, tz, Jeta(i, j, k, 0), dl_eta(i, j, k, 0));

//                     unit_t_eta(i + 1, j, k, tx, ty, tz);
//                     push(tx, ty, tz, Jeta(i + 1, j, k, 0), dl_eta(i + 1, j, k, 0));

//                     unit_t_eta(i, j, k + 1, tx, ty, tz);
//                     push(tx, ty, tz, Jeta(i, j, k + 1, 0), dl_eta(i, j, k + 1, 0));

//                     unit_t_eta(i + 1, j, k + 1, tx, ty, tz);
//                     push(tx, ty, tz, Jeta(i + 1, j, k + 1, 0), dl_eta(i + 1, j, k + 1, 0));

//                     // =====================================================
//                     // 4 zeta-edges
//                     // =====================================================
//                     unit_t_zeta(i, j, k, tx, ty, tz);
//                     push(tx, ty, tz, Jzeta(i, j, k, 0), dl_zeta(i, j, k, 0));

//                     unit_t_zeta(i + 1, j, k, tx, ty, tz);
//                     push(tx, ty, tz, Jzeta(i + 1, j, k, 0), dl_zeta(i + 1, j, k, 0));

//                     unit_t_zeta(i, j + 1, k, tx, ty, tz);
//                     push(tx, ty, tz, Jzeta(i, j + 1, k, 0), dl_zeta(i, j + 1, k, 0));

//                     unit_t_zeta(i + 1, j + 1, k, tx, ty, tz);
//                     push(tx, ty, tz, Jzeta(i + 1, j + 1, k, 0), dl_zeta(i + 1, j + 1, k, 0));

//                     // =====================================================
//                     // Weighted least squares:
//                     //   minimize sum w | t·Jcell - J_edge/|dl| |^2
//                     // =====================================================
//                     double N00 = 0.0, N01 = 0.0, N02 = 0.0;
//                     double N11 = 0.0, N12 = 0.0, N22 = 0.0;
//                     double r0 = 0.0, r1 = 0.0, r2 = 0.0;

//                     for (int n = 0; n < K; ++n)
//                     {
//                         const double w = eqs[n].w;
//                         const double tx = eqs[n].tx;
//                         const double ty = eqs[n].ty;
//                         const double tz = eqs[n].tz;
//                         const double b = eqs[n].rhs;

//                         N00 += w * tx * tx;
//                         N01 += w * tx * ty;
//                         N02 += w * tx * tz;
//                         N11 += w * ty * ty;
//                         N12 += w * ty * tz;
//                         N22 += w * tz * tz;

//                         r0 += w * tx * b;
//                         r1 += w * ty * b;
//                         r2 += w * tz * b;
//                     }

//                     auto det3 = [&](double a, double b, double c,
//                                     double d, double e, double f) -> double
//                     {
//                         // | a b c |
//                         // | b d e |
//                         // | c e f |
//                         return a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d);
//                     };

//                     double det = det3(N00, N01, N02, N11, N12, N22);
//                     const double reg = 1e-14 * (N00 + N11 + N22 + 1.0);

//                     if (std::abs(det) < reg)
//                     {
//                         N00 += reg;
//                         N11 += reg;
//                         N22 += reg;
//                         det = det3(N00, N01, N02, N11, N12, N22);
//                     }

//                     // inverse of symmetric 3x3 normal matrix
//                     const double C00 = (N11 * N22 - N12 * N12);
//                     const double C01 = (N02 * N12 - N01 * N22);
//                     const double C02 = (N01 * N12 - N02 * N11);
//                     const double C11 = (N00 * N22 - N02 * N02);
//                     const double C12 = (N01 * N02 - N00 * N12);
//                     const double C22 = (N00 * N11 - N01 * N01);

//                     const double invdet = 1.0 / det;

//                     const double Jx =
//                         invdet * (C00 * r0 + C01 * r1 + C02 * r2);
//                     const double Jy =
//                         invdet * (C01 * r0 + C11 * r1 + C12 * r2);
//                     const double Jz =
//                         invdet * (C02 * r0 + C12 * r1 + C22 * r2);

//                     Jcell(i, j, k, 0) = Jx;
//                     Jcell(i, j, k, 1) = Jy;
//                     Jcell(i, j, k, 2) = Jz;
//                 }
//     }
//     cb_.sync_dJcell();
// }

void ImplicitHallSolver::BuildLinearHallCellEMF_()
{
    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &dJc = fld_->field(fid_.fid_dJcell, ib);

        // 下面这两个 field id 请替换成你实际 frozen scratch 的名字
        auto &Bflat = (*hall_face_scratch_)[ib].Bflat;      // 3-comp cell
        auto &alpha = (*hall_face_scratch_)[ib].alpha_flat; // scalar cell

        // 下面这个也是你的线性 Hall cell EMF scratch
        auto &dEhc = (*hall_face_scratch_)[ib].dEhc; // 3-comp cell

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

        auto &dEhc = scratch.dEhc;  // 3-comp cell
        auto &Pxi = scratch.P_xi;   // 6-comp xi-face projector
        auto &Pet = scratch.P_eta;  // 6-comp eta-face projector
        auto &Pze = scratch.P_zeta; // 6-comp zeta-face projector

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

// void ImplicitHallSolver::BuildLinearHallFaceEMF_()
// {
//     constexpr double eps = 1e-14;

//     auto norm3 = [](double x, double y, double z) -> double
//     {
//         return std::sqrt(x * x + y * y + z * z);
//     };

//     ClearFaceTriplet_(fid_.fid_Eface);

//     const int nb = fld_->num_blocks();
//     for (int ib = 0; ib < nb; ++ib)
//     {
//         auto &dEhc = (*hall_face_scratch_)[ib].dEhc; // 3-comp cell

//         auto &Efxi = fld_->field(fid_.fid_Eface.xi, ib); // 3-comp face
//         auto &Efet = fld_->field(fid_.fid_Eface.eta, ib);
//         auto &Efze = fld_->field(fid_.fid_Eface.zeta, ib);

//         auto &JDxi = fld_->field(fid_.fid_metric.xi, ib); // face normal vectors
//         auto &JDet = fld_->field(fid_.fid_metric.eta, ib);
//         auto &JDze = fld_->field(fid_.fid_metric.zeta, ib);

//         if (!Efxi.is_allocated() || !Efet.is_allocated() || !Efze.is_allocated())
//             continue;

//         // ============================================================
//         // xi-face : average adjacent cells, then remove normal component
//         // ============================================================
//         {
//             Int3 lo = Efxi.inner_lo();
//             Int3 hi = Efxi.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         const int iL = i - 1;
//                         const int iR = i;

//                         double Ex = 0.5 * (dEhc(iL, j, k, 0) + dEhc(iR, j, k, 0));
//                         double Ey = 0.5 * (dEhc(iL, j, k, 1) + dEhc(iR, j, k, 1));
//                         double Ez = 0.5 * (dEhc(iL, j, k, 2) + dEhc(iR, j, k, 2));

//                         const double Sx = JDxi(i, j, k, 0);
//                         const double Sy = JDxi(i, j, k, 1);
//                         const double Sz = JDxi(i, j, k, 2);

//                         const double Smag = norm3(Sx, Sy, Sz) + eps;
//                         const double nx = Sx / Smag;
//                         const double ny = Sy / Smag;
//                         const double nz = Sz / Smag;

//                         const double En = Ex * nx + Ey * ny + Ez * nz;
//                         Ex -= En * nx;
//                         Ey -= En * ny;
//                         Ez -= En * nz;

//                         Efxi(i, j, k, 0) = Ex;
//                         Efxi(i, j, k, 1) = Ey;
//                         Efxi(i, j, k, 2) = Ez;
//                     }
//         }

//         // ============================================================
//         // eta-face
//         // ============================================================
//         {
//             Int3 lo = Efet.inner_lo();
//             Int3 hi = Efet.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         const int jL = j - 1;
//                         const int jR = j;

//                         double Ex = 0.5 * (dEhc(i, jL, k, 0) + dEhc(i, jR, k, 0));
//                         double Ey = 0.5 * (dEhc(i, jL, k, 1) + dEhc(i, jR, k, 1));
//                         double Ez = 0.5 * (dEhc(i, jL, k, 2) + dEhc(i, jR, k, 2));

//                         const double Sx = JDet(i, j, k, 0);
//                         const double Sy = JDet(i, j, k, 1);
//                         const double Sz = JDet(i, j, k, 2);

//                         const double Smag = norm3(Sx, Sy, Sz) + eps;
//                         const double nx = Sx / Smag;
//                         const double ny = Sy / Smag;
//                         const double nz = Sz / Smag;

//                         const double En = Ex * nx + Ey * ny + Ez * nz;
//                         Ex -= En * nx;
//                         Ey -= En * ny;
//                         Ez -= En * nz;

//                         Efet(i, j, k, 0) = Ex;
//                         Efet(i, j, k, 1) = Ey;
//                         Efet(i, j, k, 2) = Ez;
//                     }
//         }

//         // ============================================================
//         // zeta-face
//         // ============================================================
//         {
//             Int3 lo = Efze.inner_lo();
//             Int3 hi = Efze.inner_hi();

//             for (int i = lo.i; i < hi.i; ++i)
//                 for (int j = lo.j; j < hi.j; ++j)
//                     for (int k = lo.k; k < hi.k; ++k)
//                     {
//                         const int kL = k - 1;
//                         const int kR = k;

//                         double Ex = 0.5 * (dEhc(i, j, kL, 0) + dEhc(i, j, kR, 0));
//                         double Ey = 0.5 * (dEhc(i, j, kL, 1) + dEhc(i, j, kR, 1));
//                         double Ez = 0.5 * (dEhc(i, j, kL, 2) + dEhc(i, j, kR, 2));

//                         const double Sx = JDze(i, j, k, 0);
//                         const double Sy = JDze(i, j, k, 1);
//                         const double Sz = JDze(i, j, k, 2);

//                         const double Smag = norm3(Sx, Sy, Sz) + eps;
//                         const double nx = Sx / Smag;
//                         const double ny = Sy / Smag;
//                         const double nz = Sz / Smag;

//                         const double En = Ex * nx + Ey * ny + Ez * nz;
//                         Ex -= En * nx;
//                         Ey -= En * ny;
//                         Ez -= En * nz;

//                         Efze(i, j, k, 0) = Ex;
//                         Efze(i, j, k, 1) = Ey;
//                         Efze(i, j, k, 2) = Ez;
//                     }
//         }
//     }

//     cb_.sync_dEface();
// }

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

// Jacobi

PetscErrorCode ImplicitHallSolver::FormJacobian_(SNES, Vec X, Mat, Mat, void *ctx)
{
    auto *S = static_cast<ImplicitHallSolver *>(ctx);
    // 最简单版本：直接沿用当前 time-step 冻结状态
    // 如果还没准备，就准备一次
    if (!S->p0_frozen_ready_)
        S->PrepareWhistlerP0FrozenState_();
    return 0;

    // auto *S = static_cast<ImplicitHallSolver *>(ctx);
    // // 用当前 X 形成 trial E / trial B
    // S->UnpackVecToEhallField_(X);
    // S->BuildTrialBfaceFromUnknownE_();
    // // 基于当前 trial B 冻结线性化状态
    // S->cb_.calc_Bcell_from_current_Bface();
    // S->cb_.FillFrozenBflatFromCurrentBcell_();
    // S->cb_.FillFrozenAlphaFlatCell_();
    // S->p0_frozen_ready_ = true;
    // return 0;
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
    // 先直接把你现有的 P0 operator 拿来当 Jacobian-vector
    return ApplyWhistlerP0Operator_(in, out);
}
#endif
