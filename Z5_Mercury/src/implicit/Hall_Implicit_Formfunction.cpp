#include "4_Hall_Implicit.h"

#ifdef HALL_IMPLICIT

#include "3_field/2_MPCNS_Field.h"

PetscErrorCode ImplicitHallSolver::FormFunction_(SNES, Vec X, Vec F, void *ctx)
{
    auto *S = static_cast<ImplicitHallSolver *>(ctx);

    // X -> Ehall(unknown)
    S->UnpackVecToEhallField_(X);

    // Btrial = B* + theta*dt*(-curl(Ehall))
    S->BuildTrialBfaceFromUnknownE_();

    // 当前 trial Bface -> predicted Ehall
    S->EvaluatePredictedEhallFromTrialB_();

    // F = X - Ehall_pred
    S->WriteResidual_(X, F);

    return 0;
}

void ImplicitHallSolver::UnpackVecToEhallField_(Vec X)
{
    const PetscScalar *xarr = nullptr;
    VecGetArrayRead(X, &xarr);

    for (size_t i = 0; i < x_local_.size(); ++i)
        x_local_[i] = static_cast<double>(xarr[i]);

    VecRestoreArrayRead(X, &xarr);

    HALO_OWNER::unpack_owner_edge_1form_local(
        x_local_, *fld_, fid_.fid_Ehall, equiv_, owner_edges_sorted_, owner_pat_);
}

void ImplicitHallSolver::BuildTrialBfaceFromUnknownE_()
{
    // 从 snapshot 的 B* 开始
    RestoreCurrentBfaceFromSnapshot_();

    cb_.sync_Ehalledge();

    // RHS_B = -curl(E)
    // ClearFaceTriplet_(fid_.fid_RHS_b);

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

    // midpoint: Btrial = B* + theta*dt*RHS
    AddFaceInnerFromRHS_(fid_.fid_B.xi, fid_.fid_RHS_b.xi, theta_ * dt_);
    AddFaceInnerFromRHS_(fid_.fid_B.eta, fid_.fid_RHS_b.eta, theta_ * dt_);
    AddFaceInnerFromRHS_(fid_.fid_B.zeta, fid_.fid_RHS_b.zeta, theta_ * dt_);

    cb_.sync_Bface();
}

void ImplicitHallSolver::EvaluatePredictedEhallFromTrialB_()
{
    // cb_.calc_PV();
    // cb_.calc_Uplus();

    // 约定：这个 callback 会基于“当前 fld_->B_xi/eta/zeta”
    // 重新写出 fid_.fid_Ehall
    cb_.build_Ehall_from_current_B();
}

void ImplicitHallSolver::PackPredictedEhallToLocal_()
{
    HALO_OWNER::pack_owner_edge_1form_local(
        *fld_, fid_.fid_Ehall, equiv_, owner_edges_sorted_, eh_pred_local_);
}

void ImplicitHallSolver::WriteResidual_(Vec X, Vec F)
{
    PackPredictedEhallToLocal_();

    const PetscScalar *xarr = nullptr;
    PetscScalar *farr = nullptr;

    VecGetArrayRead(X, &xarr);
    VecGetArray(F, &farr);

    for (size_t i = 0; i < eh_pred_local_.size(); ++i)
        farr[i] = xarr[i] - eh_pred_local_[i];

    VecRestoreArrayRead(X, &xarr);
    VecRestoreArray(F, &farr);
}
#endif