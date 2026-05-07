#pragma once

#include "00_Mercury_Const.h"

#if HALL_IMPLICIT == 1

#include <functional>
#include <vector>
#include <stdexcept>
#include <petscsnes.h>

#include "0_SolverFields.h"
#include "2_topology/2_MPCNS_Topology_Equiv.h"
#include "4_halo/1_MPCNS_Halo_EdgeOwner.h"
#include "operators/CTOperators.h"
#include "4_Hall_Implicit_Type.h"

// forward decl
class Grid;
namespace TOPO
{
    class Topology;
}
class Field;
class Halo;
class Param;
class MercuryBoundary;
class FieldBlock;

class ImplicitHallSolver
{
public:
    struct Callbacks
    {
        // 把当前 Bface 补齐到可用于 stencil / curl(B)
        std::function<void()> sync_Bface;

        // 把当前 Eedge 补齐到可用于 curl(E)
        std::function<void()> sync_Ehalledge;

        // 从当前 Bface 计算派生量
        std::function<void()> calc_PV;
        std::function<void()> calc_Uplus;

        // 当前 Bface -> Jedge/Jcell -> Ehall
        // 约定：执行后，fid_.fid_Ehall 中存放预测的 Hall edge EMF
        std::function<void()> build_Ehall_from_current_B;

        std::function<void()> calc_Bcell_from_current_Bface;
        std::function<void()> FillFrozenBflatFromCurrentBcell_;
        std::function<void()> FillFrozenAlphaFlatCell_;

        std::function<void()> sync_dEedge;
        std::function<void()> sync_dBface;
        std::function<void()> sync_dJedge;
        std::function<void()> sync_dJcell;
        std::function<void()> sync_dEface;
    };

public:
    ImplicitHallSolver() = default;
    ~ImplicitHallSolver();

    void Setup(Grid *grd,
               TOPO::Topology *topo,
               Field *fld,
               Halo *halo,
               Param *par,
               MercuryBoundary *bound,
               const SolverFields &fid,
               const TOPO::TopologyEquiv &equiv,
               const HALO_OWNER::EdgeOwnerSyncPattern &owner_pat,
               std::vector<HallFaceScratchBlock_> *hall_face_scratch);

    void SetCallbacks(const Callbacks &cb) { cb_ = cb; }

    void SetTheta(double theta) { theta_ = theta; } // 1.0: BE, 0.5: midpoint
    void SetVerbose(bool x) { verbose_ = x; }

    void InitializePetsc();
    void FinalizePetsc();

    // 解 Hall 子步；输入/输出都作用在当前 fld_->B_xi/eta/zeta 上
    void SolveOneStep(double dt, bool if_outres);

private:
    Grid *grd_{nullptr};
    TOPO::Topology *topo_{nullptr};
    Field *fld_{nullptr};
    Halo *halo_{nullptr};
    Param *par_{nullptr};
    MercuryBoundary *bound_{nullptr};

    SolverFields fid_;
    TOPO::TopologyEquiv equiv_;
    HALO_OWNER::EdgeOwnerSyncPattern owner_pat_;
    Callbacks cb_;

    bool petsc_ready_{false};
    bool verbose_{true};
    double dt_{0.0};
    double theta_{0.5};

    // PETSc
    SNES snes_{nullptr};
    Vec X_{nullptr}, F_{nullptr};
    Mat Jmf_{nullptr};
    Mat Jshell_{nullptr};

    // local owner-edge buffers
    std::vector<double> x_local_;
    std::vector<double> eh_pred_local_;

    // B* snapshot + RHS scratch
    std::vector<Scalar> Bstar_xi_;
    std::vector<Scalar> Bstar_eta_;
    std::vector<Scalar> Bstar_ze_;

    std::vector<TOPO::EdgeLocalID> owner_edges_sorted_;

private:
    static PetscErrorCode FormFunction_(SNES snes, Vec X, Vec F, void *ctx);

    void CheckReady_() const;
    void CreatePetscObjects_();
    void DestroyPetscObjects_();

    void SnapshotCurrentBface_();
    void RestoreCurrentBfaceFromSnapshot_();
    void BuildTrialBfaceFromUnknownE_();

    void UnpackVecToEhallField_(Vec X);
    void PackPredictedEhallToLocal_();

    void CopyEhallToE_();
    void ClearEdgeTriplet_(const IdTriplet &fid_triplet);
    void ClearFaceTriplet_(const IdTriplet &fid_triplet);

    void EvaluatePredictedEhallFromTrialB_();
    void WriteResidual_(Vec X, Vec F);

    // helpers
    void PackFaceInner_(int fid, std::vector<Scalar> &buf);
    void RestoreFaceInner_(int fid, std::vector<Scalar> &buf);
    void AddFaceInnerFromRHS_(int fid_B, int fid_RHS, double factor);

    std::vector<HallFaceScratchBlock_> *hall_face_scratch_;
    PetscErrorCode ApplyWhistlerP0ApproxInverse_(Vec in, Vec out);
    PetscErrorCode ApplyWhistlerP0Operator_(Vec in, Vec out);
    static PetscErrorCode PCApplyWhistlerP0_Shell_(PC pc, Vec in, Vec out);

    void PrepareWhistlerP0FrozenState_()
    {
        // 1) 恢复 snapshot B* 到当前 Bface
        RestoreCurrentBfaceFromSnapshot_();
        cb_.sync_Bface();

        // 2) 从当前 Bface 重新构造当前 Bcell / Bindcell
        //    这一步要复用你真实路径里已有的 Bface -> Bcell 逻辑
        cb_.calc_Bcell_from_current_Bface(); // 名字你自己对应真实代码

        // 3) 把当前 Bcell 拷到 scratch.Bflat
        cb_.FillFrozenBflatFromCurrentBcell_();

        // 4) 用当前 UH / UNa 计算 alpha_flat(cell)
        cb_.FillFrozenAlphaFlatCell_();

        p0_frozen_ready_ = true;
    }

    void WriteP0Action_(Vec in, Vec out);
    void BuildDeltaBFromDeltaE_(Vec in);
    void Calc_DeltaJ_Edge_FromDeltaB_();
    void Calc_DeltaJcell_FromDeltaJedge_Frozen_();
    void BuildLinearHallCellEMF_();
    void BuildLinearHallFaceEMF_();
    void AssembleLinearHallEdgeEMF_();
    void SubtractPackedTempDEpreFromVec_(Vec out);

    void UnpackVecToTempDEdgeField_(Vec X);

    // Jacobi
    static PetscErrorCode FormJacobian_(SNES snes, Vec X, Mat A, Mat P, void *ctx);
    static PetscErrorCode MatMult_WhistlerShell_(Mat A, Vec in, Vec out);
    PetscErrorCode ApplyWhistlerJv_(Vec in, Vec out);

private:
    // 非 owning，只是为了方便，不单独 destroy
    KSP ksp_{nullptr};
    PC pc_{nullptr};
    bool use_shell_pc_{true};

    bool p0_frozen_ready_ = false;
    // Richardson 参数（先从最保守的开始）
    int p0_nsweeps_ = 1;    // 先做 1 步
    double p0_omega_ = 0.8; // 先取 0.8
    // PCShell 内部工作向量
    Vec pc_q_ = nullptr;   // q = P0 z
    Vec pc_res_ = nullptr; // residual = r - q

public:
    double MaxAbsTriplet_(const IdTriplet &fid_triplet);
};

#endif
