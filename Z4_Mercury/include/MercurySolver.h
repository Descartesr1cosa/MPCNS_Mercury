#pragma once

#include "00_Mercury_Const.h"
#include "6_io/IOModule.h"

#include "1_Boundary.h"
#include "0_SolverFields.h"
#include "2_Initial.h"
#include "3_Control.h"
#include "4_halo/HaloEdgeOwner.h"
#include "4_Hall_Implicit_Type.h"
#include "7_metric/SingularEdgeRegistry.h"

#include <petscksp.h>

#if HALL_IMPLICIT == 1
#include "4_Hall_Implicit.h"
#endif

struct NumInfo
{
    // true number densities, nondimensional
    double nH_true{0.0};
    double nNa_true{0.0};
    double ne_true{0.0};

    // Hall denominator regularization
    double ne_eff{0.0};

    // composition fractions: chiH + chiNa = 1 if ne_true > tiny
    double chiH{1.0};
    double chiNa{0.0};

    // optional low-density MHD source weights:
    // wH_mhd + wNa_mhd = ne_true / ne_eff <= 1
    double wH_mhd{0.0};
    double wNa_mhd{0.0};
    double mhd_taper{0.0};
};

struct ResistiveEdgeEMFControl
{
    bool is_Mercury_resistance = false;
    bool use_implicit_mercury_resistance = false;
    int n_subcycles = 1;

    double implicit_ksp_rtol = 1.0e-8;
    double implicit_ksp_atol = 1.0e-12;
    int implicit_ksp_max_it = 200;
};

struct ArtificialResistivityControl
{
    double eta_max = 0.0;
    double J_range_start = 0.0;
    double J_range_on = 0.0;
    bool local_enabled = false;
    double local_eta_max = 0.0;
    double local_center[3] = {0.0, 0.0, 0.0};
    double local_r_decay = 1.0;
    double local_r_cutoff = 0.0;
};

struct AmbipolarEdgeEMFControl { bool enabled = false; };

// ---- forward declarations (avoid heavy includes in header) ----
class Grid;
namespace TOPO
{
    struct Topology;
}
namespace HALO_OWNER
{
    struct EdgeOwnerSyncPattern;
}
class Field;
class Halo;
class Param;
class FieldBlock;

class MercurySolver
{
public:
    MercurySolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo, Param *par,
                  TOPO::Topology *topo_equiv = nullptr,
                  HALO_OWNER::EdgeOwnerSyncPattern *edge_owner_pat = nullptr,
                  METRIC::SingularEdgeRegistry *singular_edges = nullptr);
    ~MercurySolver();

    static void RegisterFields(Field *fld, int ngg);
    static void RegisterCouplingChannels(Field *fld, const TOPO::Topology &topology, int dimension, int ngg);
    static void RegisterHaloFields(Field *fld, Halo *halo);

    void Advance();

private:
    // --- core pointers ---
    Grid *grd_{nullptr};
    TOPO::Topology *topo_{nullptr};
    Field *fld_{nullptr};
    Halo *halo_{nullptr};
    Param *par_{nullptr};

    TOPO::Topology *topo_equiv_{nullptr};
    HALO_OWNER::EdgeOwnerSyncPattern *edge_owner_pat_{nullptr};
    METRIC::SingularEdgeRegistry *singular_edges_{nullptr};
#if HALL_IMPLICIT == 1
    ImplicitHallSolver hall_implicit_;
#endif
    // --- components ---
    Control control_;
    MercuryBoundary mercury_bound_;
    IOModule io_;
    Mercury_Initial initial_;
    RunData *run_data_;
    RuntimeMonitor *runtime_data_;

    // --- cached field ids  ---
    SolverFields fid_;

    // --- constants ---
    double gamma_{0.0};
    double NA{0.0};
    double R_uni{0.0};
    double k_Boltz{0.0};
    double q_e{0.0};
    double mu0{0.0};
    double dt{0.0};
    double dt_hall{0.0};
    double dt_sub{0.0};
    double ne_hall_floor{0.0};
    double ne_hall_floor_dimensional{0.0};

    double U_ref{0.0};
    double L_ref{0.0};
    double B_ref{0.0};
    double T_ref{0.0};
    double n_ref{0.0};
    double rho_ref{0.0};
    double M_ref{0.0};
    double M_H{0.0};
    double M_Na{0.0};
    double m_H{0.0};
    double m_Na{0.0};
    double state_coeff_H{0.0};
    double state_coeff_Na{0.0};
    double CFL{0.0};
    double hall_coef{0.0};
    double ambi_coef{0.0};
    double hall_taper_r_min{0.0};
    double hall_taper_r_max{0.0};

    double momentum_induce_coeff{0.0};
    double momentum_hall_coeff{0.0};

    double inver_MA2{0.0};
    double inver_Rem{0.0};

    ResistiveEdgeEMFControl resist_control;
    ArtificialResistivityControl arti_resist_control;
    AmbipolarEdgeEMFControl ambipolar_control;

    struct ImplicitResistiveDof
    {
        TOPO::EntityKey edge;
        double eta = 0.0;
    };

    std::vector<TOPO::EntityKey> resist_owner_edges_sorted_;
    std::vector<ImplicitResistiveDof> implicit_resistive_dofs_;
    std::vector<double> implicit_resistive_local_;
    std::vector<Scalar> resist_Bstar_xi_;
    std::vector<Scalar> resist_Bstar_eta_;
    std::vector<Scalar> resist_Bstar_ze_;
    std::vector<Scalar> local_arti_eta_xi_, local_arti_eta_eta_, local_arti_eta_ze_;
    bool local_arti_eta_ready_{false};

    KSP implicit_resistive_ksp_{nullptr};
    Mat implicit_resistive_A_{nullptr};
    Vec implicit_resistive_x_{nullptr};
    Vec implicit_resistive_b_{nullptr};
    double implicit_resistive_dt_{0.0};
    bool implicit_resistive_ready_{false};

    std::vector<HallFaceScratchBlock_> hall_face_scratch_;
    void SetupHallFaceScratch_();

#if HALL_IMPLICIT == 1

    void FillFrozenBflatFromCurrentBcell_()
    {
        const int nb = fld_->num_blocks();

        for (int ib = 0; ib < nb; ++ib)
        {
            auto &Bc = fld_->field(fid_.fid_Bcell, ib);
            if (!Bc.is_allocated())
                continue;

            auto &buf = hall_face_scratch_[ib];
            const Int3 clo = buf.clo;
            const Int3 chi = buf.chi;

            for (int i = clo.i; i < chi.i; ++i)
                for (int j = clo.j; j < chi.j; ++j)
                    for (int k = clo.k; k < chi.k; ++k)
                    {
                        buf.Bflat(i, j, k, 0) = Bc(i, j, k, 0);
                        buf.Bflat(i, j, k, 1) = Bc(i, j, k, 1);
                        buf.Bflat(i, j, k, 2) = Bc(i, j, k, 2);
                    }
        }
    }
    void FillFrozenAlphaFlatCell_()
    {
        constexpr double eps = 1e-14;
        const int nb = fld_->num_blocks();

        for (int ib = 0; ib < nb; ++ib)
        {
            auto &UH = fld_->field(fid_.fid_U_H, ib);
            auto &UNa = fld_->field(fid_.fid_U_Na, ib);

            if (!UH.is_allocated() || !UNa.is_allocated())
                continue;

            auto &buf = hall_face_scratch_[ib];
            const Int3 clo = buf.clo;
            const Int3 chi = buf.chi;

            for (int i = clo.i; i < chi.i; ++i)
                for (int j = clo.j; j < chi.j; ++j)
                    for (int k = clo.k; k < chi.k; ++k)
                    {
                        // double num[3];
                        // Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0), num);
                        // const double ne = std::max(num[2], eps);

                        NumInfo num = Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0));
                        const double ne = num.ne_eff;

                        buf.alpha_flat(i, j, k) = hall_coef / ne;
                    }
        }
    }
#endif

private:
    //=========================================================================
    // TOOLS
    void calc_Bcell();
    void calc_Jcell();
    void calc_divB();
    void calc_PV();
    void calc_Uplus();
    void UpdateFluidDerivedFields_();
    void UpdateMagneticDerivedFields_();
    void UpdateDerivedFields_();
    void calc_physical_constant(Param *par);
    void PrintMinMaxDiagnostics_();
    // void Hall_Num_Limiter(double rhoH, double rhoNa, double *num);
    NumInfo Hall_Num_Limiter(double rhoH, double rhoNa);
    double HallRadialTaper_(double x, double y, double z);
    double HallRadialTaperEdge_(int ib, StaggerLocation loc, int i, int j, int k);

    //=========================================================================

    bool StepOnce();
    void Compute_Timestep();
    bool UpdateControlAndOutput();

    //=========================================================================
    void ZeroRHS_();
    void AssembleRHS_Induction_CT_();
    void ApplyUpdate_Euler_();
    //---------------------------------------------------------------
    // For Fluid
    void Scheme_U_();
    void AddSourceToRHS_Fluid();
    //---------------------------------------------------------------
    // For Magnetic
    void AddResistiveEdgeEMF_To_(const IdTriplet &fid_Etarget);
    void AddArtificialResistivityToEdgeEMF_();
    void AddLocalArtificialResistivityToEdgeEMF_();
    void AddIdealEdgeEMF_();
    void AddHallEdgeEMF_();
    void AddAmbipolarEdgeEMF_();
    void Calc_J_Edge();

    // 只更新 Bface: Bface += dt_sub * RHS_b
    void ApplyUpdate_Euler_BfaceOnly_(double dt_sub, const IdTriplet &fid_RHSB);
    void ResistiveDiffusionSubcycles_();

    void SetupImplicitResistiveDiffusion_();
    void DestroyImplicitResistiveDiffusion_();
    void BuildImplicitResistiveEdgeDofMap_();
    void SolveImplicitResistiveDiffusion_(double dt_step);
    void ApplyImplicitResistiveUpdate_(double dt_step);
    void SnapshotImplicitResistiveBstar_();
    void RestoreImplicitResistiveBstar_();
    void UnpackVecToImplicitEres_(Vec X);
    void PackImplicitJedgeToVec_(Vec v, const IdTriplet &fid_Jedge,
                                 bool multiply_eta, double x_shift, const PetscScalar *x_extra);
    void CalcImplicitDeltaJedgeFromDeltaB_();
    double ImplicitResistiveEtaAtEdge_(const TOPO::EntityKey &e) const;
    static PetscErrorCode MatMultImplicitResistive_(Mat A, Vec X, Vec Y);

    void BuildHallFaceEMF_Rusanov_diff_();
    //--------------------------------
    //  For Ideal
    void AssembleOneDirectionEMF_(int dir, FieldBlock &E_face,
                                  FieldBlock &Bxi, FieldBlock &Beta, FieldBlock &Bzeta,
                                  FieldBlock &Badd_xi, FieldBlock &Badd_eta, FieldBlock &Badd_zeta,
                                  FieldBlock &Jac,
                                  FieldBlock &JDxi, FieldBlock &JDet, FieldBlock &JDze,
                                  FieldBlock &Uplus);
    void AssembleEdgeEMF_FromFaceE_Ideal_();
    void AssembleSingularEdgeEMF_NonHall_();
    void AssembleSingularEdgeEMF_HallExplicit_();
    //---------------------------------------------------------------
    double ComputeMagEnergy_Cell_()
    {
        double E = 0.0;

        const int nb = fld_->num_blocks();
        for (int ib = 0; ib < nb; ++ib)
        {
            auto &Bc = fld_->field(fid_.fid_Bcell, ib);
            auto &Jdet = fld_->field(fid_.fid_Jac, ib); // 这里换成你的 cell volume / Jacobian 字段

            if (!Bc.is_allocated() || !Jdet.is_allocated())
                continue;

            Int3 lo = Bc.inner_lo();
            Int3 hi = Bc.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double Bx = Bc(i, j, k, 0);
                        const double By = Bc(i, j, k, 1);
                        const double Bz = Bc(i, j, k, 2);

                        const double vol = Jdet(i, j, k, 0); // 按你的存法改
                        E += 0.5 * (Bx * Bx + By * By + Bz * Bz) * vol;
                    }
        }

        double Eglob = 0.0;
        MPI_Allreduce(&E, &Eglob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return Eglob;
    }

    struct DebugItem
    {
        double val = 0.0;
        double xyz[3] = {0.0, 0.0, 0.0};

        int rank = -1;
        int blk = -1;
        int i = 0, j = 0, k = 0;

        int src_rank = -1;
        int src_blk = -1;
        int src_i = 0, src_j = 0, src_k = 0;

        int fid = -1;
        int comp = -1;

        // relation:
        //   0 = self
        //   1 = topo-equivalent
        //   2 = halo-peer
        //   3 = physical-touch
        int relation = 0;

        // patch tag:
        //   0 = self
        //   1 = face-inner
        //   2 = face-par
        //   3 = edge-inner
        //   4 = edge-par
        //   5 = vertex-inner
        //   6 = vertex-par
        //   7 = physical
        int patch_tag = 0;

        // aux[0], aux[1], aux[2], aux[3] 预留:
        //   aux[0] 可作 valid flag / physical dir / 其他整数信息
        int aux[4] = {0, 0, 0, 0};
    };

    void DebugFindExtremaInner(const std::vector<int> &fids = {},
                               const std::vector<std::string> &names = {},
                               bool print_min = true,
                               bool print_max = true);

    void DebugDumpPointFields(int query_rank,
                              int blk, int i, int j, int k,
                              const std::vector<int> &fids = {},
                              const std::vector<std::string> &names = {});

    void DebugDumpPointPartners(int query_rank,
                                int blk, int i, int j, int k,
                                const std::vector<int> &fids = {},
                                const std::vector<std::string> &names = {},
                                int ngh = 0,
                                bool include_topo = true,
                                bool include_halo = true,
                                bool include_physical = true);

    void Debug_TestJOperator_Manufactured(int test_id = -1);
    //=========================================================================
};
