#pragma once

#include "00_Mercury_Const.h"
#include "6_io/IOModule.h"

#include "1_Boundary.h"
#include "0_SolverFields.h"
#include "2_Initial.h"
#include "3_Control.h"
#include "MercuryPhysicsTypes.h"
#include "4_halo/HaloEdgeOwner.h"
#include "7_metric/SingularEdgeRegistry.h"

#include <petscksp.h>

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

    // Topology is immutable after solver setup.  Reuse the sorted shared-edge
    // view and collective work arrays in every explicit/implicit alias
    // reduction instead of rebuilding and allocating them in each RHS call.
    bool edge_alias_reduce_cache_ready_{false};
    std::vector<const TOPO::EquivClass *> edge_alias_reduce_classes_;
    struct EdgeAliasReduceLocalTerm_
    {
        std::size_t class_index=0;
        TOPO::EntityKey edge;
        int orient_sign=1;
    };
    struct EdgeAliasReduceOwnerWrite_
    {
        std::size_t class_index=0;
        TOPO::EntityKey edge;
        int orient_sign=1;
        double member_count=1.0;
    };
    std::vector<EdgeAliasReduceLocalTerm_> edge_alias_reduce_local_terms_;
    std::vector<EdgeAliasReduceOwnerWrite_> edge_alias_reduce_owner_writes_;
    std::vector<double> edge_alias_reduce_local_sum_;
    std::vector<double> edge_alias_reduce_global_sum_;

    KSP implicit_resistive_ksp_{nullptr};
    Mat implicit_resistive_A_{nullptr};
    Vec implicit_resistive_x_{nullptr};
    Vec implicit_resistive_b_{nullptr};
    double implicit_resistive_dt_{0.0};
    bool implicit_resistive_ready_{false};
    bool implicit_resistive_has_guess_{false};

    std::vector<HallFaceScratchBlock> hall_face_scratch_;
    void SetupHallFaceScratch_();
    void SetupCellReconstructionWeights_();

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
    double MercuryResistivityShape_(double radius) const;
    void AddArtificialResistivityToEdgeEMF_();
    void AddLocalArtificialResistivityToEdgeEMF_();
    void AddIdealEdgeEMF_();
    void AddHallEdgeEMF_();
    void AddAmbipolarEdgeEMF_();
    void Calc_J_Edge();
    void AssembleSingularEdgeCurrent_(const IdTriplet &fid_Bface,
                                      const IdTriplet &fid_Jedge);
    void ReduceEdgeAliasCandidatesToOwners_(const IdTriplet &fid_edge);

    // 只更新 Bface: Bface += dt_sub * RHS_b
    void ApplyUpdate_Euler_BfaceOnly_(double dt_sub, const IdTriplet &fid_RHSB);

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
    void ApplyStationaryWallNonResistiveEMF_();
    //---------------------------------------------------------------
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
