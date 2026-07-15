#pragma once

#include "00_Lunar_Const.h"
#include "6_io/IOModule.h"

#include "1_Boundary.h"
#include "0_SolverFields.h"
#include "2_Initial.h"
#include "3_Control.h"
#include "LunarPhysicsTypes.h"
#include "4_halo/HaloEdgeOwner.h"
#include "7_metric/SingularEdgeRegistry.h"

// Compile-time density-floor fraction relative to the nondimensional inflow
// density. Override with -DZ3_LUNAR_DENSITY_FLOOR_FRACTION=<value>.
#ifndef Z3_LUNAR_DENSITY_FLOOR_FRACTION
#define Z3_LUNAR_DENSITY_FLOOR_FRACTION 1.0e-3
#endif


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

class LunarSolver
{
public:
    LunarSolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo, Param *par,
                  TOPO::Topology *topo_equiv = nullptr,
                  HALO_OWNER::EdgeOwnerSyncPattern *edge_owner_pat = nullptr,
                  METRIC::SingularEdgeRegistry *singular_edges = nullptr);
    ~LunarSolver();

    static void RegisterFields(Field *fld, int ngg);
    static void RegisterCouplingChannels(Field *fld, const TOPO::Topology &topology, int dimension, int ngg);
    static void RegisterHaloFields(Field *fld, Halo *halo);

    // Explicit, opt-in post-processing export interface.
    void WritePostStaticData(const std::string &output_directory,
                             POST::WriteOptions options = {}) const;

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
    LunarBoundary lunar_bound_;
    IOModule io_;
    Lunar_Initial initial_;
    RunData *run_data_;
    RuntimeMonitor *runtime_data_;

    // Optional post-data hooks. They are controlled by setup parameters and
    // remain disabled when those parameters are absent.
    bool post_static_output_enabled_{false};
    std::string post_output_path_{"./DATA_bin"};
    POST::WriteOptions post_write_options_;

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
    double m_H{0.0};
    double state_coeff_H{0.0};
    double density_floor_{0.0};
    double CFL{0.0};
    double hall_coef{0.0};
    double ambi_coef{0.0};
    double hall_taper_r_min{0.0};
    double hall_taper_r_max{0.0};
    bool hall_enabled_{true};
    bool consistent_m2_enabled_{true};
    std::string singular_current_mode_{"polynomial"};
    std::string singular_emf_mode_{"multisector_uct"};

    double momentum_induce_coeff{0.0};
    double momentum_hall_coeff{0.0};

    double inver_MA2{0.0};
    ArtificialResistivityControl arti_resist_control;
    AmbipolarEdgeEMFControl ambipolar_control;

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

    // A consistent cell M2 contributes independently from both cells sharing
    // a quotient face. Sum the oriented contributions by quotient-face class
    // and write the assembled value to every local representative.
    bool face_contribution_reduce_cache_ready_{false};
    std::vector<const TOPO::EquivClass *> face_contribution_reduce_classes_;
    struct FaceContributionLocalTerm_
    {
        std::size_t class_index=0;
        TOPO::EntityKey face;
        int orient_sign=1;
    };
    std::vector<FaceContributionLocalTerm_> face_contribution_local_terms_;
    std::vector<double> face_contribution_local_sum_;
    std::vector<double> face_contribution_global_sum_;

    std::vector<HallFaceScratchBlock> hall_face_scratch_;
    void SetupHallFaceScratch_();
    void SetupCellReconstructionWeights_();
    void CanonicalizeSharedFaceGeometry_();

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
    void UpdateTecplotNodeFields_();
    void calc_physical_constant(Param *par);
    void PrintMinMaxDiagnostics_();
    NumInfo Hall_Num_Limiter(double rhoH);
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
    void ApplyDensityFloor_();
    //---------------------------------------------------------------
    // For Fluid
    void Scheme_U_();
    void AddSourceToRHS_Fluid();
    //---------------------------------------------------------------
    // For Magnetic
    void AddArtificialResistivityToEdgeEMF_();
    void AddLocalArtificialResistivityToEdgeEMF_();
    void AddIdealEdgeEMF_();
    void AddHallEdgeEMF_();
    void AddAmbipolarEdgeEMF_();
    void Calc_J_Edge();
    void AssembleSingularEdgeCurrent_(const IdTriplet &fid_Bface,
                                      const IdTriplet &fid_Jedge);
    void ReduceEdgeAliasCandidatesToOwners_(const IdTriplet &fid_edge);
    void ApplyConsistentM2ToBface_();
    void ReduceFaceContributionsToOwners_(const IdTriplet &fid_face);
    void ExtrapolatePhysicalBoundaryTangentialCurrent_();

    // 只更新 Bface: Bface += dt_sub * RHS_b
    void ApplyUpdate_Euler_BfaceOnly_(double dt_sub, const IdTriplet &fid_RHSB);

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
