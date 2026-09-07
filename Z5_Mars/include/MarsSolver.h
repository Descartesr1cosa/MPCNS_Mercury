#pragma once

#include "MarsAtmosphereModel.h"
#include "MarsChemistry.h"
#include "MarsSpecies.h"
#include "5_boundary/Boundary.h"
#include "6_io/IOModule.h"
#include "6_io/RunControl.h"
#include "8_dec/CTOperators.h"
#include "7_metric/SingularEdgeRegistry.h"
#include "2_topology/Topology.h"

#include <array>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

class Grid;
class Field;
class Halo;
class Param;
namespace TOPO { struct Topology; }

class MarsSolver
{
public:
    MarsSolver(Grid *grid, TOPO::Topology *topology, Field *field, Halo *halo, Param *param,
               METRIC::SingularEdgeRegistry *singular_edges = nullptr);
    static void RegisterFields(Field *field, int nghost);
    static void RegisterHaloFields(Field *field, Halo *halo);
    void Advance();

private:
    struct Triplet { int xi = -1, eta = -1, zeta = -1; };

    Grid *grid_ = nullptr;
    TOPO::Topology *topology_ = nullptr;
    Field *field_ = nullptr;
    Halo *halo_ = nullptr;
    Param *param_ = nullptr;
    METRIC::SingularEdgeRegistry *singular_edges_ = nullptr;
    BoundaryCore boundary_;
    MARS::AtmosphereModel atmosphere_;
    IOModule io_;
    Control control_;

    std::array<int, MARS::NumSpecies> U_{};
    std::array<int, MARS::NumSpecies> PV_{};
    std::array<int, MARS::NumSpecies> RHS_{};
    std::vector<std::string> ucell_names_, bface_names_, uplus_names_, bcell_names_,
                             jcell_names_, ehall_names_, eedge_names_;
    int neutral_ = -1, photo_ = -1, uplus_ = -1, bcell_ = -1, jcell_ = -1, divb_ = -1;
    int hall_ecell_ = -1, hall_beta_ = -1;
    int jac_ = -1;
    Triplet metric_, face_area_, edge_dr_, edge_dl_, B_, E_, Eface_, Ehall_, RHSB_, flux_;

    double gamma_ = 1.4;
    double cfl_ = 0.5;
    double length_ref_ = 3.396e6;
    double velocity_ref_ = 1.0;
    double density_ref_ = 1.0;
    double number_density_ref_ = 1.0;
    double pressure_ref_ = 1.0;
    double magnetic_ref_ = 1.0;
    double temperature_ref_ = 1.0;
    double boltzmann_ = 1.3807e-23;
    double charge_ = 1.602e-19;
    double mu0_ = 1.256637061e-6;
    double dt_ = 0.0;
    double dt_fluid_ = 0.0;
    double dt_mhd_ = 0.0;
    double hall_scale_ = 0.0;
    double pressure_scale_ = 0.0;
    double magnetic_force_scale_ = 0.0;
    double electron_density_floor_nd_ = 1.0e-2;
    bool hall_enabled_ = true;
    bool ambipolar_enabled_ = true;
    double solar_wind_h_density_ = 0.0;
    double solar_wind_pressure_ = 0.0;
    std::array<double, 3> solar_wind_velocity_{};
    std::array<double, 3> solar_wind_imf_{};
    bool repair_restart_farfield_bface_ = false;
    double wall_temperature_K_ = 1500.0;
    std::map<std::array<int, 4>, std::array<double, MARS::NumSpecies>> wall_ion_cache_;

    struct ElectronDensityInfo
    {
        std::array<double, MARS::NumSpecies> species_nd{};
        double ne_true = 0.0;
        double ne_eff = 0.0;
        double mhd_taper = 0.0;
    };

    bool edge_alias_reduce_cache_ready_ = false;
    std::vector<const TOPO::EquivClass *> edge_alias_reduce_classes_;
    struct EdgeAliasTerm
    {
        std::size_t class_index = 0;
        TOPO::EntityKey edge;
        int orientation = 1;
    };
    struct EdgeAliasOwnerWrite
    {
        std::size_t class_index = 0;
        TOPO::EntityKey edge;
        int orientation = 1;
        double member_count = 1.0;
    };
    std::vector<EdgeAliasTerm> edge_alias_local_terms_;
    std::vector<EdgeAliasOwnerWrite> edge_alias_owner_writes_;
    std::vector<double> edge_alias_local_sum_, edge_alias_global_sum_;
    bool overridden_edge_cache_ready_ = false;
    int overridden_edge_rank_ = 0;
    std::unordered_set<TOPO::EntityKey, TOPO::EntityKey::Hash> overridden_edges_;
    std::unordered_set<TOPO::EntityKey, TOPO::EntityKey::Hash> hall_pole_edges_;

    void CacheFields();
    void ReadConstants();
    void SetupBoundary();
    void ApplyIonFarfield(FieldBlock &block, const BOUND::PhysicalRegion &region,
                          int nghost, std::size_t species);
    void ApplyMagneticFarfieldFace(FieldBlock &block,
                                   const BOUND::PhysicalRegion &region, int nghost);
    void ApplyMagneticFarfieldCell(FieldBlock &block,
                                   const BOUND::PhysicalRegion &region, int nghost);
    void ApplyFarfieldElectricEdge(FieldBlock &block,
                                   const BOUND::PhysicalRegion &region, int nghost);
    double CheckedCellJacobian(const FieldBlock &jac, int block, int i, int j, int k,
                               const char *operation) const;
    double CheckedPositiveGeometry(const FieldBlock &geometry, int block,
                                   int i, int j, int k, const char *operation) const;
    ElectronDensityInfo ElectronDensityAtCell(int block, int i, int j, int k) const;
    void SyncGroup(const std::string &name, const std::vector<std::string> &fields,
                   bool owner_sync = false);
    void Initialize();
    void PopulateAnalyticAtmosphere();
    void ApplyIonWall(FieldBlock &block, const BOUND::PhysicalRegion &region,
                      int nghost, std::size_t species);
    void UpdatePrimitive();
    void UpdateIonVelocity();
    void UpdateMagneticCell();
    void UpdateCurrentCell();
    void ComputeTimestep();
    void ZeroRHS();
    void AssembleFluidFlux();
    void AddFortranSources();
    void ApplyPointImplicitChemistry();
    void AssembleElectricAndCT();
    void ZeroEdgeElectricFields();
    void AssembleIdealEdgeEMF();
    void AssembleIdealUCTFromFaces();
    void CacheOverriddenNonResistiveEdges();
    bool NonResistiveEdgeIsOverridden(int block, int axis, int i, int j, int k);
    bool HallPoleEdge(int block, int axis, int i, int j, int k);
    void AssembleAmbipolarEdgeEMF();
    void AssembleHallRusanovEdgeEMF();
    void ReduceEdgeAliasCandidatesToOwners(const Triplet &edge_fields);
    void AssembleSingularEdgeEMFNonHall();
    void AssembleSingularEdgeEMFHall();
    void ApplyStationaryWallEMF();
    void ZeroPoleHallEdges(int block);
    void ApplyUpdate();
    void ApplyFloors();
    void UpdateDivB();
    void PrintDiagnostics();
    bool UpdateControlAndOutput();
};
