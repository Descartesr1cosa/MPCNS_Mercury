#pragma once

#include "MarsAtmosphereModel.h"
#include "MarsChemistry.h"
#include "MarsSpecies.h"
#include "3_Control.h"
#include "5_boundary/Boundary.h"
#include "6_io/IOModule.h"
#include "operators/CTOperators.h"

#include <array>
#include <map>
#include <string>
#include <vector>

class Grid;
class Field;
class Halo;
class Param;
namespace TOPO { struct Topology; }

class MarsSolver
{
public:
    MarsSolver(Grid *grid, TOPO::Topology *topology, Field *field, Halo *halo, Param *param);
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
    BoundaryCore boundary_;
    MARS::AtmosphereModel atmosphere_;
    IOModule io_;
    Control control_;

    std::array<int, MARS::NumSpecies> U_{};
    std::array<int, MARS::NumSpecies> PV_{};
    std::array<int, MARS::NumSpecies> RHS_{};
    int neutral_ = -1, photo_ = -1, uplus_ = -1, bcell_ = -1, jcell_ = -1, divb_ = -1;
    int jac_ = -1;
    Triplet metric_, edge_dr_, B_, E_, RHSB_, flux_;

    double gamma_ = 1.4;
    double cfl_ = 0.5;
    double length_ref_ = 3.396e6;
    double velocity_ref_ = 1.0;
    double density_ref_ = 1.0;
    double pressure_ref_ = 1.0;
    double magnetic_ref_ = 1.0;
    double temperature_ref_ = 1.0;
    double boltzmann_ = 1.3807e-23;
    double charge_ = 1.602e-19;
    double mu0_ = 1.256637061e-6;
    double dt_ = 0.0;
    double hall_scale_ = 0.0;
    double pressure_scale_ = 0.0;
    bool hall_enabled_ = true;
    bool ambipolar_enabled_ = true;
    double solar_wind_h_density_ = 0.0;
    double wall_temperature_K_ = 1500.0;
    std::map<std::array<int, 4>, std::array<double, MARS::NumSpecies>> wall_ion_cache_;

    void CacheFields();
    void ReadConstants();
    void SetupBoundary();
    double CheckedCellJacobian(const FieldBlock &jac, int block, int i, int j, int k,
                               const char *operation) const;
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
    void AssembleElectricAndCT();
    void ApplyUpdate();
    void ApplyFloors();
    void PrintDiagnostics() const;
};
