#include "MarsSolver.h"

#include "MarsFieldCatalog.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
constexpr double DensityFloor = 1.0e-10;
constexpr double PressureFloor = 1.0e-10;

double dot3(double ax, double ay, double az, double bx, double by, double bz)
{
    return ax * bx + ay * by + az * bz;
}

void cross3(double ax, double ay, double az, double bx, double by, double bz,
            double &cx, double &cy, double &cz)
{
    cx = ay * bz - az * by;
    cy = az * bx - ax * bz;
    cz = ax * by - ay * bx;
}

bool invert3(const double a[9], double inv[9])
{
    const double det = a[0] * (a[4] * a[8] - a[5] * a[7]) -
                       a[1] * (a[3] * a[8] - a[5] * a[6]) +
                       a[2] * (a[3] * a[7] - a[4] * a[6]);
    if (!std::isfinite(det) || std::abs(det) < 1.0e-30)
        return false;
    const double d = 1.0 / det;
    inv[0] = (a[4] * a[8] - a[5] * a[7]) * d;
    inv[1] = (a[2] * a[7] - a[1] * a[8]) * d;
    inv[2] = (a[1] * a[5] - a[2] * a[4]) * d;
    inv[3] = (a[5] * a[6] - a[3] * a[8]) * d;
    inv[4] = (a[0] * a[8] - a[2] * a[6]) * d;
    inv[5] = (a[2] * a[3] - a[0] * a[5]) * d;
    inv[6] = (a[3] * a[7] - a[4] * a[6]) * d;
    inv[7] = (a[1] * a[6] - a[0] * a[7]) * d;
    inv[8] = (a[0] * a[4] - a[1] * a[3]) * d;
    return true;
}
}

void MarsSolver::RegisterFields(Field *field, int nghost)
{
    MARS::RegisterMarsFields(*field, nghost);
}

void MarsSolver::RegisterHaloFields(Field *field, Halo *halo)
{
    MARS::RegisterMarsHalo(*field, *halo);
}

MarsSolver::MarsSolver(Grid *grid, TOPO::Topology *topology, Field *field, Halo *halo,
                       Param *param, METRIC::SingularEdgeRegistry *singular_edges)
    : grid_(grid), topology_(topology), field_(field), halo_(halo), param_(param),
      singular_edges_(singular_edges)
{
    if (!grid_ || !topology_ || !field_ || !halo_ || !param_)
        throw std::invalid_argument("MarsSolver: null construction argument");
    CacheFields();
    ReadConstants();
    SetupBoundary();

    io_.Setup(param_, grid_, field_, static_cast<int>(MARS::NumSpecies * 5 + 3));
    io_.SetTecplotSingularNodeContext(topology_, singular_edges_);
    std::vector<std::string> restart;
    for (std::string_view name : MARS::SpeciesNames)
        restart.push_back("U_" + std::string(name));
    restart.insert(restart.end(), {"B_xi", "B_eta", "B_zeta"});
    io_.ClearRestartFields();
    io_.SetRestartFields(restart);
    io_.SetTecplotOutputMode(IOModule::TecplotMode::AllNode);
    std::vector<std::string> plot;
    for (std::string_view name : MARS::SpeciesNames)
        plot.push_back("PV_" + std::string(name));
    plot.insert(plot.end(), {"B_cell", "J_cell", "neutral_density", "photo_rate", "divB"});
    io_.SetTecplotPhysicalOutputs(plot);

    if (param_->GetBoo("continue_calc"))
    {
        io_.ReadRestartBinFile();
        io_.ReadRunDataFile();
    }
    Initialize();
    control_.Setup(param_);

    int64_t local_cells=0;
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        const FieldBlock &u=field_->field(U_[MARS::H],ib);
        if(!u.is_allocated()) continue;
        const Int3 lo=u.inner_lo(),hi=u.inner_hi();
        local_cells+=int64_t(hi.i-lo.i)*int64_t(hi.j-lo.j)*int64_t(hi.k-lo.k);
    }
    double local_cells_d=static_cast<double>(local_cells);
    double global_cells_d=0.0;
    PARALLEL::mpi_sum(&local_cells_d,&global_cells_d,1);
    io_.Runtime().Begin(io_.Run(),param_,static_cast<int64_t>(std::llround(global_cells_d)));
}

void MarsSolver::CacheFields()
{
    ucell_names_.reserve(MARS::NumSpecies);
    for (std::size_t s = 0; s < MARS::NumSpecies; ++s)
    {
        const std::string suffix(MARS::SpeciesNames[s]);
        U_[s] = field_->field_id("U_" + suffix);
        PV_[s] = field_->field_id("PV_" + suffix);
        RHS_[s] = field_->field_id("RHS_" + suffix);
        ucell_names_.push_back("U_" + suffix);
    }
    neutral_ = field_->field_id("neutral_density");
    photo_ = field_->field_id("photo_rate");
    uplus_ = field_->field_id("U_plus");
    bcell_ = field_->field_id("B_cell");
    jcell_ = field_->field_id("J_cell");
    hall_ecell_ = field_->field_id("E_hall_cell");
    hall_beta_ = field_->field_id("Hall_beta");
    divb_ = field_->field_id("divB");
    jac_ = field_->field_id("Jac");
    metric_ = {field_->field_id("JDxi"), field_->field_id("JDet"), field_->field_id("JDze")};
    face_area_ = {field_->field_id("Area_xi"), field_->field_id("Area_eta"),
                  field_->field_id("Area_zeta")};
    edge_dr_ = {field_->field_id("dr_xi"), field_->field_id("dr_eta"), field_->field_id("dr_zeta")};
    edge_dl_ = {field_->field_id("dl_xi"), field_->field_id("dl_eta"),
                field_->field_id("dl_zeta")};
    B_ = {field_->field_id("B_xi"), field_->field_id("B_eta"), field_->field_id("B_zeta")};
    E_ = {field_->field_id("E_xi"), field_->field_id("E_eta"), field_->field_id("E_zeta")};
    Eface_ = {field_->field_id("Eface_xi"), field_->field_id("Eface_eta"),
              field_->field_id("Eface_zeta")};
    Ehall_ = {field_->field_id("Ehall_xi"), field_->field_id("Ehall_eta"),
              field_->field_id("Ehall_zeta")};
    RHSB_ = {field_->field_id("RHS_B_xi"), field_->field_id("RHS_B_eta"), field_->field_id("RHS_B_zeta")};
    flux_ = {field_->field_id("F_xi"), field_->field_id("F_eta"), field_->field_id("F_zeta")};
    bface_names_={"B_xi","B_eta","B_zeta"};
    uplus_names_={"U_plus"};bcell_names_={"B_cell"};jcell_names_={"J_cell"};
    ehall_names_={"Ehall_xi","Ehall_eta","Ehall_zeta"};
    eedge_names_={"E_xi","E_eta","E_zeta"};
}

void MarsSolver::ReadConstants()
{
    const auto constants = param_->GetDou_List("constant").data;
    const auto ref = param_->GetDou_List("REF").data;
    gamma_ = constants.at("gamma");
    const double avogadro = constants.at("NA");
    mu0_ = constants.at("mu_mag");
    charge_ = constants.at("q_e");
    boltzmann_ = constants.at("R_uni") / avogadro;
    length_ref_ = ref.at("L_ref");
    velocity_ref_ = ref.at("U");
    temperature_ref_ = ref.at("T");
    magnetic_ref_ = ref.at("B_ref");
    number_density_ref_ = ref.at("n");
    density_ref_ = ref.at("Molecular_mass") / avogadro * number_density_ref_;
    pressure_ref_ = density_ref_ * velocity_ref_ * velocity_ref_;
    cfl_ = param_->GetDou("CFL");
    hall_enabled_ = !param_->HasBoo("is_Hall_Efield") || param_->GetBoo("is_Hall_Efield");
    ambipolar_enabled_ = !param_->HasBoo("is_Ambipolar_Efield") || param_->GetBoo("is_Ambipolar_Efield");
    hall_scale_ = magnetic_ref_ / (mu0_ * length_ref_ * charge_ * number_density_ref_ * velocity_ref_);
    pressure_scale_ = pressure_ref_ / (charge_ * number_density_ref_ * length_ref_ * velocity_ref_ * magnetic_ref_);
    magnetic_force_scale_ = magnetic_ref_ * magnetic_ref_ /
                            (mu0_ * density_ref_ * velocity_ref_ * velocity_ref_);
    // Identical smooth Hall/ambipolar density regularization to Mercury.
    electron_density_floor_nd_ = 1.0e-2;

    MARS::AtmosphereConfig atmosphere_config;
    atmosphere_config.mars_radius_m = param_->HasDou("mars_radius_m")
                                          ? param_->GetDou("mars_radius_m") : length_ref_;
    atmosphere_config.lower_boundary_altitude_km = param_->HasDou("mars_lower_boundary_km")
                                                        ? param_->GetDou("mars_lower_boundary_km") : 100.0;
    atmosphere_config.heliocentric_distance_au = param_->HasDou("mars_orbit_au")
                                                     ? param_->GetDou("mars_orbit_au") : 1.524;
    atmosphere_config.sun_direction = {
        param_->HasDou("mars_sun_x") ? param_->GetDou("mars_sun_x") : 1.0,
        param_->HasDou("mars_sun_y") ? param_->GetDou("mars_sun_y") : 0.0,
        param_->HasDou("mars_sun_z") ? param_->GetDou("mars_sun_z") : 0.0};
    const auto wall = param_->GetDou_List("Boundary_Solid_Surface").data;
    wall_temperature_K_ = wall.at("T");
    atmosphere_config.ion_temperature_K = wall_temperature_K_;
    atmosphere_config.electron_temperature_K = wall_temperature_K_;
    atmosphere_ = MARS::AtmosphereModel(atmosphere_config);
}

void MarsSolver::SetupBoundary()
{
    const auto names = field_->boundary_field_names();
    boundary_.SetUp(grid_, field_, topology_, param_, names);
    auto copy = [](FieldBlock &block, Field *field, const BOUND::PhysicalRegion &region, int ngh)
    { BoundaryCore::DefaultPhysicalCopy(block, field, region, ngh); };
    const char *regions[] = {"Outflow", "Farfield", "Pole", "Solid_Surface",
                             "Coupled-Solid", "Coupled-Fluid"};
    for (const std::string &name : names)
        for (const char *region : regions)
            boundary_.RegisterPhysical(name, region, copy);
    for (std::size_t s = 0; s < MARS::NumSpecies; ++s)
    {
        const std::string field_name = "U_" + std::string(MARS::SpeciesNames[s]);
        auto ion_wall = [this, s](FieldBlock &block, Field *,
                                  const BOUND::PhysicalRegion &region, int nghost)
        { ApplyIonWall(block, region, nghost, s); };
        boundary_.RegisterPhysical(field_name, "Solid_Surface", ion_wall);
        boundary_.RegisterPhysical(field_name, "Coupled-Solid", ion_wall);
    }
    boundary_.CheckPhysicalHandlers(names);
}

double MarsSolver::CheckedCellJacobian(const FieldBlock &jac, int block,
                                       int i, int j, int k,
                                       const char *operation) const
{
    const double raw = jac(i, j, k, 0);
    if (std::isfinite(raw) && std::abs(raw) > 0.0)
        return raw;

    std::ostringstream message;
    message << "MarsSolver: invalid cell geometry in " << operation
            << " rank=" << param_->GetInt("myid")
            << " block=" << block
            << " cell=(" << i << ',' << j << ',' << k << ')'
            << " jac=" << raw
            << ". A zero/non-finite Jacobian must not enter a physical operator.";
    throw std::runtime_error(message.str());
}

double MarsSolver::CheckedPositiveGeometry(const FieldBlock &geometry, int block,
                                           int i, int j, int k,
                                           const char *operation) const
{
    const double value = geometry(i, j, k, 0);
    if (std::isfinite(value) && value > 0.0)
        return value;

    std::ostringstream message;
    message << "MarsSolver: invalid positive geometry in " << operation
            << " rank=" << param_->GetInt("myid")
            << " block=" << block
            << " index=(" << i << ',' << j << ',' << k << ')'
            << " value=" << value;
    throw std::runtime_error(message.str());
}

MarsSolver::ElectronDensityInfo
MarsSolver::ElectronDensityAtCell(int block, int i, int j, int k) const
{
    ElectronDensityInfo result;
    for (std::size_t s = 0; s < MARS::NumSpecies; ++s)
    {
        const FieldBlock &u = field_->field(U_[s], block);
        const double rho = std::max(u(i, j, k, 0), 0.0);
        result.species_nd[s] = rho * density_ref_ /
                               (MARS::ParticleMass[s] * number_density_ref_);
        result.ne_true += result.species_nd[s];
    }
    result.ne_eff = std::sqrt(result.ne_true * result.ne_true +
                              electron_density_floor_nd_ * electron_density_floor_nd_);
    result.mhd_taper = result.ne_true / result.ne_eff;
    if (!std::isfinite(result.ne_true) || !std::isfinite(result.ne_eff) ||
        !(result.ne_eff > 0.0))
    {
        std::ostringstream message;
        message << "MarsSolver: invalid electron density"
                << " rank=" << param_->GetInt("myid")
                << " block=" << block
                << " cell=(" << i << ',' << j << ',' << k << ')'
                << " ne_true=" << result.ne_true
                << " ne_eff=" << result.ne_eff;
        throw std::runtime_error(message.str());
    }
    return result;
}

void MarsSolver::SyncGroup(const std::string &name, const std::vector<std::string> &fields,
                           bool owner_sync)
{
    boundary_.ApplyPhysical(fields);
    boundary_.ApplyPhysicalCornerDefault(fields);
    halo_->sync_group(name);
    if (owner_sync)
        halo_->sync_owner_alias_group(name);
}

void MarsSolver::Initialize()
{
    const bool fresh = !param_->GetBoo("continue_calc");
    const auto ref = param_->GetDou_List("REF").data;
    const auto ini = param_->GetDou_List("INITIAL").data;
    const double nref = ref.at("n");
    const double molref = ref.at("Molecular_mass");
    const double rho_sw = (ini.at("n") / nref) * (ini.at("Molecular_mass") / molref);
    solar_wind_h_density_ = rho_sw;
    const double cy = ini.at("c_y"), cz = ini.at("c_z");
    const double cx = -std::sqrt(std::max(0.0, 1.0 - cy * cy - cz * cz));
    const double ux = cx * ini.at("U") / velocity_ref_;
    const double uy = cy * ini.at("U") / velocity_ref_;
    const double uz = cz * ini.at("U") / velocity_ref_;
    const double p = ini.at("n") * boltzmann_ * ini.at("T") / pressure_ref_;
    const double b0[3] = {ini.at("Bx") / magnetic_ref_, ini.at("By") / magnetic_ref_,
                          ini.at("Bz") / magnetic_ref_};

    if (fresh)
    {
        for (int ib = 0; ib < field_->num_blocks(); ++ib)
        {
            for (std::size_t s = 0; s < MARS::NumSpecies; ++s)
            {
                FieldBlock &u = field_->field(U_[s], ib);
                if (!u.is_allocated()) continue;
                const bool solar_wind = s == MARS::H;
                const double rho = solar_wind ? rho_sw : DensityFloor;
                const double us = solar_wind ? ux : 0.0;
                const double vs = solar_wind ? uy : 0.0;
                const double ws = solar_wind ? uz : 0.0;
                const double ps = solar_wind ? p : std::max(PressureFloor, p * 1.0e-6);
                const Int3 lo = u.get_lo(), hi = u.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            u(i,j,k,0)=rho; u(i,j,k,1)=rho*us; u(i,j,k,2)=rho*vs; u(i,j,k,3)=rho*ws;
                            u(i,j,k,4)=ps/(gamma_-1.0)+0.5*rho*(us*us+vs*vs+ws*ws);
                        }
            }
            FieldBlock &neutral = field_->field(neutral_, ib);
            FieldBlock &photo = field_->field(photo_, ib);
            if (neutral.is_allocated())
            {
                const Int3 lo=neutral.get_lo(), hi=neutral.get_hi();
                for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
                {
                    for(std::size_t n=0;n<MARS::NumNeutral;++n) neutral(i,j,k,n)=0.0;
                    for(std::size_t n=0;n<MARS::NumPhoto;++n) photo(i,j,k,n)=0.0;
                }
            }
            const int bf[3]={B_.xi,B_.eta,B_.zeta}, mf[3]={metric_.xi,metric_.eta,metric_.zeta};
            for(int d=0;d<3;++d)
            {
                FieldBlock &face=field_->field(bf[d],ib), &area=field_->field(mf[d],ib);
                if(!face.is_allocated()) continue;
                const Int3 lo=face.get_lo(), hi=face.get_hi();
                for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
                    face(i,j,k,0)=dot3(b0[0],b0[1],b0[2],area(i,j,k,0),area(i,j,k,1),area(i,j,k,2));
            }
        }
    }
    PopulateAnalyticAtmosphere();
    SyncGroup("Ucell",ucell_names_);
    SyncGroup("Bface",bface_names_,true);
    UpdatePrimitive();
    UpdateIonVelocity();
    UpdateMagneticCell();
    UpdateCurrentCell();
}

void MarsSolver::PopulateAnalyticAtmosphere()
{
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &neutral=field_->field(neutral_,ib), &photo=field_->field(photo_,ib);
        if(!neutral.is_allocated()) continue;
        Block &grid_block=grid_->grids(ib);
        const Int3 lo=neutral.get_lo(), hi=neutral.get_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const std::array<double,3> position={
                grid_block.dual_x(i+1,j+1,k+1)*length_ref_,
                grid_block.dual_y(i+1,j+1,k+1)*length_ref_,
                grid_block.dual_z(i+1,j+1,k+1)*length_ref_};
            const auto state=atmosphere_.Evaluate(position);
            for(std::size_t n=0;n<MARS::NumNeutral;++n) neutral(i,j,k,n)=state.neutral_cm3[n];
            for(std::size_t n=0;n<MARS::NumPhoto;++n) photo(i,j,k,n)=state.photo_cm3_s[n];
        }
    }
}

void MarsSolver::ApplyIonWall(FieldBlock &block, const BOUND::PhysicalRegion &region,
                              int nghost, std::size_t species)
{
    if (!block.is_allocated()) return;
    const Box3 ghost=BoundaryCore::MakeGhostSlabFromInner(region.inner_slab,region.direction,nghost);
    const Box3 &inner=region.inner_slab;
    const int axis=std::abs(region.direction);
    const int sign=region.direction>0?1:-1;
    const int iref=axis==1?(sign<0?inner.lo.i:inner.hi.i-1):0;
    const int jref=axis==2?(sign<0?inner.lo.j:inner.hi.j-1):0;
    const int kref=axis==3?(sign<0?inner.lo.k:inner.hi.k-1):0;
    Block &grid_block=grid_->grids(region.this_block);
    for(int i=ghost.lo.i;i<ghost.hi.i;++i)
        for(int j=ghost.lo.j;j<ghost.hi.j;++j)
            for(int k=ghost.lo.k;k<ghost.hi.k;++k)
            {
                const int ii=axis==1?iref:i;
                const int jj=axis==2?jref:j;
                const int kk=axis==3?kref:k;
                const std::array<double,3> position={
                    grid_block.dual_x(ii+1,jj+1,kk+1)*length_ref_,
                    grid_block.dual_y(ii+1,jj+1,kk+1)*length_ref_,
                    grid_block.dual_z(ii+1,jj+1,kk+1)*length_ref_};
                const std::array<int,4> key={region.this_block,ii,jj,kk};
                auto cached=wall_ion_cache_.find(key);
                if(cached==wall_ion_cache_.end())
                    cached=wall_ion_cache_.emplace(key,atmosphere_.WallIonDensity(position).ion_cm3).first;
                double rho=cached->second[species]*1.0e6*MARS::ParticleMass[species]/density_ref_;
                if(species==MARS::H) rho=0.3*solar_wind_h_density_;
                rho=std::max(rho,DensityFloor);
                const double number_m3=rho*density_ref_/MARS::ParticleMass[species];
                const double pressure=std::max(number_m3*boltzmann_*wall_temperature_K_/pressure_ref_,PressureFloor);
                block(i,j,k,0)=rho;
                block(i,j,k,1)=0.0; block(i,j,k,2)=0.0; block(i,j,k,3)=0.0;
                block(i,j,k,4)=pressure/(gamma_-1.0);
            }
}

void MarsSolver::UpdatePrimitive()
{
    for(int ib=0;ib<field_->num_blocks();++ib) for(std::size_t s=0;s<MARS::NumSpecies;++s)
    {
        FieldBlock &u=field_->field(U_[s],ib), &pv=field_->field(PV_[s],ib);
        if(!u.is_allocated()) continue;
        const Int3 lo=u.get_lo(), hi=u.get_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double rho=std::max(u(i,j,k,0),DensityFloor);
            const double vx=u(i,j,k,1)/rho, vy=u(i,j,k,2)/rho, vz=u(i,j,k,3)/rho;
            const double kinetic=0.5*rho*(vx*vx+vy*vy+vz*vz);
            const double p=std::max((gamma_-1.0)*(u(i,j,k,4)-kinetic),PressureFloor);
            const double n_m3=rho*density_ref_/MARS::ParticleMass[s];
            pv(i,j,k,0)=vx; pv(i,j,k,1)=vy; pv(i,j,k,2)=vz; pv(i,j,k,3)=p;
            pv(i,j,k,4)=p*pressure_ref_/(std::max(n_m3,1.0)*boltzmann_);
        }
    }
}

void MarsSolver::UpdateIonVelocity()
{
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &up=field_->field(uplus_,ib);
        if(!up.is_allocated()) continue;
        std::array<FieldBlock *,MARS::NumSpecies> pv{};
        for(std::size_t s=0;s<MARS::NumSpecies;++s)pv[s]=&field_->field(PV_[s],ib);
        const Int3 lo=up.get_lo(), hi=up.get_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const ElectronDensityInfo density=ElectronDensityAtCell(ib,i,j,k);
            double v[3]={0,0,0};
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                for(int d=0;d<3;++d)
                    v[d]+=density.species_nd[s]*(*pv[s])(i,j,k,d);
            }
            const double inv=density.ne_true>1.0e-300?1.0/density.ne_true:0.0;
            for(int d=0;d<3;++d) up(i,j,k,d)=v[d]*inv;
        }
    }
    SyncGroup("Uplus",uplus_names_);
}

void MarsSolver::UpdateMagneticCell()
{
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &bc=field_->field(bcell_,ib), &bx=field_->field(B_.xi,ib), &by=field_->field(B_.eta,ib), &bz=field_->field(B_.zeta,ib);
        FieldBlock &ax=field_->field(metric_.xi,ib), &ay=field_->field(metric_.eta,ib), &az=field_->field(metric_.zeta,ib);
        if(!bc.is_allocated()) continue;
        const Int3 lo=bc.inner_lo(), hi=bc.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            double a[9]={ax(i+1,j,k,0),ax(i+1,j,k,1),ax(i+1,j,k,2),
                         ay(i,j+1,k,0),ay(i,j+1,k,1),ay(i,j+1,k,2),
                         az(i,j,k+1,0),az(i,j,k+1,1),az(i,j,k+1,2)}, inv[9];
            const double phi[3]={bx(i+1,j,k,0),by(i,j+1,k,0),bz(i,j,k+1,0)};
            if(!invert3(a,inv))
            {
                std::ostringstream message;
                message << "MarsSolver: singular/non-finite face metric matrix in UpdateMagneticCell"
                        << " rank=" << param_->GetInt("myid")
                        << " block=" << ib
                        << " cell=(" << i << ',' << j << ',' << k << ')';
                throw std::runtime_error(message.str());
            }
            for(int d=0;d<3;++d) bc(i,j,k,d)=inv[3*d]*phi[0]+inv[3*d+1]*phi[1]+inv[3*d+2]*phi[2];
        }
    }
    SyncGroup("Bcell",bcell_names_);
}

void MarsSolver::UpdateCurrentCell()
{
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &b=field_->field(bcell_,ib), &jv=field_->field(jcell_,ib), &jac=field_->field(jac_,ib);
        FieldBlock &A=field_->field(metric_.xi,ib), &B=field_->field(metric_.eta,ib), &C=field_->field(metric_.zeta,ib);
        if(!b.is_allocated()) continue;
        const Int3 lo=jv.inner_lo(), hi=jv.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double V=std::abs(CheckedCellJacobian(jac,ib,i,j,k,"UpdateCurrentCell"));
            double g[9];
            for(int c=0;c<3;++c)
            {
                const double di=0.5*(b(i+1,j,k,c)-b(i-1,j,k,c));
                const double dj=0.5*(b(i,j+1,k,c)-b(i,j-1,k,c));
                const double dk=0.5*(b(i,j,k+1,c)-b(i,j,k-1,c));
                for(int d=0;d<3;++d) g[3*c+d]=(0.5*(A(i,j,k,d)+A(i+1,j,k,d))*di+0.5*(B(i,j,k,d)+B(i,j+1,k,d))*dj+0.5*(C(i,j,k,d)+C(i,j,k+1,d))*dk)/V;
            }
            jv(i,j,k,0)=g[7]-g[5]; jv(i,j,k,1)=g[2]-g[6]; jv(i,j,k,2)=g[3]-g[1];
        }
    }
    SyncGroup("Jcell",jcell_names_);
}
