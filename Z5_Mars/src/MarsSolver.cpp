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
#include <fstream>
#include <iomanip>
#include <limits>
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

MarsSolver::MarsSolver(Grid *grid, TOPO::Topology *topology, Field *field, Halo *halo, Param *param)
    : grid_(grid), topology_(topology), field_(field), halo_(halo), param_(param)
{
    if (!grid_ || !topology_ || !field_ || !halo_ || !param_)
        throw std::invalid_argument("MarsSolver: null construction argument");
    CacheFields();
    ReadConstants();
    SetupBoundary();

    io_.Setup(param_, grid_, field_, static_cast<int>(MARS::NumSpecies * 5 + 3));
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
}

void MarsSolver::CacheFields()
{
    for (std::size_t s = 0; s < MARS::NumSpecies; ++s)
    {
        const std::string suffix(MARS::SpeciesNames[s]);
        U_[s] = field_->field_id("U_" + suffix);
        PV_[s] = field_->field_id("PV_" + suffix);
        RHS_[s] = field_->field_id("RHS_" + suffix);
    }
    neutral_ = field_->field_id("neutral_density");
    photo_ = field_->field_id("photo_rate");
    uplus_ = field_->field_id("U_plus");
    bcell_ = field_->field_id("B_cell");
    jcell_ = field_->field_id("J_cell");
    divb_ = field_->field_id("divB");
    jac_ = field_->field_id("Jac");
    metric_ = {field_->field_id("JDxi"), field_->field_id("JDet"), field_->field_id("JDze")};
    edge_dr_ = {field_->field_id("dr_xi"), field_->field_id("dr_eta"), field_->field_id("dr_zeta")};
    B_ = {field_->field_id("B_xi"), field_->field_id("B_eta"), field_->field_id("B_zeta")};
    E_ = {field_->field_id("E_xi"), field_->field_id("E_eta"), field_->field_id("E_zeta")};
    RHSB_ = {field_->field_id("RHS_B_xi"), field_->field_id("RHS_B_eta"), field_->field_id("RHS_B_zeta")};
    flux_ = {field_->field_id("F_xi"), field_->field_id("F_eta"), field_->field_id("F_zeta")};
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
    const double nref = ref.at("n");
    density_ref_ = ref.at("Molecular_mass") / avogadro * nref;
    pressure_ref_ = density_ref_ * velocity_ref_ * velocity_ref_;
    cfl_ = param_->GetDou("CFL");
    hall_enabled_ = !param_->HasBoo("is_Hall_Efield") || param_->GetBoo("is_Hall_Efield");
    ambipolar_enabled_ = !param_->HasBoo("is_Ambipolar_Efield") || param_->GetBoo("is_Ambipolar_Efield");
    hall_scale_ = magnetic_ref_ / (mu0_ * length_ref_ * charge_ * nref * velocity_ref_);
    pressure_scale_ = pressure_ref_ / (charge_ * nref * length_ref_ * velocity_ref_ * magnetic_ref_);
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
    boundary_.CheckPhysicalHandlers(names);
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
    LoadAtmosphereRankFile();
    std::vector<std::string> ufields;
    for(auto n:MARS::SpeciesNames) ufields.push_back("U_"+std::string(n));
    SyncGroup("Ucell",ufields);
    SyncGroup("Bface",{"B_xi","B_eta","B_zeta"},true);
    UpdatePrimitive();
    UpdateIonVelocity();
    UpdateMagneticCell();
    UpdateCurrentCell();
}

bool MarsSolver::LoadAtmosphereRankFile()
{
    const int rank=param_->GetInt("myid");
    char path[128];
    std::snprintf(path,sizeof(path),"./CASE/input/atmosphere_rank%06d.dat",rank);
    std::ifstream input(path);
    if(!input) return false;
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &neutral=field_->field(neutral_,ib), &photo=field_->field(photo_,ib);
        if(!neutral.is_allocated()) continue;
        const Int3 lo=neutral.inner_lo(), hi=neutral.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            for(std::size_t n=0;n<MARS::NumNeutral;++n) input>>neutral(i,j,k,n);
            for(std::size_t n=0;n<MARS::NumPhoto;++n) input>>photo(i,j,k,n);
            if(!input) throw std::runtime_error(std::string("short Mars atmosphere file: ")+path);
        }
    }
    return true;
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
        FieldBlock &up=field_->field(uplus_,ib); if(!up.is_allocated()) continue;
        const Int3 lo=up.get_lo(), hi=up.get_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            double ne=0.0,v[3]={0,0,0};
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                FieldBlock &u=field_->field(U_[s],ib), &pv=field_->field(PV_[s],ib);
                const double ns=std::max(u(i,j,k,0),0.0)/MARS::ParticleMass[s];
                ne+=ns; for(int d=0;d<3;++d) v[d]+=ns*pv(i,j,k,d);
            }
            const double inv=1.0/std::max(ne,1.0e-30);
            for(int d=0;d<3;++d) up(i,j,k,d)=v[d]*inv;
        }
    }
    SyncGroup("Uplus",{"U_plus"});
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
            if(invert3(a,inv)) for(int d=0;d<3;++d) bc(i,j,k,d)=inv[3*d]*phi[0]+inv[3*d+1]*phi[1]+inv[3*d+2]*phi[2];
            else for(int d=0;d<3;++d) bc(i,j,k,d)=0.0;
        }
    }
    SyncGroup("Bcell",{"B_cell"});
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
            const double V=std::max(std::abs(jac(i,j,k,0)),1e-30);
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
    SyncGroup("Jcell",{"J_cell"});
}
