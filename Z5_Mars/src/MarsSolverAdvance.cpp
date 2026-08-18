#include "MarsSolver.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace
{
constexpr double DensityFloor = 1.0e-10;
constexpr double PressureFloor = 1.0e-10;

void cross3(double ax, double ay, double az, double bx, double by, double bz,
            double &cx, double &cy, double &cz)
{
    cx = ay * bz - az * by;
    cy = az * bx - ax * bz;
    cz = ax * by - ay * bx;
}
}

void MarsSolver::ComputeTimestep()
{
    double local=std::numeric_limits<double>::max();
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &jac=field_->field(jac_,ib), &ax=field_->field(metric_.xi,ib), &ay=field_->field(metric_.eta,ib), &az=field_->field(metric_.zeta,ib);
        if(!jac.is_allocated()) continue;
        const Int3 lo=jac.inner_lo(), hi=jac.inner_hi();
        for(std::size_t s=0;s<MARS::NumSpecies;++s)
        {
            FieldBlock &u=field_->field(U_[s],ib), &pv=field_->field(PV_[s],ib);
            if(!u.is_allocated()) continue;
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                const double rho=std::max(u(i,j,k,0),1e-10), p=std::max(pv(i,j,k,3),1e-10);
                const double cs=std::sqrt(gamma_*p/rho);
                auto term=[&](FieldBlock &a,int ii,int jj,int kk)
                {
                    const double x=a(ii,jj,kk,0),y=a(ii,jj,kk,1),z=a(ii,jj,kk,2);
                    const double area=std::sqrt(x*x+y*y+z*z);
                    return std::abs(pv(i,j,k,0)*x+pv(i,j,k,1)*y+pv(i,j,k,2)*z)+cs*area;
                };
                const double denom=term(ax,i,j,k)+term(ax,i+1,j,k)+term(ay,i,j,k)+term(ay,i,j+1,k)+term(az,i,j,k)+term(az,i,j,k+1);
                if(denom>0.0) local=std::min(local,cfl_*std::abs(jac(i,j,k,0))/denom);
            }
        }

        // The original Fortran source is explicit. Keep its loss time scale in
        // the global CFL decision so chemistry cannot remove more than 20% of
        // a species in one step.
        FieldBlock &neutral=field_->field(neutral_,ib), &photo=field_->field(photo_,ib);
        if(neutral.is_allocated()) for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            MARS::ChemistryInput input;
            for(std::size_t n=0;n<MARS::NumNeutral;++n) input.neutral_cm3[n]=std::max(neutral(i,j,k,n),0.0);
            for(std::size_t n=0;n<MARS::NumPhoto;++n) input.photo_cm3_s[n]=std::max(photo(i,j,k,n),0.0);
            double pe=0.0,ne=0.0;
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                FieldBlock &u=field_->field(U_[s],ib), &pv=field_->field(PV_[s],ib);
                input.ion_cm3[s]=std::max(u(i,j,k,0),0.0)*density_ref_/MARS::ParticleMass[s]*1e-6;
                input.ion_temperature_K[s]=pv(i,j,k,4);
                pe+=pv(i,j,k,3)*pressure_ref_; ne+=input.ion_cm3[s]*1e6;
            }
            input.electron_temperature_K=pe/(std::max(ne,1.0)*boltzmann_);
            const auto rates=MARS::EvaluateChemistry(input);
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
                if(rates.loss_cm3_s[s]>0.0 && input.ion_cm3[s]>0.0)
                    local=std::min(local,0.2*input.ion_cm3[s]/rates.loss_cm3_s[s]*velocity_ref_/length_ref_);
        }
    }
    PARALLEL::mpi_min(&local,&dt_,1);
    if(!std::isfinite(dt_) || dt_<=0.0)
        throw std::runtime_error("MarsSolver: non-positive/invalid timestep");
}

void MarsSolver::ZeroRHS()
{
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        for(std::size_t s=0;s<MARS::NumSpecies;++s)
        {
            FieldBlock &r=field_->field(RHS_[s],ib); if(!r.is_allocated()) continue;
            const Int3 lo=r.inner_lo(),hi=r.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k) for(int m=0;m<5;++m) r(i,j,k,m)=0.0;
        }
        for(int fid:{RHSB_.xi,RHSB_.eta,RHSB_.zeta})
        {
            FieldBlock &r=field_->field(fid,ib); if(!r.is_allocated()) continue;
            const Int3 lo=r.inner_lo(),hi=r.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k) r(i,j,k,0)=0.0;
        }
    }
}

void MarsSolver::AssembleFluidFlux()
{
    const int metrics[3]={metric_.xi,metric_.eta,metric_.zeta};
    const int fluxes[3]={flux_.xi,flux_.eta,flux_.zeta};
    for(int ib=0;ib<field_->num_blocks();++ib) for(std::size_t s=0;s<MARS::NumSpecies;++s)
    {
        FieldBlock &u=field_->field(U_[s],ib),&pv=field_->field(PV_[s],ib),&rhs=field_->field(RHS_[s],ib),&jac=field_->field(jac_,ib);
        if(!u.is_allocated()) continue;
        for(int d=0;d<3;++d)
        {
            FieldBlock &area=field_->field(metrics[d],ib),&f=field_->field(fluxes[d],ib);
            const Int3 lo=f.inner_lo(),hi=f.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                int il=i,jl=j,kl=k,ir=i,jr=j,kr=k;
                if(d==0) --il; else if(d==1) --jl; else --kl;
                double UL[5],UR[5],FL[5],FR[5];
                auto state=[&](int ii,int jj,int kk,double out[5],double flux[5],double &rad)
                {
                    const double rho=std::max(u(ii,jj,kk,0),1e-10),vx=pv(ii,jj,kk,0),vy=pv(ii,jj,kk,1),vz=pv(ii,jj,kk,2),p=std::max(pv(ii,jj,kk,3),1e-10);
                    const double ax=area(i,j,k,0),ay=area(i,j,k,1),az=area(i,j,k,2),vn=vx*ax+vy*ay+vz*az;
                    out[0]=rho;out[1]=rho*vx;out[2]=rho*vy;out[3]=rho*vz;out[4]=p/(gamma_-1.0)+0.5*rho*(vx*vx+vy*vy+vz*vz);
                    flux[0]=rho*vn;flux[1]=out[1]*vn+p*ax;flux[2]=out[2]*vn+p*ay;flux[3]=out[3]*vn+p*az;flux[4]=(out[4]+p)*vn;
                    rad=std::abs(vn)+std::sqrt(gamma_*p/rho)*std::sqrt(ax*ax+ay*ay+az*az);
                };
                double rl,rr;state(il,jl,kl,UL,FL,rl);state(ir,jr,kr,UR,FR,rr);
                for(int m=0;m<5;++m) f(i,j,k,m)=0.5*(FL[m]+FR[m])-0.5*std::max(rl,rr)*(UR[m]-UL[m]);
            }
        }
        FieldBlock &fx=field_->field(flux_.xi,ib),&fy=field_->field(flux_.eta,ib),&fz=field_->field(flux_.zeta,ib);
        const Int3 lo=rhs.inner_lo(),hi=rhs.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k) for(int m=0;m<5;++m)
            rhs(i,j,k,m)-=((fx(i+1,j,k,m)-fx(i,j,k,m))+(fy(i,j+1,k,m)-fy(i,j,k,m))+(fz(i,j,k+1,m)-fz(i,j,k,m)))/jac(i,j,k,0);
    }
}

void MarsSolver::AddFortranSources()
{
    const double time_ref=length_ref_/velocity_ref_;
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &neutral=field_->field(neutral_,ib),&photo=field_->field(photo_,ib),&up=field_->field(uplus_,ib),&b=field_->field(bcell_,ib),&jv=field_->field(jcell_,ib),&jac=field_->field(jac_,ib);
        FieldBlock &ax=field_->field(metric_.xi,ib),&ay=field_->field(metric_.eta,ib),&az=field_->field(metric_.zeta,ib);
        if(!neutral.is_allocated()) continue;
        const Int3 lo=jac.inner_lo(),hi=jac.inner_hi();
        auto electron_pressure=[&](int i,int jj,int k){double p=0;for(std::size_t s=0;s<MARS::NumSpecies;++s)p+=field_->field(PV_[s],ib)(i,jj,k,3);return p;};
        for(int i=lo.i;i<hi.i;++i) for(int jj=lo.j;jj<hi.j;++jj) for(int k=lo.k;k<hi.k;++k)
        {
            MARS::ChemistryInput input; double ne_cm3=0.0;
            for(std::size_t n=0;n<MARS::NumNeutral;++n) input.neutral_cm3[n]=std::max(neutral(i,jj,k,n),0.0);
            for(std::size_t n=0;n<MARS::NumPhoto;++n) input.photo_cm3_s[n]=std::max(photo(i,jj,k,n),0.0);
            double pe_dim=0.0;
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                FieldBlock &u=field_->field(U_[s],ib),&pv=field_->field(PV_[s],ib);
                input.ion_cm3[s]=std::max(u(i,jj,k,0),0.0)*density_ref_/MARS::ParticleMass[s]*1e-6;
                input.ion_temperature_K[s]=pv(i,jj,k,4);ne_cm3+=input.ion_cm3[s];pe_dim+=pv(i,jj,k,3)*pressure_ref_;
            }
            input.electron_temperature_K=pe_dim/(std::max(ne_cm3*1e6,1.0)*boltzmann_);
            const auto rates=MARS::EvaluateChemistry(input);

            const double V=std::max(std::abs(jac(i,jj,k,0)),1e-30);
            double gradp[3];
            const double dpi=0.5*(electron_pressure(i+1,jj,k)-electron_pressure(i-1,jj,k));
            const double dpj=0.5*(electron_pressure(i,jj+1,k)-electron_pressure(i,jj-1,k));
            const double dpk=0.5*(electron_pressure(i,jj,k+1)-electron_pressure(i,jj,k-1));
            for(int d=0;d<3;++d) gradp[d]=(0.5*(ax(i,jj,k,d)+ax(i+1,jj,k,d))*dpi+0.5*(ay(i,jj,k,d)+ay(i,jj+1,k,d))*dpj+0.5*(az(i,jj,k,d)+az(i,jj,k+1,d))*dpk)/V;
            double jxb[3];cross3(jv(i,jj,k,0),jv(i,jj,k,1),jv(i,jj,k,2),b(i,jj,k,0),b(i,jj,k,1),b(i,jj,k,2),jxb[0],jxb[1],jxb[2]);

            std::array<double, MARS::NumReactiveNeutral> collision_frequency{};
            double collision=0.0;
            for (std::size_t n = 0; n < MARS::NumReactiveNeutral; ++n)
            {
                collision_frequency[n] = input.neutral_cm3[n] * 2.0e-10;
                collision += collision_frequency[n];
            }
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                FieldBlock &u=field_->field(U_[s],ib),&pv=field_->field(PV_[s],ib),&rhs=field_->field(RHS_[s],ib);
                const double rho=std::max(u(i,jj,k,0),0.0),ns=input.ion_cm3[s],fraction=ns/std::max(ne_cm3,1e-30);
                const double mass_source=(rates.production_cm3_s[s]-rates.loss_cm3_s[s])*1e6*MARS::ParticleMass[s]/density_ref_*time_ref;
                rhs(i,jj,k,0)+=mass_source;
                double relxb[3];cross3(pv(i,jj,k,0)-up(i,jj,k,0),pv(i,jj,k,1)-up(i,jj,k,1),pv(i,jj,k,2)-up(i,jj,k,2),b(i,jj,k,0),b(i,jj,k,1),b(i,jj,k,2),relxb[0],relxb[1],relxb[2]);
                const double gyro=charge_*magnetic_ref_*length_ref_/(MARS::ParticleMass[s]*velocity_ref_);
                const double lorentz=magnetic_ref_*magnetic_ref_/(mu0_*density_ref_*velocity_ref_*velocity_ref_);
                const double loss_mass=rates.loss_cm3_s[s]*1e6*MARS::ParticleMass[s]/density_ref_*time_ref;
                double electromagnetic_force[3];
                double force[3];
                for(int d=0;d<3;++d)
                {
                    electromagnetic_force[d]=rho*gyro*relxb[d]+fraction*lorentz*jxb[d]-fraction*gradp[d];
                    force[d]=electromagnetic_force[d]-rho*time_ref*collision*pv(i,jj,k,d)-loss_mass*pv(i,jj,k,d);
                    rhs(i,jj,k,d+1)+=force[d];
                }
                const double speed2=pv(i,jj,k,0)*pv(i,jj,k,0)+pv(i,jj,k,1)*pv(i,jj,k,1)+pv(i,jj,k,2)*pv(i,jj,k,2);
                double b1=0.0,b2=0.0;
                for(std::size_t n=0;n<MARS::NumReactiveNeutral;++n)
                {
                    const double denominator=MARS::ParticleMass[s]+MARS::ParticleMass[n];
                    b1+=MARS::ParticleMass[s]*collision_frequency[n]/denominator;
                    b2+=MARS::ParticleMass[n]*collision_frequency[n]/denominator;
                }
                b1*=(185.0-std::max(pv(i,jj,k,4),50.0));
                const double a5=3.0e6*length_ref_*boltzmann_/(density_ref_*velocity_ref_*velocity_ref_*velocity_ref_);
                const double a6=1.0e6*length_ref_*boltzmann_/(density_ref_*velocity_ref_*velocity_ref_*velocity_ref_*(gamma_-1.0));
                rhs(i,jj,k,4)+=electromagnetic_force[0]*pv(i,jj,k,0)+electromagnetic_force[1]*pv(i,jj,k,1)+electromagnetic_force[2]*pv(i,jj,k,2);
                rhs(i,jj,k,4)+=a5*ns*b1+time_ref*rho*speed2*b2;
                rhs(i,jj,k,4)+=a6*(rates.production_cm3_s[s]*185.0-rates.loss_cm3_s[s]*std::max(pv(i,jj,k,4),50.0));
                rhs(i,jj,k,4)-=0.5*loss_mass*speed2;
                rhs(i,jj,k,4)+=a6*ns*(rates.electron_production_cm3_s*185.0-rates.electron_loss_cm3_s*input.electron_temperature_K)/std::max(ne_cm3,1.0e-30);
            }
        }
    }
}

void MarsSolver::AssembleElectricAndCT()
{
    const int ef[3]={E_.xi,E_.eta,E_.zeta},drf[3]={edge_dr_.xi,edge_dr_.eta,edge_dr_.zeta};
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &b=field_->field(bcell_,ib),&jv=field_->field(jcell_,ib),&up=field_->field(uplus_,ib),&jac=field_->field(jac_,ib);
        FieldBlock &ax=field_->field(metric_.xi,ib),&ay=field_->field(metric_.eta,ib),&az=field_->field(metric_.zeta,ib);
        if(!b.is_allocated()) continue;
        auto pe=[&](int i,int j,int k){double value=0;for(std::size_t s=0;s<MARS::NumSpecies;++s)value+=field_->field(PV_[s],ib)(i,j,k,3);return value;};
        auto ecell=[&](int i,int j,int k,double e[3])
        {
            double uxB[3],jxB[3];cross3(up(i,j,k,0),up(i,j,k,1),up(i,j,k,2),b(i,j,k,0),b(i,j,k,1),b(i,j,k,2),uxB[0],uxB[1],uxB[2]);
            cross3(jv(i,j,k,0),jv(i,j,k,1),jv(i,j,k,2),b(i,j,k,0),b(i,j,k,1),b(i,j,k,2),jxB[0],jxB[1],jxB[2]);
            double ne=0;for(std::size_t s=0;s<MARS::NumSpecies;++s) ne+=std::max(field_->field(U_[s],ib)(i,j,k,0),0.0)/MARS::ParticleMass[s];
            const double nref=param_->GetDou_List("REF").data.at("n"),ne_nd=ne*density_ref_/nref;
            const double V=std::max(std::abs(jac(i,j,k,0)),1e-30),dpi=0.5*(pe(i+1,j,k)-pe(i-1,j,k)),dpj=0.5*(pe(i,j+1,k)-pe(i,j-1,k)),dpk=0.5*(pe(i,j,k+1)-pe(i,j,k-1));
            for(int d=0;d<3;++d)
            {
                const double gp=(0.5*(ax(i,j,k,d)+ax(i+1,j,k,d))*dpi+0.5*(ay(i,j,k,d)+ay(i,j+1,k,d))*dpj+0.5*(az(i,j,k,d)+az(i,j,k+1,d))*dpk)/V;
                e[d]=-uxB[d]+(hall_enabled_?hall_scale_*jxB[d]/std::max(ne_nd,1e-12):0.0)-(ambipolar_enabled_?pressure_scale_*gp/std::max(ne_nd,1e-12):0.0);
            }
        };
        for(int d=0;d<3;++d)
        {
            FieldBlock &edge=field_->field(ef[d],ib),&dr=field_->field(drf[d],ib);const Int3 lo=edge.inner_lo(),hi=edge.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                double avg[3]={0,0,0};
                for(int a=0;a<2;++a) for(int c=0;c<2;++c)
                {
                    int ii=i,jj=j,kk=k;if(d==0){jj=j-1+a;kk=k-1+c;}else if(d==1){ii=i-1+a;kk=k-1+c;}else{ii=i-1+a;jj=j-1+c;}
                    double e[3];ecell(ii,jj,kk,e);for(int q=0;q<3;++q)avg[q]+=0.25*e[q];
                }
                edge(i,j,k,0)=avg[0]*dr(i,j,k,0)+avg[1]*dr(i,j,k,1)+avg[2]*dr(i,j,k,2);
            }
        }
    }
    SyncGroup("Eedge",{"E_xi","E_eta","E_zeta"},true);
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &ex=field_->field(E_.xi,ib),&ey=field_->field(E_.eta,ib),&ez=field_->field(E_.zeta,ib);
        if(!ex.is_allocated())continue;
        CTOperators::CurlEdgeToFace(ib,ex,ey,ez,field_->field(RHSB_.xi,ib),field_->field(RHSB_.eta,ib),field_->field(RHSB_.zeta,ib),-1.0);
    }
}

void MarsSolver::ApplyUpdate()
{
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        for(std::size_t s=0;s<MARS::NumSpecies;++s){FieldBlock &u=field_->field(U_[s],ib),&r=field_->field(RHS_[s],ib);if(!u.is_allocated())continue;const Int3 lo=u.inner_lo(),hi=u.inner_hi();for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)for(int m=0;m<5;++m)u(i,j,k,m)+=dt_*r(i,j,k,m);}
        for(auto pair:{std::pair<int,int>{B_.xi,RHSB_.xi},{B_.eta,RHSB_.eta},{B_.zeta,RHSB_.zeta}}){FieldBlock &u=field_->field(pair.first,ib),&r=field_->field(pair.second,ib);if(!u.is_allocated())continue;const Int3 lo=u.inner_lo(),hi=u.inner_hi();for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)u(i,j,k,0)+=dt_*r(i,j,k,0);}
    }
}

void MarsSolver::ApplyFloors()
{
    for(int ib=0;ib<field_->num_blocks();++ib)for(std::size_t s=0;s<MARS::NumSpecies;++s){FieldBlock &u=field_->field(U_[s],ib);if(!u.is_allocated())continue;const Int3 lo=u.inner_lo(),hi=u.inner_hi();for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k){double rho=u(i,j,k,0);if(!std::isfinite(rho)||rho<DensityFloor){const double scale=DensityFloor/std::max(rho,1e-30);for(int m=1;m<4;++m)u(i,j,k,m)*=scale;u(i,j,k,0)=DensityFloor;rho=DensityFloor;}const double ke=0.5*(u(i,j,k,1)*u(i,j,k,1)+u(i,j,k,2)*u(i,j,k,2)+u(i,j,k,3)*u(i,j,k,3))/rho;u(i,j,k,4)=std::max(u(i,j,k,4),ke+PressureFloor/(gamma_-1.0));}}
}

void MarsSolver::PrintDiagnostics() const
{
    double local=0.0;for(int ib=0;ib<field_->num_blocks();++ib)for(std::size_t s=0;s<MARS::NumSpecies;++s){const FieldBlock &u=field_->field(U_[s],ib);if(!u.is_allocated())continue;const Int3 lo=u.inner_lo(),hi=u.inner_hi();for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)local=std::max(local,std::abs(u(i,j,k,0)));}double global=0;PARALLEL::mpi_max(&local,&global,1);if(param_->GetInt("myid")==0)std::printf("[MarsZ5] step=%lld time=%.8e dt=%.3e max(rho)=%.6e\n",static_cast<long long>(io_.Run().step),io_.Run().time,dt_,global);
}

void MarsSolver::Advance()
{
    auto &run=io_.Run();control_.UpdateSwitches(run);
    while(!control_.if_stop)
    {
        ComputeTimestep();ZeroRHS();AssembleFluidFlux();AddFortranSources();AssembleElectricAndCT();ApplyUpdate();ApplyFloors();
        run.dt=dt_;run.time+=dt_;++run.step;
        std::vector<std::string> names;for(auto n:MARS::SpeciesNames)names.push_back("U_"+std::string(n));SyncGroup("Ucell",names);SyncGroup("Bface",{"B_xi","B_eta","B_zeta"},true);
        UpdatePrimitive();UpdateIonVelocity();UpdateMagneticCell();UpdateCurrentCell();
        control_.UpdateSwitches(run);
        if(control_.if_outres)PrintDiagnostics();
        if(control_.if_outfile){io_.WriteTecplotBinFile(run.step,run.time);io_.WriteRestartBinFile(run.step,run.time);io_.WriteRunDataFile();}
    }
}
