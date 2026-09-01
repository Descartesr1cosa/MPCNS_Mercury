#include "MarsSolver.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>

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
    double fluid_local=std::numeric_limits<double>::max();
    double mhd_local=std::numeric_limits<double>::max();
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &jac=field_->field(jac_,ib), &ax=field_->field(metric_.xi,ib), &ay=field_->field(metric_.eta,ib), &az=field_->field(metric_.zeta,ib);
        if(!jac.is_allocated()) continue;
        auto metric_area=[&](FieldBlock &metric,int i,int j,int k,const char *direction)
        {
            const double x=metric(i,j,k,0),y=metric(i,j,k,1),z=metric(i,j,k,2);
            const double area=std::sqrt(x*x+y*y+z*z);
            if(std::isfinite(area) && area>0.0) return area;
            std::ostringstream message;
            message << "MarsSolver: invalid face metric in ComputeTimestep"
                    << " rank=" << param_->GetInt("myid")
                    << " block=" << ib
                    << " face=" << direction
                    << " index=(" << i << ',' << j << ',' << k << ')'
                    << " area=" << area;
            throw std::runtime_error(message.str());
        };
        const Int3 lo=jac.inner_lo(), hi=jac.inner_hi();
        for(std::size_t s=0;s<MARS::NumSpecies;++s)
        {
            FieldBlock &u=field_->field(U_[s],ib), &pv=field_->field(PV_[s],ib);
            if(!u.is_allocated()) continue;
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                const double volume=std::abs(CheckedCellJacobian(jac,ib,i,j,k,"ComputeTimestep"));
                const double rho=std::max(u(i,j,k,0),1e-10), p=std::max(pv(i,j,k,3),1e-10);
                const double cs=std::sqrt(gamma_*p/rho);
                auto term=[&](FieldBlock &a,int ii,int jj,int kk,const char *direction)
                {
                    const double x=a(ii,jj,kk,0),y=a(ii,jj,kk,1),z=a(ii,jj,kk,2);
                    const double area=metric_area(a,ii,jj,kk,direction);
                    return std::abs(pv(i,j,k,0)*x+pv(i,j,k,1)*y+pv(i,j,k,2)*z)+cs*area;
                };
                const double denom=term(ax,i,j,k,"xi-")+term(ax,i+1,j,k,"xi+")+
                                   term(ay,i,j,k,"eta-")+term(ay,i,j+1,k,"eta+")+
                                   term(az,i,j,k,"zeta-")+term(az,i,j,k+1,"zeta+");
                if(denom>0.0) fluid_local=std::min(fluid_local,cfl_*volume/denom);
            }
        }

        FieldBlock &b=field_->field(bcell_,ib),&up=field_->field(uplus_,ib);
        if(!b.is_allocated()) continue;
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double volume=std::abs(CheckedCellJacobian(jac,ib,i,j,k,"ComputeTimestep/EM"));
            double rho=0.0,p=0.0;
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                rho+=std::max(field_->field(U_[s],ib)(i,j,k,0),0.0);
                p+=std::max(field_->field(PV_[s],ib)(i,j,k,3),0.0);
            }
            rho=std::max(rho,DensityFloor);
            p=std::max(p,PressureFloor);
            const double bx=b(i,j,k,0),by=b(i,j,k,1),bz=b(i,j,k,2);
            const double b2=bx*bx+by*by+bz*bz;
            const double cs2=gamma_*p/rho;
            const double va2=magnetic_force_scale_*b2/rho;
            auto fast_term=[&](FieldBlock &metric,int ii,int jj,int kk,const char *direction)
            {
                const double sx=metric(ii,jj,kk,0),sy=metric(ii,jj,kk,1),sz=metric(ii,jj,kk,2);
                const double area=metric_area(metric,ii,jj,kk,direction);
                const double bn=(bx*sx+by*sy+bz*sz)/area;
                const double van2=magnetic_force_scale_*bn*bn/rho;
                const double sum=cs2+va2;
                const double disc=std::max(0.0,sum*sum-4.0*cs2*van2);
                const double cf=std::sqrt(0.5*(sum+std::sqrt(disc)));
                return std::abs(up(i,j,k,0)*sx+up(i,j,k,1)*sy+up(i,j,k,2)*sz)+cf*area;
            };
            const double mhd_denom=fast_term(ax,i,j,k,"xi-")+fast_term(ax,i+1,j,k,"xi+")+
                                   fast_term(ay,i,j,k,"eta-")+fast_term(ay,i,j+1,k,"eta+")+
                                   fast_term(az,i,j,k,"zeta-")+fast_term(az,i,j,k+1,"zeta+");
            if(mhd_denom>0.0) mhd_local=std::min(mhd_local,cfl_*volume/mhd_denom);

        }
    }
    PARALLEL::mpi_min(&fluid_local,&dt_fluid_,1);
    PARALLEL::mpi_min(&mhd_local,&dt_mhd_,1);
    dt_=std::min(dt_fluid_,dt_mhd_);
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
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double cell_jacobian=CheckedCellJacobian(jac,ib,i,j,k,"AssembleFluidFlux");
            for(int m=0;m<5;++m)
                rhs(i,j,k,m)-=((fx(i+1,j,k,m)-fx(i,j,k,m))+(fy(i,j+1,k,m)-fy(i,j,k,m))+(fz(i,j,k+1,m)-fz(i,j,k,m)))/cell_jacobian;
        }
    }
}

void MarsSolver::AddFortranSources()
{
    const double time_ref=length_ref_/velocity_ref_;
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &neutral=field_->field(neutral_,ib),&up=field_->field(uplus_,ib),&b=field_->field(bcell_,ib),&jv=field_->field(jcell_,ib),&jac=field_->field(jac_,ib);
        FieldBlock &ax=field_->field(metric_.xi,ib),&ay=field_->field(metric_.eta,ib),&az=field_->field(metric_.zeta,ib);
        if(!neutral.is_allocated()) continue;
        const Int3 lo=jac.inner_lo(),hi=jac.inner_hi();
        auto electron_pressure=[&](int i,int jj,int k){double p=0;for(std::size_t s=0;s<MARS::NumSpecies;++s)p+=field_->field(PV_[s],ib)(i,jj,k,3);return p;};
        for(int i=lo.i;i<hi.i;++i) for(int jj=lo.j;jj<hi.j;++jj) for(int k=lo.k;k<hi.k;++k)
        {
            MARS::ChemistryInput input;
            const ElectronDensityInfo density=ElectronDensityAtCell(ib,i,jj,k);
            for(std::size_t n=0;n<MARS::NumNeutral;++n) input.neutral_cm3[n]=std::max(neutral(i,jj,k,n),0.0);
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                FieldBlock &u=field_->field(U_[s],ib);
                input.ion_cm3[s]=std::max(u(i,jj,k,0),0.0)*density_ref_/MARS::ParticleMass[s]*1e-6;
            }

            const double V=std::abs(CheckedCellJacobian(jac,ib,i,jj,k,"AddFortranSources"));
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
                const double rho=std::max(u(i,jj,k,0),0.0),ns=input.ion_cm3[s];
                // Mercury's w_s=n_s/ne_eff: preserve normal-region charge
                // fractions while smoothly removing MHD forcing in vacuum.
                const double fraction=density.species_nd[s]/density.ne_eff;
                double relxb[3];cross3(pv(i,jj,k,0)-up(i,jj,k,0),pv(i,jj,k,1)-up(i,jj,k,1),pv(i,jj,k,2)-up(i,jj,k,2),b(i,jj,k,0),b(i,jj,k,1),b(i,jj,k,2),relxb[0],relxb[1],relxb[2]);
                const double gyro=charge_*magnetic_ref_*length_ref_/(MARS::ParticleMass[s]*velocity_ref_);
                double electromagnetic_force[3];
                double force[3];
                for(int d=0;d<3;++d)
                {
                    electromagnetic_force[d]=rho*gyro*relxb[d]+
                        fraction*magnetic_force_scale_*jxb[d]-
                        (ambipolar_enabled_?fraction*gradp[d]:0.0);
                    force[d]=electromagnetic_force[d]-rho*time_ref*collision*pv(i,jj,k,d);
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
                rhs(i,jj,k,4)+=electromagnetic_force[0]*pv(i,jj,k,0)+electromagnetic_force[1]*pv(i,jj,k,1)+electromagnetic_force[2]*pv(i,jj,k,2);
                rhs(i,jj,k,4)+=a5*ns*b1+time_ref*rho*speed2*b2;
            }
        }
    }
}

void MarsSolver::ApplyPointImplicitChemistry()
{
    const double time_ref=length_ref_/velocity_ref_;
    const double dt_seconds=dt_*time_ref;
    const double energy_source_scale=1.0e6*length_ref_*boltzmann_/
        (density_ref_*velocity_ref_*velocity_ref_*velocity_ref_*(gamma_-1.0));

    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &neutral=field_->field(neutral_,ib),&photo=field_->field(photo_,ib);
        if(!neutral.is_allocated()) continue;
        const Int3 lo=neutral.inner_lo(),hi=neutral.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            MARS::ChemistryInput input;
            std::array<double,MARS::NumSpecies> density_floor_cm3{};
            std::array<double,MARS::NumSpecies> energy_star{};
            std::array<double,MARS::NumSpecies> speed2_star{};
            double electron_density_m3=0.0;
            double electron_pressure_Pa=0.0;

            for(std::size_t n=0;n<MARS::NumNeutral;++n)
                input.neutral_cm3[n]=std::max(neutral(i,j,k,n),0.0);
            for(std::size_t n=0;n<MARS::NumPhoto;++n)
                input.photo_cm3_s[n]=std::max(photo(i,j,k,n),0.0);

            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                FieldBlock &u=field_->field(U_[s],ib);
                const double rho=std::max(u(i,j,k,0),DensityFloor);
                const double mx=u(i,j,k,1),my=u(i,j,k,2),mz=u(i,j,k,3);
                const double speed2=(mx*mx+my*my+mz*mz)/(rho*rho);
                const double kinetic=0.5*rho*speed2;
                const double total_energy=std::max(u(i,j,k,4),kinetic+PressureFloor/(gamma_-1.0));
                const double pressure=std::max((gamma_-1.0)*(total_energy-kinetic),PressureFloor);
                const double number_density_m3=rho*density_ref_/MARS::ParticleMass[s];

                input.ion_cm3[s]=number_density_m3*1.0e-6;
                density_floor_cm3[s]=DensityFloor*density_ref_/MARS::ParticleMass[s]*1.0e-6;
                input.ion_temperature_K[s]=pressure*pressure_ref_/
                    (std::max(number_density_m3,1.0)*boltzmann_);
                energy_star[s]=total_energy;
                speed2_star[s]=speed2;
                electron_density_m3+=number_density_m3;
                electron_pressure_Pa+=pressure*pressure_ref_;
            }
            input.electron_temperature_K=electron_pressure_Pa/
                (std::max(electron_density_m3,1.0)*boltzmann_);

            MARS::PointImplicitChemistryResult chemistry=
                MARS::SolvePointImplicitChemistry(input,dt_seconds,density_floor_cm3,16,1.0e-8);
            for(int subdivisions=2;!chemistry.converged && subdivisions<=16;subdivisions*=2)
            {
                MARS::ChemistryInput substep_input=input;
                bool converged=true;
                int total_iterations=0;
                MARS::PointImplicitChemistryResult substep_result;
                MARS::ChemistryRates averaged_rates;
                for(int substep=0;substep<subdivisions;++substep)
                {
                    substep_input.ion_cm3=substep==0?input.ion_cm3:substep_result.ion_cm3;
                    substep_result=MARS::SolvePointImplicitChemistry(
                        substep_input,dt_seconds/static_cast<double>(subdivisions),
                        density_floor_cm3,16,1.0e-8);
                    total_iterations+=substep_result.iterations;
                    if(!substep_result.converged){converged=false;break;}
                    for(std::size_t s=0;s<MARS::NumSpecies;++s)
                    {
                        averaged_rates.production_cm3_s[s]+=substep_result.rates.production_cm3_s[s]/subdivisions;
                        averaged_rates.loss_cm3_s[s]+=substep_result.rates.loss_cm3_s[s]/subdivisions;
                    }
                    averaged_rates.electron_production_cm3_s+=
                        substep_result.rates.electron_production_cm3_s/subdivisions;
                    averaged_rates.electron_loss_cm3_s+=
                        substep_result.rates.electron_loss_cm3_s/subdivisions;
                }
                if(converged)
                {
                    chemistry=substep_result;
                    chemistry.iterations=total_iterations;
                    chemistry.rates=averaged_rates;
                }
            }
            if(!chemistry.converged)
            {
                std::ostringstream message;
                message << "MarsSolver: point-implicit chemistry failed"
                        << " rank=" << param_->GetInt("myid")
                        << " block=" << ib << " cell=(" << i << ',' << j << ',' << k << ')'
                        << " dt=" << dt_ << " dt_seconds=" << dt_seconds
                        << " iterations=" << chemistry.iterations
                        << " attempted_substeps=16";
                throw std::runtime_error(message.str());
            }

            double electron_density_new_cm3=0.0;
            for(double value:chemistry.ion_cm3) electron_density_new_cm3+=value;
            for(std::size_t s=0;s<MARS::NumSpecies;++s)
            {
                FieldBlock &u=field_->field(U_[s],ib);
                const double number_density_new=chemistry.ion_cm3[s];
                const double loss_coefficient_s=chemistry.rates.loss_cm3_s[s]/
                    std::max(number_density_new,1.0e-300);
                const double momentum_denominator=1.0+dt_seconds*loss_coefficient_s;
                for(int d=1;d<4;++d) u(i,j,k,d)/=momentum_denominator;
                const bool momentum_finite=std::isfinite(u(i,j,k,1)) &&
                    std::isfinite(u(i,j,k,2)) && std::isfinite(u(i,j,k,3));

                const double loss_mass_rhs=chemistry.rates.loss_cm3_s[s]*1.0e6*
                    MARS::ParticleMass[s]/density_ref_*time_ref;
                const double electron_fraction=number_density_new/
                    std::max(electron_density_new_cm3,1.0e-300);
                const double positive_energy_rhs=energy_source_scale*(
                    chemistry.rates.production_cm3_s[s]*185.0+
                    electron_fraction*chemistry.rates.electron_production_cm3_s*185.0);
                const double negative_energy_rhs=energy_source_scale*(
                    chemistry.rates.loss_cm3_s[s]*std::max(input.ion_temperature_K[s],50.0)+
                    electron_fraction*chemistry.rates.electron_loss_cm3_s*
                        std::max(input.electron_temperature_K,50.0))+
                    0.5*loss_mass_rhs*speed2_star[s];
                const double energy_denominator=1.0+dt_*negative_energy_rhs/
                    std::max(energy_star[s],PressureFloor/(gamma_-1.0));

                u(i,j,k,0)=number_density_new*1.0e6*MARS::ParticleMass[s]/density_ref_;
                u(i,j,k,4)=(energy_star[s]+dt_*positive_energy_rhs)/energy_denominator;
                if(!std::isfinite(u(i,j,k,0)) || !momentum_finite ||
                   !std::isfinite(u(i,j,k,4)) ||
                   !(momentum_denominator>0.0) || !(energy_denominator>0.0))
                {
                    std::ostringstream message;
                    message << "MarsSolver: invalid point-implicit chemistry state"
                            << " rank=" << param_->GetInt("myid")
                            << " block=" << ib << " cell=(" << i << ',' << j << ',' << k << ')'
                            << " species=" << MARS::SpeciesNames[s]
                            << " rho=" << u(i,j,k,0) << " E=" << u(i,j,k,4)
                            << " momentum_denominator=" << momentum_denominator
                            << " energy_denominator=" << energy_denominator;
                    throw std::runtime_error(message.str());
                }
            }
        }
    }
}

void MarsSolver::AssembleElectricAndCT()
{
    ZeroEdgeElectricFields();
    AssembleIdealEdgeEMF();
    AssembleAmbipolarEdgeEMF();

    // A rotated regular edge receives one complete UCT candidate per alias.
    // Reconcile those oriented candidates before the owner synchronization.
    ReduceEdgeAliasCandidatesToOwners(E_);
    AssembleSingularEdgeEMFNonHall();

    if(hall_enabled_)
    {
        AssembleHallRusanovEdgeEMF();
        AssembleSingularEdgeEMFHall();
        SyncGroup("Ehall",{"Ehall_xi","Ehall_eta","Ehall_zeta"},true);
        const int eids[3]={E_.xi,E_.eta,E_.zeta};
        const int hids[3]={Ehall_.xi,Ehall_.eta,Ehall_.zeta};
        for(int ib=0;ib<field_->num_blocks();++ib) for(int d=0;d<3;++d)
        {
            FieldBlock &edge=field_->field(eids[d],ib),&hall=field_->field(hids[d],ib);
            if(!edge.is_allocated()) continue;
            const Int3 lo=edge.inner_lo(),hi=edge.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
                edge(i,j,k,0)+=hall(i,j,k,0);
        }
    }

    ApplyStationaryWallEMF();
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

void MarsSolver::UpdateDivB()
{
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &divb=field_->field(divb_,ib),&bx=field_->field(B_.xi,ib),
                   &by=field_->field(B_.eta,ib),&bz=field_->field(B_.zeta,ib),
                   &jac=field_->field(jac_,ib);
        if(!divb.is_allocated()) continue;
        const Int3 lo=divb.inner_lo(),hi=divb.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double cell_jacobian=CheckedCellJacobian(jac,ib,i,j,k,"UpdateDivB");
            divb(i,j,k,0)=((bx(i+1,j,k,0)-bx(i,j,k,0))+
                           (by(i,j+1,k,0)-by(i,j,k,0))+
                           (bz(i,j,k+1,0)-bz(i,j,k,0)))/cell_jacobian;
        }
    }
}

void MarsSolver::PrintDiagnostics()
{
    constexpr std::size_t fluid_values=3*MARS::NumSpecies;
    std::array<double,fluid_values> fluid_min_l,fluid_max_l,fluid_min_g,fluid_max_g;
    fluid_min_l.fill(std::numeric_limits<double>::infinity());
    fluid_max_l.fill(-std::numeric_limits<double>::infinity());

    std::array<double,3> bmin_l,bmax_l,jmin_l,jmax_l,bmin_g,bmax_g,jmin_g,jmax_g;
    bmin_l.fill(std::numeric_limits<double>::infinity());
    jmin_l.fill(std::numeric_limits<double>::infinity());
    bmax_l.fill(-std::numeric_limits<double>::infinity());
    jmax_l.fill(-std::numeric_limits<double>::infinity());
    double pmag_min_l=std::numeric_limits<double>::infinity(),pmag_max_l=-std::numeric_limits<double>::infinity();
    double babs_max_l=0.0,jabs_max_l=0.0,divb_max_l=0.0;

    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        for(std::size_t s=0;s<MARS::NumSpecies;++s)
        {
            FieldBlock &u=field_->field(U_[s],ib),&pv=field_->field(PV_[s],ib);
            if(!u.is_allocated()) continue;
            const Int3 lo=u.inner_lo(),hi=u.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                const double speed=std::sqrt(pv(i,j,k,0)*pv(i,j,k,0)+pv(i,j,k,1)*pv(i,j,k,1)+pv(i,j,k,2)*pv(i,j,k,2));
                const double values[3]={u(i,j,k,0),pv(i,j,k,3),speed};
                for(int q=0;q<3;++q)
                {
                    const std::size_t n=3*s+static_cast<std::size_t>(q);
                    fluid_min_l[n]=std::min(fluid_min_l[n],values[q]);
                    fluid_max_l[n]=std::max(fluid_max_l[n],values[q]);
                }
            }
        }

        FieldBlock &b=field_->field(bcell_,ib),&jv=field_->field(jcell_,ib),&divb=field_->field(divb_,ib);
        if(b.is_allocated() && jv.is_allocated() && divb.is_allocated())
        {
            const Int3 lo=b.inner_lo(),hi=b.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                double b2=0.0,j2=0.0;
                for(int d=0;d<3;++d)
                {
                    const double bv=b(i,j,k,d),jvalue=jv(i,j,k,d);
                    bmin_l[d]=std::min(bmin_l[d],bv);bmax_l[d]=std::max(bmax_l[d],bv);
                    jmin_l[d]=std::min(jmin_l[d],jvalue);jmax_l[d]=std::max(jmax_l[d],jvalue);
                    b2+=bv*bv;j2+=jvalue*jvalue;
                }
                const double pmag=0.5*b2;
                pmag_min_l=std::min(pmag_min_l,pmag);pmag_max_l=std::max(pmag_max_l,pmag);
                babs_max_l=std::max(babs_max_l,std::sqrt(b2));
                jabs_max_l=std::max(jabs_max_l,std::sqrt(j2));
                divb_max_l=std::max(divb_max_l,std::abs(divb(i,j,k,0)));
            }
        }

    }

    PARALLEL::mpi_min(fluid_min_l.data(),fluid_min_g.data(),static_cast<int>(fluid_values));
    PARALLEL::mpi_max(fluid_max_l.data(),fluid_max_g.data(),static_cast<int>(fluid_values));
    PARALLEL::mpi_min(bmin_l.data(),bmin_g.data(),3);PARALLEL::mpi_max(bmax_l.data(),bmax_g.data(),3);
    PARALLEL::mpi_min(jmin_l.data(),jmin_g.data(),3);PARALLEL::mpi_max(jmax_l.data(),jmax_g.data(),3);
    double scalar_min_l[1]={pmag_min_l},scalar_min_g[1];
    double scalar_max_l[4]={pmag_max_l,babs_max_l,jabs_max_l,divb_max_l},scalar_max_g[4];
    PARALLEL::mpi_min(scalar_min_l,scalar_min_g,1);PARALLEL::mpi_max(scalar_max_l,scalar_max_g,4);

    if(param_->GetInt("myid")!=0) return;
    const RunData &run=io_.Run();
    const double cpu_time=run.GetMonitorD("Perf/CPUInterval_s",0.0);
    const double speed=run.GetMonitorD("Perf/s_per_MCells_step",0.0);
    const double mdofs=static_cast<double>(run.GetMonitorI("Perf/Cells",0))*1.0e-6;
    std::printf("[  Diag  ] Nstep   =%9d            t  =%12.5e        dt = %10.3e\n",
                run.step,run.time,run.dt);
    std::printf("           %8s=%10.4f s    %8s=%8.3fM    %8s=%10.3f   s/MDOF/step\n",
                "CPU_Time",cpu_time,"DOFs",mdofs,"speed",speed);
    for(std::size_t s=0;s<MARS::NumSpecies;++s)
    {
        const std::string name=std::string(MARS::SpeciesNames[s])+"+";
        std::printf("           %-4s rho=[%.3e, %.3e]  p=[%.3e, %.3e]  |u|=[%.3e, %.3e]\n",
                    name.c_str(),fluid_min_g[3*s],fluid_max_g[3*s],
                    fluid_min_g[3*s+1],fluid_max_g[3*s+1],fluid_min_g[3*s+2],fluid_max_g[3*s+2]);
    }
    std::printf("           B    x=[%.3e, %.3e]  y=[%.3e, %.3e]  z=[%.3e, %.3e]\n",
                bmin_g[0],bmax_g[0],bmin_g[1],bmax_g[1],bmin_g[2],bmax_g[2]);
    std::printf("           J    x=[%.3e, %.3e]  y=[%.3e, %.3e]  z=[%.3e, %.3e]\n",
                jmin_g[0],jmax_g[0],jmin_g[1],jmax_g[1],jmin_g[2],jmax_g[2]);
    std::printf("           EM   Pmag=[%.3e, %.3e]  |B|max=%.3e  |J|max=%.3e  |divB|max=%.3e\n",
                scalar_min_g[0],scalar_max_g[0],scalar_max_g[1],scalar_max_g[2],scalar_max_g[3]);
    std::printf("\n");
    std::fflush(stdout);
}

void MarsSolver::Advance()
{
    auto &run=io_.Run();control_.UpdateSwitches(run);
    if(control_.if_outres)
    {
        UpdateDivB();
        io_.Runtime().UpdateOnOutres(run);
        PrintDiagnostics();
    }
    while(!control_.if_stop)
    {
        ComputeTimestep();ZeroRHS();AssembleFluidFlux();AddFortranSources();AssembleElectricAndCT();ApplyUpdate();ApplyPointImplicitChemistry();ApplyFloors();
        run.dt=dt_;run.time+=dt_;++run.step;
        std::vector<std::string> names;for(auto n:MARS::SpeciesNames)names.push_back("U_"+std::string(n));SyncGroup("Ucell",names);SyncGroup("Bface",{"B_xi","B_eta","B_zeta"},true);
        UpdatePrimitive();UpdateIonVelocity();UpdateMagneticCell();UpdateCurrentCell();
        control_.UpdateSwitches(run);
        if(control_.if_outres || control_.if_outfile)UpdateDivB();
        if(control_.if_outres){io_.Runtime().UpdateOnOutres(run);PrintDiagnostics();}
        if(control_.if_outfile){io_.WriteTecplotBinFile(run.step,run.time);io_.WriteRestartBinFile(run.step,run.time);io_.WriteRunDataFile();}
    }
}
