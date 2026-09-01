#include "MarsSolver.h"

#include "0_basic/LayoutTraits.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace
{
void ShiftCell(int axis, int side, int &i, int &j, int &k)
{
    if(axis==0) i+=side;
    else if(axis==1) j+=side;
    else k+=side;
}

bool InBox(const Box3 &box, int i, int j, int k)
{
    return i>=box.lo.i && i<box.hi.i &&
           j>=box.lo.j && j<box.hi.j &&
           k>=box.lo.k && k<box.hi.k;
}
}

void MarsSolver::ZeroEdgeElectricFields()
{
    const int ids[6]={E_.xi,E_.eta,E_.zeta,Ehall_.xi,Ehall_.eta,Ehall_.zeta};
    for(int ib=0;ib<field_->num_blocks();++ib) for(int id:ids)
    {
        FieldBlock &edge=field_->field(id,ib);
        if(!edge.is_allocated()) continue;
        const Int3 lo=edge.get_lo(),hi=edge.get_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            edge(i,j,k,0)=0.0;
    }
}

void MarsSolver::AssembleIdealEdgeEMF()
{
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &up=field_->field(uplus_,ib);
        if(!up.is_allocated()) continue;
        FieldBlock &bx=field_->field(B_.xi,ib),&by=field_->field(B_.eta,ib),&bz=field_->field(B_.zeta,ib);
        FieldBlock &jac=field_->field(jac_,ib),&ax=field_->field(metric_.xi,ib),
                   &ay=field_->field(metric_.eta,ib),&az=field_->field(metric_.zeta,ib);
        AssembleOneDirectionIdealFaceEMF(ib,0,field_->field(Eface_.xi,ib),bx,by,bz,jac,ax,ay,az,up);
        AssembleOneDirectionIdealFaceEMF(ib,1,field_->field(Eface_.eta,ib),bx,by,bz,jac,ax,ay,az,up);
        AssembleOneDirectionIdealFaceEMF(ib,2,field_->field(Eface_.zeta,ib),bx,by,bz,jac,ax,ay,az,up);
    }

    // Eface components live in the receiver block's computational basis.  As
    // in Mercury, fill the consumed ghost layer locally and never exchange it.
    AssembleIdealUCTFromFaces();
}

void MarsSolver::AssembleOneDirectionIdealFaceEMF(
    int block,int flux_axis,FieldBlock &eface,FieldBlock &bx,FieldBlock &by,FieldBlock &bz,
    FieldBlock &jac,FieldBlock &ax,FieldBlock &ay,FieldBlock &az,FieldBlock &up)
{
    if(flux_axis<0 || flux_axis>2)
        throw std::runtime_error("MarsSolver: invalid ideal EMF direction");

    auto face_beta=[&](int axis,int i,int j,int k)->double
    {
        return axis==0?bx(i,j,k,0):(axis==1?by(i,j,k,0):bz(i,j,k,0));
    };
    auto cell_beta=[&](int axis,int i,int j,int k)->double
    {
        int ip=i,jp=j,kp=k;ShiftCell(axis,1,ip,jp,kp);
        return 0.5*(face_beta(axis,i,j,k)+face_beta(axis,ip,jp,kp));
    };
    auto metric=[&](int axis,int i,int j,int k,int component)->double
    {
        return axis==0?ax(i,j,k,component):(axis==1?ay(i,j,k,component):az(i,j,k,component));
    };
    auto cell_u_contra=[&](int axis,int i,int j,int k)->double
    {
        int ip=i,jp=j,kp=k;ShiftCell(axis,1,ip,jp,kp);
        const double kx=0.5*(metric(axis,i,j,k,0)+metric(axis,ip,jp,kp,0));
        const double ky=0.5*(metric(axis,i,j,k,1)+metric(axis,ip,jp,kp,1));
        const double kz=0.5*(metric(axis,i,j,k,2)+metric(axis,ip,jp,kp,2));
        const double volume=std::abs(CheckedCellJacobian(jac,block,i,j,k,
                                                         "AssembleOneDirectionIdealFaceEMF"));
        return (kx*up(i,j,k,0)+ky*up(i,j,k,1)+kz*up(i,j,k,2))/volume;
    };
    auto flux=[&](int beta_axis,int iL,int jL,int kL,int iR,int jR,int kR,
                  double beta_normal)->double
    {
        const double betaL=cell_beta(beta_axis,iL,jL,kL),betaR=cell_beta(beta_axis,iR,jR,kR);
        const double ufL=cell_u_contra(flux_axis,iL,jL,kL),ufR=cell_u_contra(flux_axis,iR,jR,kR);
        const double ubL=cell_u_contra(beta_axis,iL,jL,kL),ubR=cell_u_contra(beta_axis,iR,jR,kR);
        const double gL=ufL*betaL-ubL*beta_normal,gR=ufR*betaR-ubR*beta_normal;
        return 0.5*(gL+gR)-0.5*std::max(std::abs(ufL),std::abs(ufR))*(betaR-betaL);
    };

    Int3 lo=eface.inner_lo(),hi=eface.inner_hi();
    const Int3 alo=eface.get_lo(),ahi=eface.get_hi();
    lo.i=std::max(lo.i-1,alo.i+1);lo.j=std::max(lo.j-1,alo.j+1);lo.k=std::max(lo.k-1,alo.k+1);
    hi.i=std::min(hi.i+1,ahi.i-1);hi.j=std::min(hi.j+1,ahi.j-1);hi.k=std::min(hi.k+1,ahi.k-1);
    for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
    {
        int iL=i,jL=j,kL=k,iR=i,jR=j,kR=k;ShiftCell(flux_axis,-1,iL,jL,kL);
        const double beta_normal=face_beta(flux_axis,i,j,k);
        double e[3]={0.0,0.0,0.0};
        if(flux_axis==0)
        {
            e[1]=flux(2,iL,jL,kL,iR,jR,kR,beta_normal);
            e[2]=-flux(1,iL,jL,kL,iR,jR,kR,beta_normal);
        }
        else if(flux_axis==1)
        {
            e[0]=-flux(2,iL,jL,kL,iR,jR,kR,beta_normal);
            e[2]=flux(0,iL,jL,kL,iR,jR,kR,beta_normal);
        }
        else
        {
            e[0]=flux(1,iL,jL,kL,iR,jR,kR,beta_normal);
            e[1]=-flux(0,iL,jL,kL,iR,jR,kR,beta_normal);
        }
        for(int d=0;d<3;++d) eface(i,j,k,d)=e[d];
    }
}

void MarsSolver::AssembleIdealUCTFromFaces()
{
    constexpr double Cuct=0.5;
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &up=field_->field(uplus_,ib),&jac=field_->field(jac_,ib);
        if(!up.is_allocated()) continue;
        FieldBlock &bx=field_->field(B_.xi,ib),&by=field_->field(B_.eta,ib),&bz=field_->field(B_.zeta,ib);
        FieldBlock &ax=field_->field(metric_.xi,ib),&ay=field_->field(metric_.eta,ib),&az=field_->field(metric_.zeta,ib);
        FieldBlock &fx=field_->field(Eface_.xi,ib),&fy=field_->field(Eface_.eta,ib),&fz=field_->field(Eface_.zeta,ib);
        auto metric=[&](int axis,int i,int j,int k,int c)->double
        {return axis==0?ax(i,j,k,c):(axis==1?ay(i,j,k,c):az(i,j,k,c));};
        auto ucontra=[&](int axis,int i,int j,int k)->double
        {
            int ip=i,jp=j,kp=k;ShiftCell(axis,1,ip,jp,kp);
            const double kx=0.5*(metric(axis,i,j,k,0)+metric(axis,ip,jp,kp,0));
            const double ky=0.5*(metric(axis,i,j,k,1)+metric(axis,ip,jp,kp,1));
            const double kz=0.5*(metric(axis,i,j,k,2)+metric(axis,ip,jp,kp,2));
            const double volume=std::abs(CheckedCellJacobian(jac,ib,i,j,k,"AssembleIdealUCTFromFaces"));
            return (kx*up(i,j,k,0)+ky*up(i,j,k,1)+kz*up(i,j,k,2))/volume;
        };
        auto speed4=[&](int axis,int i0,int j0,int k0,int i1,int j1,int k1,
                        int i2,int j2,int k2,int i3,int j3,int k3)->double
        {
            return std::max(std::max(std::abs(ucontra(axis,i0,j0,k0)),std::abs(ucontra(axis,i1,j1,k1))),
                            std::max(std::abs(ucontra(axis,i2,j2,k2)),std::abs(ucontra(axis,i3,j3,k3))));
        };
        {
            FieldBlock &e=field_->field(E_.xi,ib);const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                const double ec=0.25*(fy(i,j,k,0)+fy(i,j,k-1,0)+fz(i,j,k,0)+fz(i,j-1,k,0));
                const double aeta=speed4(1,i,j,k,i,j-1,k,i,j,k-1,i,j-1,k-1);
                const double azeta=speed4(2,i,j,k,i,j-1,k,i,j,k-1,i,j-1,k-1);
                e(i,j,k,0)=ec+Cuct*(0.5*aeta*(bz(i,j,k,0)-bz(i,j-1,k,0))-
                                     0.5*azeta*(by(i,j,k,0)-by(i,j,k-1,0)));
            }
        }
        {
            FieldBlock &e=field_->field(E_.eta,ib);const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                const double ec=0.25*(fx(i,j,k,1)+fx(i,j,k-1,1)+fz(i,j,k,1)+fz(i-1,j,k,1));
                const double azeta=speed4(2,i,j,k,i-1,j,k,i,j,k-1,i-1,j,k-1);
                const double axi=speed4(0,i,j,k,i-1,j,k,i,j,k-1,i-1,j,k-1);
                e(i,j,k,0)=ec+Cuct*(0.5*azeta*(bx(i,j,k,0)-bx(i,j,k-1,0))-
                                     0.5*axi*(bz(i,j,k,0)-bz(i-1,j,k,0)));
            }
        }
        {
            FieldBlock &e=field_->field(E_.zeta,ib);const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                const double ec=0.25*(fx(i,j,k,2)+fx(i,j-1,k,2)+fy(i,j,k,2)+fy(i-1,j,k,2));
                const double axi=speed4(0,i,j,k,i-1,j,k,i,j-1,k,i-1,j-1,k);
                const double aeta=speed4(1,i,j,k,i-1,j,k,i,j-1,k,i-1,j-1,k);
                e(i,j,k,0)=ec+Cuct*(0.5*axi*(by(i,j,k,0)-by(i-1,j,k,0))-
                                     0.5*aeta*(bx(i,j,k,0)-bx(i,j-1,k,0)));
            }
        }
    }
}

void MarsSolver::AssembleAmbipolarEdgeEMF()
{
    if(!ambipolar_enabled_) return;
    const int eids[3]={E_.xi,E_.eta,E_.zeta};
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        auto pe_cell=[&](int i,int j,int k)
        {
            double p=0.0;for(std::size_t s=0;s<MARS::NumSpecies;++s)p+=field_->field(PV_[s],ib)(i,j,k,3);
            return p;
        };
        auto pe_node=[&](int i,int j,int k)
        {
            double p=0.0;for(int di=-1;di<=0;++di)for(int dj=-1;dj<=0;++dj)for(int dk=-1;dk<=0;++dk)
                p+=pe_cell(i+di,j+dj,k+dk);return 0.125*p;
        };
        auto ne_edge=[&](int d,int i,int j,int k)
        {
            double ne=0.0;for(int a=0;a<2;++a)for(int c=0;c<2;++c)
            {
                int ii=i,jj=j,kk=k;
                if(d==0){jj=j-1+a;kk=k-1+c;}else if(d==1){ii=i-1+a;kk=k-1+c;}else{ii=i-1+a;jj=j-1+c;}
                ne+=0.25*ElectronDensityAtCell(ib,ii,jj,kk).ne_true;
            }
            return std::sqrt(ne*ne+electron_density_floor_nd_*electron_density_floor_nd_);
        };
        for(int d=0;d<3;++d)
        {
            FieldBlock &e=field_->field(eids[d],ib);if(!e.is_allocated())continue;
            const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
            {
                const double p0=pe_node(i,j,k);
                const double p1=d==0?pe_node(i+1,j,k):(d==1?pe_node(i,j+1,k):pe_node(i,j,k+1));
                const double ne=ne_edge(d,i,j,k);
                if(!std::isfinite(p0)||!std::isfinite(p1)||!std::isfinite(ne)||!(ne>0.0))
                    throw std::runtime_error("MarsSolver: invalid state in ambipolar edge EMF");
                e(i,j,k,0)-=pressure_scale_*(p1-p0)/ne;
            }
        }
    }
}

void MarsSolver::AssembleHallRusanovEdgeEMF()
{
    constexpr double Cwh=0.1;
    const int bids[3]={B_.xi,B_.eta,B_.zeta},hids[3]={Ehall_.xi,Ehall_.eta,Ehall_.zeta};
    const int aids[3]={face_area_.xi,face_area_.eta,face_area_.zeta};
    const int dlids[3]={edge_dl_.xi,edge_dl_.eta,edge_dl_.zeta};
    const int drids[3]={edge_dr_.xi,edge_dr_.eta,edge_dr_.zeta};
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &bc=field_->field(bcell_,ib),&jc=field_->field(jcell_,ib);
        if(!bc.is_allocated())continue;
        FieldBlock &ehc=field_->field(hall_ecell_,ib),&beta=field_->field(hall_beta_,ib);
        const Int3 clo=bc.get_lo(),chi=bc.get_hi();
        for(int i=clo.i;i<chi.i;++i)for(int j=clo.j;j<chi.j;++j)for(int k=clo.k;k<chi.k;++k)
        {
            const ElectronDensityInfo density=ElectronDensityAtCell(ib,i,j,k);
            const double alpha=hall_scale_/density.ne_eff;
            const double jx=jc(i,j,k,0),jy=jc(i,j,k,1),jz=jc(i,j,k,2);
            const double bx=bc(i,j,k,0),by=bc(i,j,k,1),bz=bc(i,j,k,2);
            ehc(i,j,k,0)=alpha*(jy*bz-jz*by);
            ehc(i,j,k,1)=alpha*(jz*bx-jx*bz);
            ehc(i,j,k,2)=alpha*(jx*by-jy*bx);
            beta(i,j,k,0)=std::abs(alpha)*std::sqrt(bx*bx+by*by+bz*bz);
        }
        auto avg4=[](double a,double b,double c,double d){return 0.25*(a+b+c+d);};
        auto max4=[](double a,double b,double c,double d){return std::max(std::max(a,b),std::max(c,d));};
        auto pole_edge=[&](int d,int i,int j,int k)
        {
            const StaggerLocation loc=d==0?StaggerLocation::EdgeXi:(d==1?StaggerLocation::EdgeEt:StaggerLocation::EdgeZe);
            for(const auto &p:topology_->physical_patches)if(p.this_block==ib&&p.bc_name=="Pole")
            {
                const int normal=std::abs(p.direction)-1;
                if(d==normal)continue;
                if(InBox(LAYOUT::node_box_to_dof_box(loc,p.this_box_node),i,j,k))return true;
            }
            return false;
        };
        for(int d=0;d<3;++d)
        {
            FieldBlock &edge=field_->field(hids[d],ib),&dr=field_->field(drids[d],ib),
                       &dl=field_->field(dlids[d],ib);
            FieldBlock &b0=field_->field(bids[0],ib),&b1=field_->field(bids[1],ib),&b2=field_->field(bids[2],ib);
            FieldBlock &a0=field_->field(aids[0],ib),&a1=field_->field(aids[1],ib),&a2=field_->field(aids[2],ib);
            const Int3 lo=edge.inner_lo(),hi=edge.inner_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
            {
                if(pole_edge(d,i,j,k)){edge(i,j,k,0)=0.0;continue;}
                int ci[4],cj[4],ck[4];
                if(d==0){int ti[4]={i,i,i,i},tj[4]={j-1,j-1,j,j},tk[4]={k-1,k,k-1,k};std::copy(ti,ti+4,ci);std::copy(tj,tj+4,cj);std::copy(tk,tk+4,ck);}
                else if(d==1){int ti[4]={i-1,i-1,i,i},tj[4]={j,j,j,j},tk[4]={k-1,k,k-1,k};std::copy(ti,ti+4,ci);std::copy(tj,tj+4,cj);std::copy(tk,tk+4,ck);}
                else{int ti[4]={i-1,i-1,i,i},tj[4]={j-1,j,j-1,j},tk[4]={k,k,k,k};std::copy(ti,ti+4,ci);std::copy(tj,tj+4,cj);std::copy(tk,tk+4,ck);}
                double ec[3];for(int q=0;q<3;++q)ec[q]=avg4(ehc(ci[0],cj[0],ck[0],q),ehc(ci[1],cj[1],ck[1],q),ehc(ci[2],cj[2],ck[2],q),ehc(ci[3],cj[3],ck[3],q));
                const double central=ec[0]*dr(i,j,k,0)+ec[1]*dr(i,j,k,1)+ec[2]*dr(i,j,k,2);
                const double be=max4(beta(ci[0],cj[0],ck[0],0),beta(ci[1],cj[1],ck[1],0),beta(ci[2],cj[2],ck[2],0),beta(ci[3],cj[3],ck[3],0));
                const double L=CheckedPositiveGeometry(dl,ib,i,j,k,"AssembleHallRusanovEdgeEMF/dl");
                if(d==0)
                {
                    const double hz=0.5*(CheckedPositiveGeometry(a1,ib,i,j,k,"Hall/Area_eta")+CheckedPositiveGeometry(a1,ib,i,j,k-1,"Hall/Area_eta"))/L;
                    const double he=0.5*(CheckedPositiveGeometry(a2,ib,i,j,k,"Hall/Area_zeta")+CheckedPositiveGeometry(a2,ib,i,j-1,k,"Hall/Area_zeta"))/L;
                    const double mu_eta=Cwh*be/(he*he),mu_zeta=Cwh*be/(hz*hz);
                    edge(i,j,k,0)=central+0.5*mu_eta*(b2(i,j,k,0)-b2(i,j-1,k,0))-0.5*mu_zeta*(b1(i,j,k,0)-b1(i,j,k-1,0));
                }
                else if(d==1)
                {
                    const double hz=0.5*(CheckedPositiveGeometry(a0,ib,i,j,k,"Hall/Area_xi")+CheckedPositiveGeometry(a0,ib,i,j,k-1,"Hall/Area_xi"))/L;
                    const double hx=0.5*(CheckedPositiveGeometry(a2,ib,i,j,k,"Hall/Area_zeta")+CheckedPositiveGeometry(a2,ib,i-1,j,k,"Hall/Area_zeta"))/L;
                    const double mu_zeta=Cwh*be/(hz*hz),mu_xi=Cwh*be/(hx*hx);
                    edge(i,j,k,0)=central+0.5*mu_zeta*(b0(i,j,k,0)-b0(i,j,k-1,0))-0.5*mu_xi*(b2(i,j,k,0)-b2(i-1,j,k,0));
                }
                else
                {
                    const double hx=0.5*(CheckedPositiveGeometry(a1,ib,i,j,k,"Hall/Area_eta")+CheckedPositiveGeometry(a1,ib,i-1,j,k,"Hall/Area_eta"))/L;
                    const double he=0.5*(CheckedPositiveGeometry(a0,ib,i,j,k,"Hall/Area_xi")+CheckedPositiveGeometry(a0,ib,i,j-1,k,"Hall/Area_xi"))/L;
                    const double mu_xi=Cwh*be/(hx*hx),mu_eta=Cwh*be/(he*he);
                    edge(i,j,k,0)=central+0.5*mu_xi*(b1(i,j,k,0)-b1(i-1,j,k,0))-0.5*mu_eta*(b0(i,j,k,0)-b0(i,j-1,k,0));
                }
            }
        }
        ZeroPoleHallEdges(ib);
    }
}

void MarsSolver::ReduceEdgeAliasCandidatesToOwners(const Triplet &ids)
{
    if(!topology_)return;
    const int rank=param_->GetInt("myid");
    if(!edge_alias_reduce_cache_ready_)
    {
        for(const auto &cls:topology_->edges.classes)if(cls.global_id>=0&&cls.members.size()>1)
            edge_alias_reduce_classes_.push_back(&cls);
        std::sort(edge_alias_reduce_classes_.begin(),edge_alias_reduce_classes_.end(),
                  [](const auto *a,const auto *b){return a->global_id<b->global_id;});
        edge_alias_local_sum_.resize(edge_alias_reduce_classes_.size());
        edge_alias_global_sum_.resize(edge_alias_reduce_classes_.size());
        for(std::size_t n=0;n<edge_alias_reduce_classes_.size();++n)
        {
            const auto &cls=*edge_alias_reduce_classes_[n];
            for(const auto &m:cls.members)if(m.entity.rank==rank)
                edge_alias_local_terms_.push_back({n,m.entity,m.orient_sign});
            if(cls.owner.entity.rank==rank)
                edge_alias_owner_writes_.push_back({n,cls.owner.entity,cls.owner.orient_sign,
                                                     static_cast<double>(cls.members.size())});
        }
        edge_alias_reduce_cache_ready_=true;
    }
    if(edge_alias_reduce_classes_.empty())return;
    auto field_id=[&](TOPO::EntityAxis a){return a==TOPO::EntityAxis::Xi?ids.xi:(a==TOPO::EntityAxis::Eta?ids.eta:ids.zeta);};
    std::fill(edge_alias_local_sum_.begin(),edge_alias_local_sum_.end(),0.0);
    for(const auto &term:edge_alias_local_terms_)
    {
        const auto &e=term.edge;
        edge_alias_local_sum_[term.class_index]+=term.orientation*field_->field(field_id(e.axis),e.block)(e.i,e.j,e.k,0);
    }
    MPI_Allreduce(edge_alias_local_sum_.data(),edge_alias_global_sum_.data(),
                  static_cast<int>(edge_alias_local_sum_.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(const auto &write:edge_alias_owner_writes_)
    {
        const auto &e=write.edge;
        field_->field(field_id(e.axis),e.block)(e.i,e.j,e.k,0)=write.orientation*
            edge_alias_global_sum_[write.class_index]/write.member_count;
    }
}

void MarsSolver::AssembleSingularEdgeEMFNonHall()
{
    if(!singular_edges_||singular_edges_->empty())return;
    auto contribution=[&](const METRIC::SingularPhysicalEdge &edge,
                          const METRIC::WeightedIncidentEntity &inc)
    {
        const auto &c=inc.entity;FieldBlock &up=field_->field(uplus_,c.block),&b=field_->field(bcell_,c.block);
        const double ex=-(up(c.i,c.j,c.k,1)*b(c.i,c.j,c.k,2)-up(c.i,c.j,c.k,2)*b(c.i,c.j,c.k,1));
        const double ey=-(up(c.i,c.j,c.k,2)*b(c.i,c.j,c.k,0)-up(c.i,c.j,c.k,0)*b(c.i,c.j,c.k,2));
        const double ez=-(up(c.i,c.j,c.k,0)*b(c.i,c.j,c.k,1)-up(c.i,c.j,c.k,1)*b(c.i,c.j,c.k,0));
        double emf=ex*edge.canonical_edge_vector[0]+ey*edge.canonical_edge_vector[1]+ez*edge.canonical_edge_vector[2];
        if(ambipolar_enabled_)
        {
            int im=c.i,jm=c.j,km=c.k,ip=c.i,jp=c.j,kp=c.k;
            const int axis=static_cast<int>(inc.source_alias.axis);
            ShiftCell(axis,-1,im,jm,km);ShiftCell(axis,1,ip,jp,kp);
            auto pe=[&](int i,int j,int k){double p=0.0;for(std::size_t s=0;s<MARS::NumSpecies;++s)p+=field_->field(PV_[s],c.block)(i,j,k,3);return p;};
            const double dp=0.5*(pe(ip,jp,kp)-pe(im,jm,km));
            emf-=pressure_scale_*inc.source_orientation*dp/ElectronDensityAtCell(c.block,c.i,c.j,c.k).ne_eff;
        }
        return emf;
    };
    singular_edges_->assemble_cell_field_to_local_owners(*field_,"E_xi",contribution);
    singular_edges_->assemble_cell_field_to_local_owners(*field_,"E_eta",contribution);
    singular_edges_->assemble_cell_field_to_local_owners(*field_,"E_zeta",contribution);
}

void MarsSolver::AssembleSingularEdgeEMFHall()
{
    if(!singular_edges_||singular_edges_->empty())return;
    auto contribution=[&](const METRIC::SingularPhysicalEdge &edge,
                          const METRIC::WeightedIncidentEntity &inc)
    {
        const auto &c=inc.entity;FieldBlock &b=field_->field(bcell_,c.block),&j=field_->field(jcell_,c.block);
        const double alpha=hall_scale_/ElectronDensityAtCell(c.block,c.i,c.j,c.k).ne_eff;
        const double ex=alpha*(j(c.i,c.j,c.k,1)*b(c.i,c.j,c.k,2)-j(c.i,c.j,c.k,2)*b(c.i,c.j,c.k,1));
        const double ey=alpha*(j(c.i,c.j,c.k,2)*b(c.i,c.j,c.k,0)-j(c.i,c.j,c.k,0)*b(c.i,c.j,c.k,2));
        const double ez=alpha*(j(c.i,c.j,c.k,0)*b(c.i,c.j,c.k,1)-j(c.i,c.j,c.k,1)*b(c.i,c.j,c.k,0));
        return ex*edge.canonical_edge_vector[0]+ey*edge.canonical_edge_vector[1]+ez*edge.canonical_edge_vector[2];
    };
    singular_edges_->assemble_cell_field_to_local_owners(*field_,"Ehall_xi",contribution);
    singular_edges_->assemble_cell_field_to_local_owners(*field_,"Ehall_eta",contribution);
    singular_edges_->assemble_cell_field_to_local_owners(*field_,"Ehall_zeta",contribution);
}

void MarsSolver::ApplyStationaryWallEMF()
{
    const int ids[3]={E_.xi,E_.eta,E_.zeta};
    const StaggerLocation locs[3]={StaggerLocation::EdgeXi,StaggerLocation::EdgeEt,StaggerLocation::EdgeZe};
    auto zero=[&](int block,const Box3 &nodes,int direction)
    {
        const int normal=std::abs(direction)-1;
        for(int d=0;d<3;++d)if(d!=normal)
        {
            FieldBlock &e=field_->field(ids[d],block);if(!e.is_allocated())continue;
            Box3 box=LAYOUT::node_box_to_dof_box(locs[d],nodes);const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            box.lo.i=std::max(box.lo.i,lo.i);box.lo.j=std::max(box.lo.j,lo.j);box.lo.k=std::max(box.lo.k,lo.k);
            box.hi.i=std::min(box.hi.i,hi.i);box.hi.j=std::min(box.hi.j,hi.j);box.hi.k=std::min(box.hi.k,hi.k);
            for(int i=box.lo.i;i<box.hi.i;++i)for(int j=box.lo.j;j<box.hi.j;++j)for(int k=box.lo.k;k<box.hi.k;++k)e(i,j,k,0)=0.0;
        }
    };
    for(const auto &p:topology_->physical_patches)
        if(p.bc_name=="Solid_Surface"||p.bc_name=="Coupled-Solid"||p.bc_name=="Coupled-Fluid")zero(p.this_block,p.this_box_node,p.direction);
    for(const auto &p:topology_->inner_patches)if(p.is_coupling&&(p.this_block_name=="Solid"||p.nb_block_name=="Solid"))
    {zero(p.this_block,p.this_box_node,p.direction);zero(p.nb_block,p.nb_box_node,p.nb_direction);}
    for(const auto &p:topology_->parallel_patches)if(p.is_coupling&&(p.this_block_name=="Solid"||p.nb_block_name=="Solid"))
        zero(p.this_block,p.this_box_node,p.direction);
}

void MarsSolver::ZeroPoleHallEdges(int block)
{
    const int ids[3]={Ehall_.xi,Ehall_.eta,Ehall_.zeta};
    const StaggerLocation locs[3]={StaggerLocation::EdgeXi,StaggerLocation::EdgeEt,StaggerLocation::EdgeZe};
    for(const auto &p:topology_->physical_patches)if(p.this_block==block&&p.bc_name=="Pole")
    {
        const int normal=std::abs(p.direction)-1;
        for(int d=0;d<3;++d)if(d!=normal)
        {
            FieldBlock &e=field_->field(ids[d],block);Box3 box=LAYOUT::node_box_to_dof_box(locs[d],p.this_box_node);
            const Int3 lo=e.inner_lo(),hi=e.inner_hi();box.lo.i=std::max(box.lo.i,lo.i);box.lo.j=std::max(box.lo.j,lo.j);box.lo.k=std::max(box.lo.k,lo.k);
            box.hi.i=std::min(box.hi.i,hi.i);box.hi.j=std::min(box.hi.j,hi.j);box.hi.k=std::min(box.hi.k,hi.k);
            for(int i=box.lo.i;i<box.hi.i;++i)for(int j=box.lo.j;j<box.hi.j;++j)for(int k=box.lo.k;k<box.hi.k;++k)e(i,j,k,0)=0.0;
        }
    }
}
