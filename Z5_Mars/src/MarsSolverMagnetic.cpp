#include "MarsSolver.h"

#include "0_basic/LayoutTraits.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"

#include <algorithm>
#include <cmath>
#include <limits>
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
    // As in Mercury, Eface components live in the receiver block's
    // computational basis and are never exchanged.  AssembleIdealUCTFromFaces
    // evaluates only the face components actually requested by regular edges;
    // singular and stationary-wall edges use their dedicated assembly below.
    AssembleIdealUCTFromFaces();
}

void MarsSolver::CacheOverriddenNonResistiveEdges()
{
    if(overridden_edge_cache_ready_) return;
    const int rank=param_->GetInt("myid");
    overridden_edge_rank_=rank;
    auto insert=[&](int block,int axis,int i,int j,int k)
    {
        overridden_edges_.insert({TOPO::EntityDim::Edge,rank,block,i,j,k,
                                  static_cast<TOPO::EntityAxis>(axis)});
    };

    if(singular_edges_)
        for(const auto &physical:singular_edges_->entries())
            for(const auto &alias:physical.aliases)
                if(alias.edge.rank==rank)
                    overridden_edges_.insert(alias.edge);

    const StaggerLocation locations[3]={StaggerLocation::EdgeXi,
                                        StaggerLocation::EdgeEt,
                                        StaggerLocation::EdgeZe};
    auto mark_wall=[&](int block,const Box3 &node_box,int direction)
    {
        if(block<0 || block>=field_->num_blocks()) return;
        const int normal=std::abs(direction)-1;
        if(normal<0 || normal>2) return;
        for(int axis=0;axis<3;++axis)
        {
            if(axis==normal) continue;
            FieldBlock &edge=field_->field(axis==0?E_.xi:(axis==1?E_.eta:E_.zeta),block);
            if(!edge.is_allocated()) continue;
            Box3 box=LAYOUT::node_box_to_dof_box(locations[axis],node_box);
            const Int3 lo=edge.inner_lo(),hi=edge.inner_hi();
            box.lo.i=std::max(box.lo.i,lo.i);box.lo.j=std::max(box.lo.j,lo.j);box.lo.k=std::max(box.lo.k,lo.k);
            box.hi.i=std::min(box.hi.i,hi.i);box.hi.j=std::min(box.hi.j,hi.j);box.hi.k=std::min(box.hi.k,hi.k);
            for(int i=box.lo.i;i<box.hi.i;++i)for(int j=box.lo.j;j<box.hi.j;++j)
                for(int k=box.lo.k;k<box.hi.k;++k)insert(block,axis,i,j,k);
        }
    };

    for(const auto &patch:topology_->physical_patches)
        if(patch.bc_name=="Solid_Surface" || patch.bc_name=="Coupled-Solid" ||
           patch.bc_name=="Coupled-Fluid")
            mark_wall(patch.this_block,patch.this_box_node,patch.direction);

    // Hall imposes an additional exact zero on tangential Pole edges.  Cache
    // those aliases once instead of scanning all physical patches per edge and
    // per time step.
    for(const auto &patch:topology_->physical_patches)if(patch.bc_name=="Pole")
    {
        const int normal=std::abs(patch.direction)-1;
        if(patch.this_block<0 || patch.this_block>=field_->num_blocks() ||
           normal<0 || normal>2)continue;
        for(int axis=0;axis<3;++axis)if(axis!=normal)
        {
            FieldBlock &edge=field_->field(axis==0?E_.xi:(axis==1?E_.eta:E_.zeta),
                                           patch.this_block);
            if(!edge.is_allocated())continue;
            Box3 box=LAYOUT::node_box_to_dof_box(locations[axis],patch.this_box_node);
            const Int3 lo=edge.inner_lo(),hi=edge.inner_hi();
            box.lo.i=std::max(box.lo.i,lo.i);box.lo.j=std::max(box.lo.j,lo.j);box.lo.k=std::max(box.lo.k,lo.k);
            box.hi.i=std::min(box.hi.i,hi.i);box.hi.j=std::min(box.hi.j,hi.j);box.hi.k=std::min(box.hi.k,hi.k);
            for(int i=box.lo.i;i<box.hi.i;++i)for(int j=box.lo.j;j<box.hi.j;++j)
                for(int k=box.lo.k;k<box.hi.k;++k)
                    hall_pole_edges_.insert({TOPO::EntityDim::Edge,rank,patch.this_block,
                                             i,j,k,static_cast<TOPO::EntityAxis>(axis)});
        }
    }
    for(const auto &patch:topology_->inner_patches)
        if(patch.is_coupling &&
           (patch.this_block_name=="Solid" || patch.nb_block_name=="Solid"))
        {
            mark_wall(patch.this_block,patch.this_box_node,patch.direction);
            mark_wall(patch.nb_block,patch.nb_box_node,patch.nb_direction);
        }
    for(const auto &patch:topology_->parallel_patches)
        if(patch.is_coupling &&
           (patch.this_block_name=="Solid" || patch.nb_block_name=="Solid"))
            mark_wall(patch.this_block,patch.this_box_node,patch.direction);
    overridden_edge_cache_ready_=true;
}

bool MarsSolver::NonResistiveEdgeIsOverridden(int block,int axis,int i,int j,int k)
{
    CacheOverriddenNonResistiveEdges();
    const TOPO::EntityKey key{TOPO::EntityDim::Edge,overridden_edge_rank_,block,i,j,k,
                              static_cast<TOPO::EntityAxis>(axis)};
    return overridden_edges_.find(key)!=overridden_edges_.end();
}

bool MarsSolver::HallPoleEdge(int block,int axis,int i,int j,int k)
{
    CacheOverriddenNonResistiveEdges();
    const TOPO::EntityKey key{TOPO::EntityDim::Edge,overridden_edge_rank_,block,i,j,k,
                              static_cast<TOPO::EntityAxis>(axis)};
    return hall_pole_edges_.find(key)!=hall_pole_edges_.end();
}

void MarsSolver::AssembleIdealUCTFromFaces()
{
    constexpr double Cuct=0.5;
    const double unset=std::numeric_limits<double>::quiet_NaN();
    CacheOverriddenNonResistiveEdges();
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &up=field_->field(uplus_,ib),&jac=field_->field(jac_,ib);
        if(!up.is_allocated()) continue;
        FieldBlock *b[3]={&field_->field(B_.xi,ib),&field_->field(B_.eta,ib),&field_->field(B_.zeta,ib)};
        FieldBlock *background[3]={&field_->field(Badd_.xi,ib),&field_->field(Badd_.eta,ib),&field_->field(Badd_.zeta,ib)};
        FieldBlock *a[3]={&field_->field(metric_.xi,ib),&field_->field(metric_.eta,ib),&field_->field(metric_.zeta,ib)};
        FieldBlock *face[3]={&field_->field(Eface_.xi,ib),&field_->field(Eface_.eta,ib),&field_->field(Eface_.zeta,ib)};

        // The diagonal Eface component is not a physical face candidate.  Use
        // it as a block-local cell cache for u^axis; this removes dozens of
        // identical metric/Jacobian projections per cell without adding a
        // field or changing any arithmetic expression.
        const Int3 cell_lo=jac.get_lo(),cell_hi=jac.get_hi();
        for(int axis=0;axis<3;++axis)
            for(int i=cell_lo.i;i<cell_hi.i;++i)for(int j=cell_lo.j;j<cell_hi.j;++j)
                for(int k=cell_lo.k;k<cell_hi.k;++k)(*face[axis])(i,j,k,axis)=unset;

        // Each face component is consumed only by edges parallel to that
        // component.  Initialise exactly its one-transverse-layer demand slab;
        // the normal direction is never extended.  NaN is an internal cache
        // marker, not a numerical fallback.
        for(int flux_axis=0;flux_axis<3;++flux_axis)for(int component=0;component<3;++component)
        {
            if(component==flux_axis) continue;
            FieldBlock &f=*face[flux_axis];Int3 lo=f.inner_lo(),hi=f.inner_hi();
            const Int3 alo=f.get_lo(),ahi=f.get_hi();
            const int transverse=3-flux_axis-component;
            if(transverse==0){lo.i=std::max(lo.i-1,alo.i);hi.i=std::min(hi.i+1,ahi.i);}
            else if(transverse==1){lo.j=std::max(lo.j-1,alo.j);hi.j=std::min(hi.j+1,ahi.j);}
            else{lo.k=std::max(lo.k-1,alo.k);hi.k=std::min(hi.k+1,ahi.k);}
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
                f(i,j,k,component)=unset;
        }

        auto face_beta=[&](int axis,int i,int j,int k){return (*b[axis])(i,j,k,0);};
        auto cell_beta=[&](int axis,int i,int j,int k)
        {
            int ip=i,jp=j,kp=k;ShiftCell(axis,1,ip,jp,kp);
            return 0.5*(face_beta(axis,i,j,k)+face_beta(axis,ip,jp,kp));
        };
        // Mercury convention: central flux uses total B, dissipation uses Bind.
        auto face_total=[&](int axis,int i,int j,int k)
        {return face_beta(axis,i,j,k)+(*background[axis])(i,j,k,0);};
        auto cell_total=[&](int axis,int i,int j,int k)
        {
            int ip=i,jp=j,kp=k;ShiftCell(axis,1,ip,jp,kp);
            return 0.5*(face_total(axis,i,j,k)+face_total(axis,ip,jp,kp));
        };
        auto metric=[&](int axis,int i,int j,int k,int component)
        {return (*a[axis])(i,j,k,component);};
        auto ucontra=[&](int axis,int i,int j,int k,const char *operation)
        {
            double &cached=(*face[axis])(i,j,k,axis);
            if(!std::isnan(cached))return cached;
            int ip=i,jp=j,kp=k;ShiftCell(axis,1,ip,jp,kp);
            const double kx=0.5*(metric(axis,i,j,k,0)+metric(axis,ip,jp,kp,0));
            const double ky=0.5*(metric(axis,i,j,k,1)+metric(axis,ip,jp,kp,1));
            const double kz=0.5*(metric(axis,i,j,k,2)+metric(axis,ip,jp,kp,2));
            const double volume=std::abs(CheckedCellJacobian(jac,ib,i,j,k,operation));
            cached=(kx*up(i,j,k,0)+ky*up(i,j,k,1)+kz*up(i,j,k,2))/volume;
            return cached;
        };
        auto face_component=[&](int flux_axis,int component,int i,int j,int k)
        {
            FieldBlock &f=*face[flux_axis];
            if(!InBox(Box3{f.get_lo(),f.get_hi()},i,j,k))
                throw std::runtime_error("MarsSolver: ideal UCT requested an unallocated face candidate");
            double &cached=f(i,j,k,component);
            if(!std::isnan(cached)) return cached;
            int iL=i,jL=j,kL=k,iR=i,jR=j,kR=k;ShiftCell(flux_axis,-1,iL,jL,kL);
            const double beta_normal=face_total(flux_axis,i,j,k);
            auto flux=[&](int beta_axis)
            {
                const double betaL=cell_total(beta_axis,iL,jL,kL),betaR=cell_total(beta_axis,iR,jR,kR);
                const double jump_ind=cell_beta(beta_axis,iR,jR,kR)-cell_beta(beta_axis,iL,jL,kL);
                const double ufL=ucontra(flux_axis,iL,jL,kL,"IdealFaceEMF"),
                             ufR=ucontra(flux_axis,iR,jR,kR,"IdealFaceEMF");
                const double ubL=ucontra(beta_axis,iL,jL,kL,"IdealFaceEMF"),
                             ubR=ucontra(beta_axis,iR,jR,kR,"IdealFaceEMF");
                const double gL=ufL*betaL-ubL*beta_normal,gR=ufR*betaR-ubR*beta_normal;
                return 0.5*(gL+gR)-0.5*std::max(std::abs(ufL),std::abs(ufR))*jump_ind;
            };
            if(flux_axis==0) cached=component==1?flux(2):-flux(1);
            else if(flux_axis==1) cached=component==0?-flux(2):flux(0);
            else cached=component==0?flux(1):-flux(0);
            if(!std::isfinite(cached))
            {
                std::ostringstream message;message<<"MarsSolver: non-finite ideal face EMF rank="
                    <<param_->GetInt("myid")<<" block="<<ib<<" face_axis="<<flux_axis
                    <<" component="<<component<<" face=("<<i<<','<<j<<','<<k<<')';
                throw std::runtime_error(message.str());
            }
            return cached;
        };
        auto speed4=[&](int axis,int i0,int j0,int k0,int i1,int j1,int k1,
                        int i2,int j2,int k2,int i3,int j3,int k3)->double
        {
            return std::max(std::max(std::abs(ucontra(axis,i0,j0,k0,"IdealUCT")),std::abs(ucontra(axis,i1,j1,k1,"IdealUCT"))),
                            std::max(std::abs(ucontra(axis,i2,j2,k2,"IdealUCT")),std::abs(ucontra(axis,i3,j3,k3,"IdealUCT"))));
        };
        {
            FieldBlock &e=field_->field(E_.xi,ib);const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                if(NonResistiveEdgeIsOverridden(ib,0,i,j,k)) continue;
                const double ec=0.25*(face_component(1,0,i,j,k)+face_component(1,0,i,j,k-1)+
                                      face_component(2,0,i,j,k)+face_component(2,0,i,j-1,k));
                const double aeta=speed4(1,i,j,k,i,j-1,k,i,j,k-1,i,j-1,k-1);
                const double azeta=speed4(2,i,j,k,i,j-1,k,i,j,k-1,i,j-1,k-1);
                e(i,j,k,0)=ec+Cuct*(0.5*aeta*(face_beta(2,i,j,k)-face_beta(2,i,j-1,k))-
                                     0.5*azeta*(face_beta(1,i,j,k)-face_beta(1,i,j,k-1)));
            }
        }
        {
            FieldBlock &e=field_->field(E_.eta,ib);const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                if(NonResistiveEdgeIsOverridden(ib,1,i,j,k)) continue;
                const double ec=0.25*(face_component(0,1,i,j,k)+face_component(0,1,i,j,k-1)+
                                      face_component(2,1,i,j,k)+face_component(2,1,i-1,j,k));
                const double azeta=speed4(2,i,j,k,i-1,j,k,i,j,k-1,i-1,j,k-1);
                const double axi=speed4(0,i,j,k,i-1,j,k,i,j,k-1,i-1,j,k-1);
                e(i,j,k,0)=ec+Cuct*(0.5*azeta*(face_beta(0,i,j,k)-face_beta(0,i,j,k-1))-
                                     0.5*axi*(face_beta(2,i,j,k)-face_beta(2,i-1,j,k)));
            }
        }
        {
            FieldBlock &e=field_->field(E_.zeta,ib);const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
            {
                if(NonResistiveEdgeIsOverridden(ib,2,i,j,k)) continue;
                const double ec=0.25*(face_component(0,2,i,j,k)+face_component(0,2,i,j-1,k)+
                                      face_component(1,2,i,j,k)+face_component(1,2,i-1,j,k));
                const double axi=speed4(0,i,j,k,i-1,j,k,i,j-1,k,i-1,j-1,k);
                const double aeta=speed4(1,i,j,k,i-1,j,k,i,j-1,k,i-1,j-1,k);
                e(i,j,k,0)=ec+Cuct*(0.5*axi*(face_beta(1,i,j,k)-face_beta(1,i-1,j,k))-
                                     0.5*aeta*(face_beta(0,i,j,k)-face_beta(0,i,j-1,k)));
            }
        }

        // The diagonal work components were zero before this optimization and
        // are not part of the UCT tensor.  Restore them so diagnostic/output
        // state remains unchanged.
        for(int axis=0;axis<3;++axis)
            for(int i=cell_lo.i;i<cell_hi.i;++i)for(int j=cell_lo.j;j<cell_hi.j;++j)
                for(int k=cell_lo.k;k<cell_hi.k;++k)(*face[axis])(i,j,k,axis)=0.0;
    }
}

void MarsSolver::AssembleAmbipolarEdgeEMF()
{
    if(!ambipolar_enabled_) return;
    const int eids[3]={E_.xi,E_.eta,E_.zeta};
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        FieldBlock &electron=field_->field(hall_ecell_,ib);
        auto pe_cell=[&](int i,int j,int k)
        {
            return electron(i,j,k,2);
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
                ne+=0.25*electron(ii,jj,kk,0);
            }
            return std::sqrt(ne*ne+electron_density_floor_nd_*electron_density_floor_nd_);
        };
        for(int d=0;d<3;++d)
        {
            FieldBlock &e=field_->field(eids[d],ib);if(!e.is_allocated())continue;
            const Int3 lo=e.inner_lo(),hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
            {
                if(NonResistiveEdgeIsOverridden(ib,d,i,j,k)) continue;
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
        Int3 clo=bc.inner_lo(),chi=bc.inner_hi();
        const Int3 alloc_lo=bc.get_lo(),alloc_hi=bc.get_hi();
        clo.i=std::max(clo.i-1,alloc_lo.i);clo.j=std::max(clo.j-1,alloc_lo.j);clo.k=std::max(clo.k-1,alloc_lo.k);
        chi.i=std::min(chi.i+1,alloc_hi.i);chi.j=std::min(chi.j+1,alloc_hi.j);chi.k=std::min(chi.k+1,alloc_hi.k);
        for(int i=clo.i;i<chi.i;++i)for(int j=clo.j;j<chi.j;++j)for(int k=clo.k;k<chi.k;++k)
        {
            const double alpha=hall_scale_/ehc(i,j,k,1);
            const double jx=jc(i,j,k,0),jy=jc(i,j,k,1),jz=jc(i,j,k,2);
            const double bx=bc(i,j,k,0),by=bc(i,j,k,1),bz=bc(i,j,k,2);
            ehc(i,j,k,0)=alpha*(jy*bz-jz*by);
            ehc(i,j,k,1)=alpha*(jz*bx-jx*bz);
            ehc(i,j,k,2)=alpha*(jx*by-jy*bx);
            beta(i,j,k,0)=std::abs(alpha)*std::sqrt(bx*bx+by*by+bz*bz);
        }
        auto avg4=[](double a,double b,double c,double d){return 0.25*(a+b+c+d);};
        auto max4=[](double a,double b,double c,double d){return std::max(std::max(a,b),std::max(c,d));};
        for(int d=0;d<3;++d)
        {
            FieldBlock &edge=field_->field(hids[d],ib),&dr=field_->field(drids[d],ib),
                       &dl=field_->field(dlids[d],ib);
            FieldBlock &b0=field_->field(bids[0],ib),&b1=field_->field(bids[1],ib),&b2=field_->field(bids[2],ib);
            FieldBlock &a0=field_->field(aids[0],ib),&a1=field_->field(aids[1],ib),&a2=field_->field(aids[2],ib);
            const Int3 lo=edge.inner_lo(),hi=edge.inner_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
            {
                if(NonResistiveEdgeIsOverridden(ib,d,i,j,k) || HallPoleEdge(ib,d,i,j,k))
                {edge(i,j,k,0)=0.0;continue;}
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
            emf-=pressure_scale_*inc.source_orientation*dp/
                 field_->field(hall_ecell_,c.block)(c.i,c.j,c.k,1);
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
