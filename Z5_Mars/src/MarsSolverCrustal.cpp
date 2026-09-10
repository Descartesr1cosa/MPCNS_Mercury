#include "MarsSolver.h"
#include "MarsDipoles.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

void MarsSolver::InitializeCrustalField(bool fresh)
{
    const int aids[3]={Aadd_.xi,Aadd_.eta,Aadd_.zeta};
    const int bids[3]={Badd_.xi,Badd_.eta,Badd_.zeta};
    const int total[3]={B_.xi,B_.eta,B_.zeta};
    for(int ib=0;ib<field_->num_blocks();++ib)
        for(int id:{aids[0],aids[1],aids[2],bids[0],bids[1],bids[2]})
        {
            auto &f=field_->field(id,ib); if(!f.is_allocated()) continue;
            const Int3 lo=f.get_lo(),hi=f.get_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
                f(i,j,k,0)=0;
        }
    if(!param_->HasBoo("is_Crustal_Bfield") || !param_->GetBoo("is_Crustal_Bfield")) return;

    MARS::DipoleModel model;
    const std::string path=param_->HasStr("Crustal_Dipoles_File")
        ?param_->GetStr("Crustal_Dipoles_File"):"CASE/setup/Dipoles";
    // Rank-local filesystem/evaluation failures must terminate all ranks before
    // any peer enters the collective edge-owner exchange.
    try
    {
        std::ifstream input(path);
        if(!input) throw std::runtime_error("cannot open "+path);
        model.Read(input,length_ref_,magnetic_ref_,mu0_,path);
        for(int ib=0;ib<field_->num_blocks();++ib)
        {
            auto &g=grid_->grids(ib);
            for(int axis=0;axis<3;++axis)
            {
                auto &a=field_->field(aids[axis],ib); if(!a.is_allocated()) continue;
                const Int3 lo=a.get_lo(),hi=a.get_hi();
                for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
                {
                    const int ip=i+(axis==0),jp=j+(axis==1),kp=k+(axis==2);
                    a(i,j,k,0)=model.EdgeIntegral({g.x(i,j,k),g.y(i,j,k),g.z(i,j,k)},
                                                 {g.x(ip,jp,kp),g.y(ip,jp,kp),g.z(ip,jp,kp)});
                }
            }
        }
    }
    catch(const std::exception &error)
    {
        std::cerr<<"Mars crustal field rank="<<param_->GetInt("myid")<<": "<<error.what()<<'\n';
        MPI_Abort(MPI_COMM_WORLD,1);
        throw;
    }
    // Includes rotated interfaces, MPI aliases and the singular-edge classes.
    // Only topology synchronization: electric wall/farfield BCs do not apply to A.
    halo_->sync_owner_alias_group("Aadd");
    halo_->sync_group("Aadd");
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        if(!field_->field(bids[0],ib).is_allocated()) continue;
        CTOperators::CurlEdgeToFace(ib,field_->field(aids[0],ib),field_->field(aids[1],ib),
            field_->field(aids[2],ib),field_->field(bids[0],ib),field_->field(bids[1],ib),
            field_->field(bids[2],ib),1.0,true);
    }
    halo_->sync_group("Badd");
    double local[2]={0,0},global[2]={0,0};
    for(int ib=0;ib<field_->num_blocks();++ib)
    {
        auto &cell=field_->field(bcell_,ib); if(!cell.is_allocated()) continue;
        const auto &x=field_->field(bids[0],ib),&y=field_->field(bids[1],ib),&z=field_->field(bids[2],ib);
        const Int3 lo=cell.inner_lo(),hi=cell.inner_hi();
        for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
        {
            const double f[6]={x(i+1,j,k,0),-x(i,j,k,0),y(i,j+1,k,0),-y(i,j,k,0),z(i,j,k+1,0),-z(i,j,k,0)};
            double sum=0,scale=0;for(double v:f){sum+=v;scale+=std::abs(v);}
            local[0]=std::max(local[0],std::abs(sum));local[1]=std::max(local[1],scale);
        }
    }
    MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(global[0]>1.0e-11*global[1])
    {
        std::cerr<<"Mars Badd: nonzero cell flux balance after topology synchronization\n";
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    if(param_->GetInt("myid")==0)
        std::cout<<"[Mars Badd] "<<model.size()<<" dipoles from "<<path
                 <<", max |sum outward flux|="<<global[0]
                 <<", relative="<<(global[1]>0?global[0]/global[1]:0)
                 <<(fresh?", added to initial total field\n":", restart total field retained\n");
    if(fresh)
    {
        // B_* remains TOTAL magnetic flux throughout MarsZ5, including restart,
        // ideal UCT, Hall dissipation, cell reconstruction and divB diagnostics.
        // Badd itself is never advanced and must never be added again during
        // reconstruction. Apply ordinary total-field ghost BCs after addition
        // so fresh initialization and restart have the same ghost contract.
        for(int ib=0;ib<field_->num_blocks();++ib)for(int d=0;d<3;++d)
        {
            auto &b=field_->field(total[d],ib);if(!b.is_allocated())continue;
            const auto &add=field_->field(bids[d],ib);
            const Int3 lo=b.get_lo(),hi=b.get_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
                b(i,j,k,0)+=add(i,j,k,0);
        }
        SyncGroup("Bface",bface_names_,true);
    }
}
