#include "MarsSolver.h"
#include "MarsDipoles.h"
#include "MarsDipoleCache.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

void MarsSolver::InitializeMagneticBackground(bool fresh)
{
    const int background[3]={Badd_.xi,Badd_.eta,Badd_.zeta};
    const int metrics[3]={metric_.xi,metric_.eta,metric_.zeta};
    // Saved Badd is for postprocessing only. Reset/rebuild it from current
    // setup/cache without modifying the checkpoint's induced field.
    InitializeCrustalField();
    for(int ib=0;ib<field_->num_blocks();++ib)for(int d=0;d<3;++d)
    {
        auto &add=field_->field(background[d],ib);if(!add.is_allocated())continue;
        const auto &area=field_->field(metrics[d],ib);
        const Int3 lo=add.get_lo(),hi=add.get_hi();
        for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
            for(int c=0;c<3;++c)add(i,j,k,0)+=solar_wind_imf_[c]*area(i,j,k,c);
    }
    halo_->sync_group("Badd");
    CheckBackgroundFlux("IMF + dipoles");
    if(param_->GetInt("myid")==0)
        std::cout<<"[Mars Badd] IMF + current dipoles ready; "
                 <<(fresh?"fresh Bind=0":"restart Bind retained")<<'\n';
}

void MarsSolver::InitializeCrustalField()
{
    const int aids[3]={Aadd_.xi,Aadd_.eta,Aadd_.zeta};
    const int bids[3]={Badd_.xi,Badd_.eta,Badd_.zeta};
    for(int ib=0;ib<field_->num_blocks();++ib)
        for(int id:{aids[0],aids[1],aids[2],bids[0],bids[1],bids[2]})
        {
            auto &f=field_->field(id,ib); if(!f.is_allocated()) continue;
            const Int3 lo=f.get_lo(),hi=f.get_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
                f(i,j,k,0)=0;
        }
    if(!param_->HasBoo("is_Crustal_Bfield") || !param_->GetBoo("is_Crustal_Bfield"))
    {
        if(param_->GetInt("myid")==0)
            std::cout<<"[Mars Dipoles] disabled; IMF is still included in Badd\n";
        return;
    }

    const double started=MPI_Wtime();
    MARS::DipoleModel model;
    MARS::DipoleFingerprint signature;
    bool cached=false;
    std::size_t count=0;
    for(int ib=0;ib<field_->num_blocks();++ib)for(int id:bids)
    {
        const auto &b=field_->field(id,ib);if(!b.is_allocated())continue;
        const Int3 lo=b.get_lo(),hi=b.get_hi();
        count+=std::size_t(hi.i-lo.i)*std::size_t(hi.j-lo.j)*std::size_t(hi.k-lo.k);
    }
    std::vector<double> cache_values(count);
    const int rank=param_->GetInt("myid");int ranks=0;MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    std::ostringstream cache_name;cache_name<<"DATA/Badd_Dipoles/dipole"<<std::setw(4)
        <<std::setfill('0')<<rank<<".bin";
    const bool read_cache=param_->HasBoo("Read_Badd_Dipoles") && param_->GetBoo("Read_Badd_Dipoles");
    const std::string path=param_->HasStr("Crustal_Dipoles_File")
        ?param_->GetStr("Crustal_Dipoles_File"):"CASE/setup/Dipoles";
    // Rank-local filesystem/evaluation failures must terminate all ranks before
    // any peer enters the collective edge-owner exchange.
    try
    {
        std::ifstream input(path);
        if(!input) throw std::runtime_error("cannot open "+path);
        const std::string source((std::istreambuf_iterator<char>(input)),{});
        std::istringstream records(source);
        model.Read(records,length_ref_,magnetic_ref_,mu0_,path);
        const double cutoff_km=param_->HasDou("Crustal_Dipoles_Cutoff_km")
            ?param_->GetDou("Crustal_Dipoles_Cutoff_km"):0.0;
        model.SetCutoff(cutoff_km*1000.0/length_ref_);
        if(rank==0)std::cout<<"[Mars Dipoles] cutoff_km="<<cutoff_km
                           <<" (0=all sources; applied to vector potential)\n";
        signature.Add(cutoff_km);
        signature.Bytes(source.data(),source.size());
        signature.Add(length_ref_);signature.Add(magnetic_ref_);signature.Add(mu0_);
        signature.Add(rank);signature.Add(ranks);
        const int revision=2;signature.Add(revision); // quadrature/CT cache revision
        // Hash physical endpoint coordinates and extents, not just file names.
        for(int ib=0;ib<field_->num_blocks();++ib)for(int axis=0;axis<3;++axis)
        {
            const auto &a=field_->field(aids[axis],ib);const bool active=a.is_allocated();
            signature.Add(active);if(!active)continue;
            const Int3 lo=a.get_lo(),hi=a.get_hi();
            for(int v:{lo.i,lo.j,lo.k,hi.i,hi.j,hi.k})signature.Add(v);
            auto &g=grid_->grids(ib);
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
            {
                const int ip=i+(axis==0),jp=j+(axis==1),kp=k+(axis==2);
                for(double v:{g.x(i,j,k),g.y(i,j,k),g.z(i,j,k),g.x(ip,jp,kp),g.y(ip,jp,kp),g.z(ip,jp,kp)})
                    signature.Add(v);
            }
        }
        // Include boundary connectivity: identical positions alone do not
        // establish an identical orientation/owner synchronization contract.
        for(const std::string name:{"inner","parallel","boundary","boundary_name"})
        {
            std::ostringstream filename;filename<<"CASE/geometry/boundary_condition/"<<name;
            if(name=="boundary_name")filename<<".fvbnd";
            else filename<<std::setw(4)<<rank<<".txt";
            std::ifstream connectivity(filename.str(),std::ios::binary);
            if(!connectivity)throw std::runtime_error("cannot fingerprint "+filename.str());
            const std::string data((std::istreambuf_iterator<char>(connectivity)),{});
            signature.Bytes(data.data(),data.size());
        }
        int hit=read_cache && MARS::ReadDipoleCache(cache_name.str(),signature.value,cache_values);
        int all_hit=0;MPI_Allreduce(&hit,&all_hit,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
        cached=all_hit!=0;
        if(rank==0)std::cout<<"[Mars Dipole cache] "<<(cached?"hit":(read_cache?"miss; recomputing":"read disabled; recomputing"))<<'\n';
        if(!cached)
        {
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
    }
    catch(const std::exception &error)
    {
        std::cerr<<"Mars crustal field rank="<<param_->GetInt("myid")<<": "<<error.what()<<'\n';
        MPI_Abort(MPI_COMM_WORLD,1);
        throw;
    }
    // Includes rotated interfaces, MPI aliases and the singular-edge classes.
    // Only topology synchronization: electric wall/farfield BCs do not apply to A.
    if(!cached)
    {
        halo_->sync_owner_alias_group("Aadd");
        halo_->sync_group("Aadd");
        for(int ib=0;ib<field_->num_blocks();++ib)
        {
            if(!field_->field(bids[0],ib).is_allocated()) continue;
            CTOperators::CurlEdgeToFace(ib,field_->field(aids[0],ib),field_->field(aids[1],ib),
                field_->field(aids[2],ib),field_->field(bids[0],ib),field_->field(bids[1],ib),
                field_->field(bids[2],ib),1.0,true);
        }
    }
    else
    {
        std::size_t n=0;
        for(int ib=0;ib<field_->num_blocks();++ib)for(int id:bids)
        {
            auto &b=field_->field(id,ib);if(!b.is_allocated())continue;
            const Int3 lo=b.get_lo(),hi=b.get_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
                b(i,j,k,0)=cache_values[n++];
        }
    }
    halo_->sync_group("Badd");
    CheckBackgroundFlux("dipoles only");
    if(param_->GetInt("myid")==0)
        std::cout<<"[Mars Dipoles] "<<model.size()<<" dipoles from "<<path<<'\n';
    if(!cached)
    {
        std::size_t n=0;
        for(int ib=0;ib<field_->num_blocks();++ib)for(int id:bids)
        {
            const auto &b=field_->field(id,ib);if(!b.is_allocated())continue;
            const Int3 lo=b.get_lo(),hi=b.get_hi();
            for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)
                cache_values[n++]=b(i,j,k,0);
        }
        if(!MARS::WriteDipoleCache(cache_name.str(),signature.value,cache_values))
            std::cerr<<"[Mars Dipole cache] could not write "<<cache_name.str()<<"; computed field remains valid\n";
    }
    double elapsed=MPI_Wtime()-started,max_elapsed=0;
    MPI_Allreduce(&elapsed,&max_elapsed,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(rank==0)std::cout<<"[Mars Dipole cache] initialization seconds="<<max_elapsed<<'\n';
    // The cache intentionally excludes IMF; InitializeMagneticBackground adds
    // the current IMF exactly once after this function returns.
}

void MarsSolver::CheckBackgroundFlux(const char *label)
{
    const int bids[3]={Badd_.xi,Badd_.eta,Badd_.zeta};
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
            if(!std::isfinite(sum) || !std::isfinite(scale))
            {
                std::cerr<<"Mars Badd: non-finite background flux during "<<label<<'\n';
                MPI_Abort(MPI_COMM_WORLD,1);
            }
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
        std::cout<<"[Mars Badd check] "<<label<<", max |sum outward flux|="<<global[0]
                 <<", relative="<<(global[1]>0?global[0]/global[1]:0)<<'\n';
}
