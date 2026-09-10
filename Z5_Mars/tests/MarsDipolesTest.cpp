#include "MarsDipoles.h"
#include "3_field/FieldBlock.h"
#include "8_dec/CTOperators.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
void check(bool condition,const char *message)
{ if(!condition) throw std::runtime_error(message); }

MARS::DipoleModel read(const std::string &data)
{
    std::istringstream in(data);MARS::DipoleModel model;
    model.Read(in,1,1,4*std::acos(-1.0));return model;
}

// Exercise the production CT operator on a genuinely skewed grid.
double ct_test(const MARS::DipoleModel &model)
{
    Block block;block.mx=4;block.my=3;block.mz=2;
    FieldBlock a[3],b[3];
    const StaggerLocation el[]={StaggerLocation::EdgeXi,StaggerLocation::EdgeEt,StaggerLocation::EdgeZe};
    const StaggerLocation fl[]={StaggerLocation::FaceXi,StaggerLocation::FaceEt,StaggerLocation::FaceZe};
    auto pos=[](int i,int j,int k)->MARS::DipoleModel::Vec
    {return {1.2+0.09*i+0.03*j,0.1+0.02*i+0.1*j+0.01*k,0.2+0.02*j+0.08*k};};
    for(int d=0;d<3;++d)
    {
        a[d].allocate(block,{"A",el[d],1,1});b[d].allocate(block,{"B",fl[d],1,1});
        auto hi=a[d].get_hi();
        for(int i=-1;i<hi.i;++i)for(int j=-1;j<hi.j;++j)for(int k=-1;k<hi.k;++k)
            a[d](i,j,k,0)=model.EdgeIntegral(pos(i,j,k),pos(i+(d==0),j+(d==1),k+(d==2)));
    }
    CTOperators::CurlEdgeToFace(0,a[0],a[1],a[2],b[0],b[1],b[2],1.0,true);
    double residual=0,scale=0;
    for(int i=-1;i<=block.mx;++i)for(int j=-1;j<=block.my;++j)for(int k=-1;k<=block.mz;++k)
    {
        const double f[]={b[0](i+1,j,k,0),-b[0](i,j,k,0),b[1](i,j+1,k,0),-b[1](i,j,k,0),b[2](i,j,k+1,0),-b[2](i,j,k,0)};
        double sum=0,s=0;for(double v:f){sum+=v;s+=std::abs(v);}
        residual=std::max(residual,std::abs(sum));scale=std::max(scale,s);
    }
    check(scale>0,"CT fixture must have nonzero magnetic flux");
    check(residual/scale<2e-14,"div(curl A) exceeds roundoff");
    return residual/scale;
}
}

int main(int argc,char **argv)
{
    auto model=read("# test\nMARS_DIPOLES_V1 1\n0 0 0 0 0 1 0 # axial\n");
    auto a=model.Potential({2,0,0});
    check(std::abs(a[1]-0.25)<1e-15 && a[0]==0 && a[2]==0,"dipole sign/scale");
    // Independent analytic curl at equator: Bz=-m/r^3.
    const double h=1e-5;
    auto xp=model.Potential({2+h,0,0}),xm=model.Potential({2-h,0,0});
    auto yp=model.Potential({2,h,0}),ym=model.Potential({2,-h,0});
    check(std::abs((xp[1]-xm[1]-yp[0]+ym[0])/(2*h)+0.125)<1e-10,"analytic B=curl A");
    const MARS::DipoleModel::Vec p={2,1,0.4},q={2.1,1.3,0.5};
    check(std::abs(model.EdgeIntegral(p,q)+model.EdgeIntegral(q,p))<1e-16,"edge orientation");
    check(model.EdgeIntegral(p,p)==0,"collapsed edge");
    auto doubled=read("MARS_DIPOLES_V1 2\n0 0 0 0 0 1 0\n0 0 0 0 0 1 0\n");
    check(std::abs(doubled.Potential(p)[0]-2*model.Potential(p)[0])<1e-16,"superposition");
    std::istringstream si("MARS_DIPOLES_V1 1\n0 0 0 0 0 80 0\n");
    MARS::DipoleModel scaled;scaled.Read(si,2,10,4*std::acos(-1.0));
    check(std::abs(scaled.Potential(p)[1]-model.Potential(p)[1])<1e-16,"SI normalization");
    for(const char *invalid:{"", "MARS_DIPOLES_V1 -1", "MARS_DIPOLES_V1 1", "MARS_DIPOLES_V1 0\n0 0 0 0 0 1 0", "MARS_DIPOLES_V1 1\n0 0 0 0 0 1 -1", "MARS_DIPOLES_V1 1\n0 0 0 0 0 nan 0", "MARS_DIPOLES_V1 1\n0 0 0 0 0 1 0 extra"})
    {
        bool rejected=false;try{read(invalid);}catch(const std::runtime_error&){rejected=true;}
        check(rejected,"invalid input accepted");
    }
    bool singular=false;try{model.Potential({0,0,0});}catch(const std::runtime_error&){singular=true;}
    check(singular,"point-source singularity accepted");
    auto core=read("MARS_DIPOLES_V1 1\n0 0 0 0 0 1 0.1");
    check(core.Potential({0,0,0})[0]==0,"softened source not finite");
    std::cout<<"single dipole CT relative flux residual="<<ct_test(model)<<'\n';
    check(argc==2,"expected example file");
    std::ifstream input(argv[1]);MARS::DipoleModel many;
    many.Read(input,3.396e6,1e-9,1.256637061e-6,argv[1]);
    check(many.size()==4840,"example source count");
    std::cout<<"4840 dipoles CT relative flux residual="<<ct_test(many)<<'\n';
}
