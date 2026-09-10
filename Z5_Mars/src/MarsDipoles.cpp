#include "MarsDipoles.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace MARS
{
void DipoleModel::Read(std::istream &input, double L, double B, double mu0,
                       const std::string &source)
{
    if (!(L > 0 && B > 0 && mu0 > 0) || !std::isfinite(L+B+mu0))
        throw std::runtime_error(source + ": invalid reference scales");
    std::vector<Dipole> parsed;
    std::string line;
    std::size_t line_number=0, count=0;
    bool header=false;
    const double scale=mu0/(4.0*std::acos(-1.0)*B*L*L*L);
    while(std::getline(input,line))
    {
        ++line_number;
        line=line.substr(0,line.find('#'));
        std::istringstream row(line);
        row >> std::ws;
        if(row.eof()) continue;
        auto fail=[&]() { throw std::runtime_error(source+":"+std::to_string(line_number)+
                                                  ": invalid Dipoles record"); };
        std::string extra;
        if(!header)
        {
            std::string magic; long long n=-1;
            if(!(row>>magic>>n) || magic!="MARS_DIPOLES_V1" || n<0 || n>1000000 || row>>extra)
                fail();
            count=static_cast<std::size_t>(n); parsed.reserve(count); header=true;
            continue;
        }
        Vec r{},m{}; double core=0;
        if(!(row>>r[0]>>r[1]>>r[2]>>m[0]>>m[1]>>m[2]>>core) || row>>extra)
            fail();
        if(parsed.size()>=count || core<0 || !std::isfinite(core)) fail();
        for(int d=0;d<3;++d)
        {
            if(!std::isfinite(r[d]) || !std::isfinite(m[d])) fail();
            r[d]/=L; m[d]*=scale;
            if(!std::isfinite(r[d]) || !std::isfinite(m[d])) fail();
        }
        core/=L;
        if(!std::isfinite(core*core)) fail();
        parsed.push_back({r,m,core*core});
    }
    if(input.bad() || !header || parsed.size()!=count)
        throw std::runtime_error(source+": missing header, unreadable input or dipole count mismatch");
    dipoles_=std::move(parsed);
}

void DipoleModel::SetCutoff(double distance_nd)
{
    if(distance_nd<0 || !std::isfinite(distance_nd) || !std::isfinite(distance_nd*distance_nd))
        throw std::runtime_error("Dipoles: cutoff must be finite and nonnegative");
    cutoff2_=distance_nd*distance_nd;
}

DipoleModel::Vec DipoleModel::Potential(const Vec &position) const
{
    Vec sum{},correction{};
    for(const auto &dip:dipoles_)
    {
        Vec r{}; double distance2=0;
        for(int d=0;d<3;++d) { r[d]=position[d]-dip.position[d]; distance2+=r[d]*r[d]; }
        if(!std::isfinite(distance2))
            throw std::runtime_error("Dipoles: non-finite evaluation distance");
        // Cut only A, before the shared-edge synchronization and CT curl.
        // Cutting face B afterward would violate the discrete divergence constraint.
        if(cutoff2_>0 && distance2>cutoff2_) continue;
        const double r2=distance2+dip.core2;
        if(!(r2>0) || !std::isfinite(r2))
            throw std::runtime_error("Dipoles: singular/non-finite evaluation point");
        const double inv=1.0/(r2*std::sqrt(r2));
        for(int d=0;d<3;++d)
        {
            const int j=(d+1)%3,k=(d+2)%3;
            const double term=(dip.moment[j]*r[k]-dip.moment[k]*r[j])*inv;
            // Compensated summation matters for thousands of cancelling sources.
            const double y=term-correction[d],t=sum[d]+y;
            correction[d]=(t-sum[d])-y; sum[d]=t;
        }
    }
    for(double v:sum) if(!std::isfinite(v))
        throw std::runtime_error("Dipoles: non-finite vector potential");
    return sum;
}

double DipoleModel::EdgeIntegral(const Vec &begin, const Vec &end) const
{
    Vec mid{},dr{};
    for(int d=0;d<3;++d) {mid[d]=0.5*(begin[d]+end[d]);dr[d]=end[d]-begin[d];}
    if(dr[0]==0 && dr[1]==0 && dr[2]==0) return 0;
    // Symmetric two-point Gauss quadrature on the straight grid edge.
    const double q=0.5/std::sqrt(3.0);
    Vec p{},m{};
    for(int d=0;d<3;++d) {p[d]=mid[d]+q*dr[d];m[d]=mid[d]-q*dr[d];}
    const Vec ap=Potential(p),am=Potential(m);
    double value=0;
    for(int d=0;d<3;++d) value+=0.5*(ap[d]+am[d])*dr[d];
    if(!std::isfinite(value)) throw std::runtime_error("Dipoles: non-finite edge integral");
    return value;
}
}
