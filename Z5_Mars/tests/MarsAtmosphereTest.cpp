#include "MarsAtmosphereModel.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace
{
std::array<double, 3> at(double radius, double altitude_km,
                         double x, double y, double z)
{
    const double n = std::sqrt(x*x+y*y+z*z);
    const double r = radius + altitude_km*1.0e3;
    return {r*x/n,r*y/n,r*z/n};
}
}

int main()
{
    MARS::AtmosphereModel model;
    const double radius=model.config().mars_radius_m;
    const auto n100=model.NeutralDensity(100.0);
    const auto n300=model.NeutralDensity(300.0);
    assert(n100[MARS::CO2]>n300[MARS::CO2]);
    assert(n100[MARS::N2]>n300[MARS::N2]);
    assert(n100[MARS::O]>0.0 && n300[MARS::H]>0.0);
    for(double v:n300) assert(std::isfinite(v) && v>=0.0);

    const auto day=model.Evaluate(at(radius,200.0,1.0,0.0,0.0));
    const auto night=model.Evaluate(at(radius,200.0,-1.0,0.0,0.0));
    double qday=0.0,qnight=0.0;
    for(double v:day.photo_cm3_s) qday+=v;
    for(double v:night.photo_cm3_s) qnight+=v;
    assert(qday>0.0);
    assert(qnight==0.0);

    // The spherical ray integral must remain regular at the terminator.
    const auto terminator=at(radius,200.0,0.0,1.0,0.0);
    const double tau=model.OpticalDepth(terminator);
    assert(std::isfinite(tau) && tau>=0.0);

    const auto wall=model.WallIonDensity(at(radius,100.0,1.0,0.0,0.0));
    std::cout<<"tau100="<<model.OpticalDepth(at(radius,100.0,1.0,0.0,0.0))
             <<" wall converged="<<wall.converged<<" residual="<<wall.relative_residual
             <<" iterations="<<wall.iterations<<'\n';
    assert(wall.converged);
    assert(wall.relative_residual<1.0e-6);
    for(double v:wall.ion_cm3) assert(std::isfinite(v) && v>0.0);
    std::cout<<"wall ions:";
    for(double v:wall.ion_cm3) std::cout<<' '<<v;
    std::cout<<'\n';
    const auto dark_wall=model.WallIonDensity(at(radius,100.0,-1.0,0.0,0.0));
    assert(dark_wall.converged && dark_wall.iterations==0);
    const auto terminator_wall=model.WallIonDensity(at(radius,100.0,0.0,1.0,0.0));
    const auto terminator_state=model.Evaluate(at(radius,100.0,0.0,1.0,0.0));
    double q_terminator_wall=0.0;
    for(double v:terminator_state.photo_cm3_s) q_terminator_wall+=v;
    std::cout<<"q_terminator_wall="<<q_terminator_wall
             <<" residual="<<terminator_wall.relative_residual<<'\n';
    assert(terminator_wall.converged);

    std::cout<<"q_day="<<qday<<" tau_terminator="<<tau
             <<" pce_residual="<<wall.relative_residual
             <<" pce_iterations="<<wall.iterations<<'\n';
    return 0;
}
