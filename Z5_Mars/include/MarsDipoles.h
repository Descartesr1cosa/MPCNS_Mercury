#pragma once

#include <array>
#include <istream>
#include <string>
#include <vector>

namespace MARS
{
// Input is Cartesian SI in the simulation frame. Internally A/(B_ref L_ref).
class DipoleModel
{
public:
    using Vec = std::array<double, 3>;
    void Read(std::istream &input, double length_ref, double magnetic_ref, double mu0,
              const std::string &source = "Dipoles");
    // Euclidean source-to-evaluation distance in grid units; zero disables cutoff.
    void SetCutoff(double distance_nd);
    Vec Potential(const Vec &position_nd) const;
    double EdgeIntegral(const Vec &begin_nd, const Vec &end_nd) const;
    std::size_t size() const { return dipoles_.size(); }
private:
    struct Dipole { Vec position, moment; double core2; };
    std::vector<Dipole> dipoles_;
    double cutoff2_ = 0;
};
}
