#include "MercuryGridValidation.h"

#include "1_grid/1_MPCNS_Grid.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace MERCURY
{
void ValidateCubicSphereGridOrAbort(Grid &grid, int nghost, int myid)
{
    if (nghost < 1)
        throw std::runtime_error("MercuryZ4 cubic-sphere grids require at least one ghost layer");

    double min_abs_j = std::numeric_limits<double>::infinity();
    std::size_t checked = 0;
    for (int ib = 0; ib < grid.nblock; ++ib)
    {
        Block &b = grid.grids(ib);
        for (int k = -nghost; k <= b.mz + nghost; ++k)
            for (int j = -nghost; j <= b.my + nghost; ++j)
                for (int i = -nghost; i <= b.mx + nghost; ++i)
                {
                    const double jac = b.jacobi(i, j, k);
                    if (!std::isfinite(jac) || jac == 0.0)
                    {
                        std::ostringstream msg;
                        msg << "MercuryZ4 degenerate cubic-sphere metric: rank=" << myid
                            << " block=" << ib << " index=(" << i << ',' << j << ',' << k
                            << ") jacobi=" << jac
                            << ". Check panel connectivity and non-degenerate ghost coordinates.";
                        throw std::runtime_error(msg.str());
                    }
                    min_abs_j = std::min(min_abs_j, std::abs(jac));
                    ++checked;
                }
    }
    if (myid == 0)
        std::cout << "Cubic-sphere grid validation: " << checked
                  << " physical/ghost metric points checked; min |J|=" << min_abs_j << '\n';
}
}
