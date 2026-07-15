#include "LunarGridValidation.h"

#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace LUNAR
{
void ValidateCubicSphereGridOrAbort(Grid &grid, int nghost, int myid)
{
    if (nghost < 1)
        throw std::runtime_error("LunarZ3 cubic-sphere grids require at least one ghost layer");

    double min_physical_j = std::numeric_limits<double>::infinity();
    double min_ghost_j = std::numeric_limits<double>::infinity();
    std::size_t physical_checked = 0, ghost_checked = 0, ghost_bad = 0;
    std::string first_bad_ghost;
    for (int ib = 0; ib < grid.nblock; ++ib)
    {
        Block &b = grid.grids(ib);
        for (int k = -nghost; k <= b.mz + nghost; ++k)
            for (int j = -nghost; j <= b.my + nghost; ++j)
                for (int i = -nghost; i <= b.mx + nghost; ++i)
                {
                    const double jac = b.jacobi(i, j, k);
                    const bool physical = i >= 0 && i <= b.mx &&
                                          j >= 0 && j <= b.my &&
                                          k >= 0 && k <= b.mz;
                    const bool bad = !std::isfinite(jac) || jac == 0.0;
                    if (physical && bad)
                    {
                        std::ostringstream msg;
                        msg << "LunarZ3 degenerate cubic-sphere metric: rank=" << myid
                            << " block=" << ib << " index=(" << i << ',' << j << ',' << k
                            << ") jacobi=" << jac
                            << ". The original physical grid is degenerate.";
                        throw std::runtime_error(msg.str());
                    }
                    if (physical)
                    {
                        min_physical_j = std::min(min_physical_j, std::abs(jac));
                        ++physical_checked;
                    }
                    else
                    {
                        ++ghost_checked;
                        if (bad)
                        {
                            ++ghost_bad;
                            if (first_bad_ghost.empty())
                            {
                                std::ostringstream msg;
                                msg << "rank=" << myid << " block=" << ib << " index=("
                                    << i << ',' << j << ',' << k << ") jacobi=" << jac;
                                first_bad_ghost = msg.str();
                            }
                        }
                        else
                            min_ghost_j = std::min(min_ghost_j, std::abs(jac));
                    }
                }
    }

    double global_min_physical = 0.0, global_min_ghost = 0.0;
    double local_counts[3] = {static_cast<double>(physical_checked),
                              static_cast<double>(ghost_checked),
                              static_cast<double>(ghost_bad)};
    double global_counts[3] = {0.0, 0.0, 0.0};
    PARALLEL::mpi_min(&min_physical_j, &global_min_physical, 1);
    PARALLEL::mpi_min(&min_ghost_j, &global_min_ghost, 1);
    PARALLEL::mpi_sum(local_counts, global_counts, 3);

    if (!first_bad_ghost.empty())
        std::cerr << "Cubic-sphere ghost metric warning: " << first_bad_ghost << '\n';
    if (myid == 0)
        std::cout << "Cubic-sphere grid validation: physical=" << global_counts[0]
                  << " min physical |J|=" << global_min_physical
                  << "; ghost=" << global_counts[1]
                  << " bad ghost=" << global_counts[2]
                  << " min valid ghost |J|=" << global_min_ghost << '\n';
}
}
