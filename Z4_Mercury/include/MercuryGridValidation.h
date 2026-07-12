#pragma once

class Grid;

namespace MERCURY
{
    void ValidateCubicSphereGridOrAbort(Grid &grid, int nghost, int myid);
}
