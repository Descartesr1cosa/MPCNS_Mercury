#pragma once

class Grid;

namespace LUNAR
{
    void ValidateCubicSphereGridOrAbort(Grid &grid, int nghost, int myid);
}
