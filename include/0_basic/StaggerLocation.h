#pragma once

enum class StaggerLocation
{
    Cell, // cell-centered DOF
    Node, // node / vertex DOF

    FaceXi, // face normal to xi
    FaceEt, // face normal to eta
    FaceZe, // face normal to zeta

    EdgeXi, // edge tangent to xi
    EdgeEt, // edge tangent to eta
    EdgeZe  // edge tangent to zeta
};

enum class HaloLevel : int
{
    Corner1D = 1, // 1D corner/face halo
    Corner2D = 2, // 2D corner halo
    Corner3D = 3  // 3D corner halo
};
