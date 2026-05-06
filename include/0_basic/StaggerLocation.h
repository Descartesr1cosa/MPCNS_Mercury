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
    FaceOnly = 1, // 只需要 1D（面）halo
    Edge = 2,     // 需要到 2D corner（棱）
    Vertex = 3    // 需要到 3D corner（角点）
};