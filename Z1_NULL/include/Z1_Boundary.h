#pragma once

#include <string>

class Field;
class Grid;
class Halo;
class Param;

namespace TOPO
{
    struct Topology;
}

class Z1_Boundary
{
public:
    Z1_Boundary(Grid *grid,
                Field *field,
                Halo *halo,
                TOPO::Topology *topology,
                TOPO::Topology *topology_equiv,
                Param *param);

    void ApplyPhysicalBoundary(const std::string &group);
    void ApplyAllPhysicalBoundaries();

    void ApplyCoupling();

    void SyncGroup(const std::string &group);
    void SyncAllRegistered();

    void PrepareStepBoundary();
    void FinishStepBoundary();

private:
    Grid *grid_ = nullptr;
    Field *field_ = nullptr;
    Halo *halo_ = nullptr;
    TOPO::Topology *topology_ = nullptr;
    TOPO::Topology *topology_equiv_ = nullptr;
    Param *param_ = nullptr;
};
