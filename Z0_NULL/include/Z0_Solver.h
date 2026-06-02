#pragma once

#include "Z0_Control.h"

class Field;
class Grid;
class Halo;
class Param;
class Z0_Boundary;

namespace TOPO
{
    struct Topology;
}

class Z0_Solver
{
public:
    Z0_Solver(Grid *grid,
              Field *field,
              Halo *halo,
              TOPO::Topology *topology,
              TOPO::Topology *topology_equiv,
              Z0_Boundary *boundary,
              Param *param);

    void Initialize();
    void Advance();

private:
    void Setup_();
    void PrepareStep_();
    void ComputeTimeStep_();
    void ComputeRHS_();
    void ComputeFlux_();
    void ComputeSource_();
    void UpdateFields_();
    void ApplyBoundaryAndSync_();
    void Diagnostics_();
    void Output_();
    void AdvanceRunState_();

    Grid *grid_ = nullptr;
    Field *field_ = nullptr;
    Halo *halo_ = nullptr;
    TOPO::Topology *topology_ = nullptr;
    TOPO::Topology *topology_equiv_ = nullptr;
    Z0_Boundary *boundary_ = nullptr;
    Param *param_ = nullptr;

    Z0::Control control_;
    bool initialized_ = false;
    bool setup_done_ = false;
};
