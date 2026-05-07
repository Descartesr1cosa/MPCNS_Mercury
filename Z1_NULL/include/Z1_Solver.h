#pragma once

#include "Z1_Control.h"

class Field;
class Grid;
class Halo;
class Param;
class Z1_Boundary;

namespace TOPO
{
    struct Topology;
    struct TopologyEquiv;
}

class Z1_Solver
{
public:
    Z1_Solver(Grid *grid,
              Field *field,
              Halo *halo,
              TOPO::Topology *topology,
              TOPO::TopologyEquiv *topology_equiv,
              Z1_Boundary *boundary,
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
    TOPO::TopologyEquiv *topology_equiv_ = nullptr;
    Z1_Boundary *boundary_ = nullptr;
    Param *param_ = nullptr;

    Z1::Control control_;
    bool initialized_ = false;
    bool setup_done_ = false;
};
