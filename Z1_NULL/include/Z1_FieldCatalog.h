#pragma once

class Field;
class Halo;
class Param;

namespace TOPO
{
    struct Topology;
}

namespace Z1
{
    void RegisterFields(Field &field, Param &param, int nghost);
    void RegisterCouplingChannels(Field &field,
                                  const TOPO::Topology &topology,
                                  Param &param,
                                  int nghost,
                                  int dimension);
    void RegisterHaloFields(Field &field, Halo &halo);
}
