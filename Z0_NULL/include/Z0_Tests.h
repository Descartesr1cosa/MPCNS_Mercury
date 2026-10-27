#pragma once

class Field;
class Grid;
class Halo;
class Param;
class Z0_Boundary;
namespace TOPO
{
    struct Topology;
}

namespace Z0_TEST
{
    bool RunTopologyTests(const TOPO::Topology &topology, int myid);
    bool RunMeshPartitionTests(const Grid &grid, const TOPO::Topology &topology, int myid);
    bool RunHaloCommunicationTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param);
    bool RunExchangePlanTests(const Grid &grid, const TOPO::Topology &topology,
                              Field &field, Halo &halo, int myid);
    bool RunDecChainTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param);
    bool RunPhysicalIoTests(Field &field, Param &param);
}
