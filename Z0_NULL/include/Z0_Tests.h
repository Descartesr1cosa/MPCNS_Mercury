#pragma once

class Field;
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
    bool RunHaloCommunicationTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param);
    bool RunDecChainTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param);
    bool RunPhysicalIoTests(Field &field, Param &param);
}
