#pragma once

class Field;
class Halo;
class Param;
class Z0_Boundary;

namespace Z0_TEST
{
    bool RunHaloCommunicationTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param);
    bool RunDecChainTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param);
    bool RunPhysicalIoTests(Field &field, Param &param);
}
