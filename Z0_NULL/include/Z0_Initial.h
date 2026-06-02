#pragma once

class Field;
class Param;
class Z0_Boundary;

namespace Z0
{
    void InitializeFields(Field &field, Param &param, Z0_Boundary &boundary);
}
