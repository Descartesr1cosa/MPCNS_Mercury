#pragma once

class Field;
class Param;
class Z1_Boundary;

namespace Z1
{
    void InitializeFields(Field &field, Param &param, Z1_Boundary &boundary);
}
