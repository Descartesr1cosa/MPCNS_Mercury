#include "6_io/IOModule.h"

#include <cstdio>
#include <cstdlib>

void IOModule::ReadRunDataFile()
{
    // expected_nref 用 run_.residual_ref.size() 做一致性检查
    int expected_nref = (int)run_.residual_ref.size();
    run_.ReadBinary(expected_nref);
}

void IOModule::WriteRunDataFile()
{
    if (myid_ == 0)
        run_.WriteBinary();
}
