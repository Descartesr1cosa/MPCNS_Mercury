#pragma once
#include <chrono>
#include <ctime>
#include <cstdint>
#include <string>
#include <algorithm>
#include <cstdio>
#include "6_io/RunData.h"

class Param;
class RuntimeMonitor
{
public:
    void Begin(RunData &run, Param *par, int64_t global_cells);
    void UpdateOnOutres(RunData &run);
    void PrintLineMinimal(const RunData &run) const;
    bool Inited() const { return inited_; }

private:
    Param *par_{nullptr};
    int64_t global_cells_{0};
    bool inited_{false};
    double cpu_last_s_{0};
    int32_t step_last_{0};
};