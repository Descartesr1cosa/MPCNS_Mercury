#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "0_basic/1_MPCNS_Parameter.h"
#include "5_io/RunData.h"

class Control
{
public:
    void Setup(Param *par)
    {
        if (!par)
        {
            std::fprintf(stderr, "[Control] Setup: par is null\n");
            std::abort();
        }
        par_ = par;
        myid_ = par_->GetInt("myid");

        maximum_Nstep_ = par_->GetInt("max_Nstep");
        // tolerance_ = par_->GetDou("tolerance");
        maximum_time_ = par_->GetDou("max_Time");

        file_step_ = par_->GetInt("output_step");
        res_step_ = par_->GetInt("output_residual");

        // 默认状态
        if_stop = false;
        if_outres = false;
        if_outfile = false;
    }

    void UpdateSwitches(const RunData &run_data)
    {
        // ---------- output switches ----------
        // 残差/监测打印
        if_outres = (res_step_ > 0) && (run_data.step % res_step_ == 0);

        // 文件输出：一般每个 rank 都要写自己的 flow_field####.bin / zone 数据
        if_outfile = (file_step_ > 0) && (run_data.step % file_step_ == 0);

        // ---------- stop switches ----------
        if_stop = false;

        // max steps
        if (run_data.step >= maximum_Nstep_)
            if_stop = true;

        // max time
        if (maximum_time_ > 0.0 && run_data.time >= maximum_time_)
            if_stop = true;
        // if_outfile = if_outfile || if_stop; // when program stops, it will not output bin file
        if_outres = if_outres || if_stop;
    }

public:
    bool if_stop = false;
    bool if_outres = false;
    bool if_outfile = false;

private:
    Param *par_ = nullptr;
    int myid_ = 0;

    int maximum_Nstep_ = 0;
    int file_step_ = 0;
    int res_step_ = 0;

    // double tolerance_ = 0.0;
    double maximum_time_ = 0.0;
};