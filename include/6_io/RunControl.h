#pragma once

#include "0_basic/1_MPCNS_Parameter.h"
#include "6_io/RunData.h"

#include <cstdio>
#include <cstdlib>

class Control
{
public:
    void Setup(Param *parameter)
    {
        if(!parameter)
        {
            std::fprintf(stderr,"[Control] Setup: parameter is null\n");
            std::abort();
        }
        maximum_steps_=parameter->GetInt("max_Nstep");
        maximum_time_=parameter->GetDou("max_Time");
        file_interval_=parameter->GetInt("output_step");
        residual_interval_=parameter->GetInt("output_residual");
        if_stop=false;
        if_outres=false;
        if_outfile=false;
    }

    void UpdateSwitches(const RunData &run)
    {
        if_outres=residual_interval_>0 && run.step%residual_interval_==0;
        if_outfile=file_interval_>0 && run.step%file_interval_==0;
        if_stop=run.step>=maximum_steps_ || (maximum_time_>0.0 && run.time>=maximum_time_);
        if_outres=if_outres || if_stop;
    }

    bool if_stop=false;
    bool if_outres=false;
    bool if_outfile=false;

private:
    int maximum_steps_=0;
    int file_interval_=0;
    int residual_interval_=0;
    double maximum_time_=0.0;
};
