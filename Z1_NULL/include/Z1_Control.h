#pragma once

class Param;

namespace Z1
{
    struct Control
    {
        int step = 0;
        double time = 0.0;
        double dt = 1.0e-3;
        int max_step = 1;
        double final_time = 1.0e-3;

        bool should_stop() const;
        void advance();
    };

    Control LoadControl(Param &param);
    void PrepareCaseWorkdirIfNeeded(int myid);
}
