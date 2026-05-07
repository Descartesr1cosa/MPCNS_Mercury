#pragma once

#include <string>

namespace Z0_NULL
{
    struct NullConfig
    {
        bool dump_registry = false;
        bool sync_test = false;
        bool write_tecplot = true;

        int output_step = 0;
        double output_time = 0.0;

        std::string title = "Z0_NULL Synchronization Framework Example";
    };

    NullConfig parse_config(int argc, char **argv);

    void use_z4_case_workdir_if_needed(int myid);
}
