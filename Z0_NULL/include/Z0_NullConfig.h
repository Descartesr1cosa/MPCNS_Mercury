#pragma once

#include <string>

namespace Z0_NULL
{
    enum class NullMode
    {
        Summary,
        Sync,
        IO,
        All
    };

    struct NullConfig
    {
        NullMode mode = NullMode::Summary;
        bool dump_registry = false;
        bool write_tecplot = false;

        int output_step = 0;
        double output_time = 0.0;

        std::string title = "Z0_NULL Synchronization Framework Example";
    };

    NullConfig parse_config(int argc, char **argv);

    const char *mode_name(NullMode mode);

    void use_z4_case_workdir_if_needed(int myid);
}
