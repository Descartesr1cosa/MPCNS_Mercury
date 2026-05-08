#pragma once

#include <string>

namespace Z0
{
    enum class Mode
    {
        Summary,
        Sync,
        Coupling,
        Output,
        All
    };

    struct Config
    {
        Mode mode = Mode::All;
        bool dump_registry = false;
        bool write_tecplot = false;
        int output_step = 0;
        double output_time = 0.0;
        std::string title = "Z0_CoreDebug Framework Test";
    };

    Config parse_config(int argc, char **argv);
    const char *mode_name(Mode mode);
    void use_z4_case_workdir_if_needed(int myid);
}
