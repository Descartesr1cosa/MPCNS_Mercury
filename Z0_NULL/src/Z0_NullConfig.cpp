#include "Z0_NullConfig.h"

#include <filesystem>
#include <iostream>

namespace
{
    bool has_arg(int argc, char **argv, const std::string &needle)
    {
        for (int i = 1; i < argc; ++i)
        {
            if (argv[i] && needle == argv[i])
                return true;
        }
        return false;
    }
}

namespace Z0_NULL
{
    NullConfig parse_config(int argc, char **argv)
    {
        NullConfig cfg;
        cfg.dump_registry = has_arg(argc, argv, "--dump-registry");
        cfg.sync_test = has_arg(argc, argv, "--sync-test");
        cfg.write_tecplot = !has_arg(argc, argv, "--no-tecplot");
        return cfg;
    }

    void use_z4_case_workdir_if_needed(int myid)
    {
        namespace fs = std::filesystem;

        if (fs::exists(fs::path("CASE") / "setup" / "filenames"))
            return;

        const fs::path z4_case = fs::path("..") / "Z4_Mercury" / "CASE" / "setup" / "filenames";
        if (!fs::exists(z4_case))
            return;

        fs::current_path(fs::path("..") / "Z4_Mercury");
        if (myid == 0)
        {
            // Z0_NULL is physics-null. It may reuse Z4_Mercury/CASE as a
            // convenient multi-block Mercury grid input, but it does not reuse
            // Mercury field catalog or solver.
            std::cout << "Z0_NULL using Z4_Mercury CASE working directory: "
                      << fs::current_path().string() << "\n"
                      << std::flush;
        }
    }
}
