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

    std::string arg_value(int argc, char **argv, const std::string &name)
    {
        const std::string prefix = name + "=";
        for (int i = 1; i < argc; ++i)
        {
            if (!argv[i])
                continue;

            const std::string arg = argv[i];
            if (arg == name && i + 1 < argc && argv[i + 1])
                return argv[i + 1];
            if (arg.rfind(prefix, 0) == 0)
                return arg.substr(prefix.size());
        }
        return "";
    }

    Z0_NULL::NullMode parse_mode_name(const std::string &mode)
    {
        if (mode == "sync")
            return Z0_NULL::NullMode::Sync;
        if (mode == "io")
            return Z0_NULL::NullMode::IO;
        if (mode == "all")
            return Z0_NULL::NullMode::All;
        return Z0_NULL::NullMode::Summary;
    }
}

namespace Z0_NULL
{
    NullConfig parse_config(int argc, char **argv)
    {
        NullConfig cfg;
        cfg.mode = parse_mode_name(arg_value(argc, argv, "--mode"));
        if (has_arg(argc, argv, "--sync-test"))
            cfg.mode = NullMode::Sync;

        cfg.dump_registry = has_arg(argc, argv, "--dump-registry");
        cfg.write_tecplot = (cfg.mode == NullMode::IO || cfg.mode == NullMode::All) &&
                            !has_arg(argc, argv, "--no-tecplot");
        return cfg;
    }

    const char *mode_name(NullMode mode)
    {
        switch (mode)
        {
        case NullMode::Summary:
            return "summary";
        case NullMode::Sync:
            return "sync";
        case NullMode::IO:
            return "io";
        case NullMode::All:
            return "all";
        }
        return "summary";
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
