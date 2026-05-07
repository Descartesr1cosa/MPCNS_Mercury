#include "Z0_Config.h"

#include <filesystem>
#include <iostream>

namespace
{
    bool has_arg(int argc, char **argv, const std::string &needle)
    {
        for (int i = 1; i < argc; ++i)
            if (argv[i] && needle == argv[i])
                return true;
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

    Z0::Mode parse_mode_name(const std::string &mode)
    {
        if (mode.empty())
            return Z0::Mode::Summary;
        if (mode == "summary")
            return Z0::Mode::Summary;
        if (mode == "sync")
            return Z0::Mode::Sync;
        if (mode == "coupling")
            return Z0::Mode::Coupling;
        if (mode == "output" || mode == "io")
            return Z0::Mode::Output;
        return Z0::Mode::All;
    }
}

namespace Z0
{
    Config parse_config(int argc, char **argv)
    {
        Config cfg;
        cfg.mode = parse_mode_name(arg_value(argc, argv, "--mode"));
        if (has_arg(argc, argv, "--sync-test"))
            cfg.mode = Mode::Sync;
        if (has_arg(argc, argv, "--summary"))
            cfg.mode = Mode::Summary;
        cfg.dump_registry = has_arg(argc, argv, "--dump-registry");
        cfg.write_tecplot = (cfg.mode == Mode::Output || cfg.mode == Mode::All) &&
                            !has_arg(argc, argv, "--no-tecplot");
        return cfg;
    }

    const char *mode_name(Mode mode)
    {
        switch (mode)
        {
        case Mode::Summary:
            return "summary";
        case Mode::Sync:
            return "sync";
        case Mode::Coupling:
            return "coupling";
        case Mode::Output:
            return "output";
        case Mode::All:
            return "all";
        }
        return "all";
    }

    void use_z4_case_workdir_if_needed(int myid)
    {
        namespace fs = std::filesystem;
        if (fs::exists(fs::path("CASE") / "setup" / "filenames"))
            return;

        fs::path base = fs::current_path();
        for (int depth = 0; depth < 6; ++depth)
        {
            const fs::path z0_case = base / "Z0_CoreDebug" / "CASE" / "setup" / "filenames";
            if (fs::exists(z0_case))
            {
                fs::current_path(base / "Z0_CoreDebug");
                break;
            }

            const fs::path z4_case = base / "Z4_Mercury" / "CASE" / "setup" / "filenames";
            if (fs::exists(z4_case))
            {
                fs::current_path(base / "Z4_Mercury");
                break;
            }

            const fs::path z4_setup = base / "Z4_Mercury" / "9999setup" / "filenames";
            if (fs::exists(z4_setup))
            {
                const fs::path z0_dir = base / "Z0_CoreDebug";
                const fs::path z0_case_dir = z0_dir / "CASE";
                const fs::path z0_setup_dir = z0_case_dir / "setup";
                fs::create_directories(z0_case_dir);
                if (!fs::exists(z0_setup_dir))
                {
                    try
                    {
                        fs::create_directory_symlink(fs::path("..") / ".." / "Z4_Mercury" / "9999setup",
                                                     z0_setup_dir);
                    }
                    catch (...)
                    {
                        fs::create_directories(z0_setup_dir);
                        for (const auto &entry : fs::directory_iterator(base / "Z4_Mercury" / "9999setup"))
                            fs::copy_file(entry.path(),
                                          z0_setup_dir / entry.path().filename(),
                                          fs::copy_options::overwrite_existing);
                    }
                }
                fs::current_path(z0_dir);
                break;
            }

            if (!base.has_parent_path() || base == base.parent_path())
                break;
            base = base.parent_path();
        }

        if (myid == 0)
        {
            std::cout << "Z0_CoreDebug using CASE working directory: "
                      << fs::current_path().string() << "\n";
        }
    }
}
