#include "MercuryCase.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace MERCURY
{
void PrepareCaseWorkdirIfNeeded(int myid)
{
    namespace fs = std::filesystem;
    if (fs::exists(fs::path("CASE") / "setup" / "filenames"))
        return;

    fs::path application;
    if (fs::exists(fs::path("Z4_Mercury") / "9999setup" / "filenames"))
        application = "Z4_Mercury";
    else if (fs::exists(fs::path("..") / "Z4_Mercury" / "9999setup" / "filenames"))
        application = fs::path("..") / "Z4_Mercury";
    else if (fs::exists(fs::path("9999setup") / "filenames"))
        application = ".";
    else
        throw std::runtime_error("MercuryZ4: cannot locate Z4_Mercury/9999setup/filenames");

    const fs::path setup = fs::absolute(application / "9999setup");
    const fs::path case_dir = application / "CASE";
    const fs::path case_setup = case_dir / "setup";
    fs::create_directories(case_dir);
    if (!fs::exists(case_setup))
    {
        try
        {
            fs::create_directory_symlink(setup, case_setup);
        }
        catch (const fs::filesystem_error &)
        {
            fs::create_directories(case_setup);
            for (const auto &entry : fs::directory_iterator(setup))
                if (entry.is_regular_file())
                    fs::copy_file(entry.path(), case_setup / entry.path().filename(),
                                  fs::copy_options::overwrite_existing);
        }
    }
    fs::current_path(application);
    if (myid == 0)
        std::cout << "MercuryZ4 CASE working directory: " << fs::current_path() << '\n';
}
}
