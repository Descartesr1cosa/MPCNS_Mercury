#include "Z1_Control.h"

#include "Z1_Const.h"
#include "0_basic/1_MPCNS_Parameter.h"

#include <filesystem>
#include <iostream>

namespace Z1
{
    bool Control::should_stop() const
    {
        return step >= max_step || time >= final_time;
    }

    void Control::advance()
    {
        ++step;
        time += dt;
    }

    Control LoadControl(Param &param)
    {
        (void)param;
        Control c;
        c.dt = DefaultDt;
        c.max_step = DefaultMaxStep;
        c.final_time = DefaultFinalTime;
        return c;
    }

    void PrepareCaseWorkdirIfNeeded(int myid)
    {
        namespace fs = std::filesystem;
        if (fs::exists(fs::path("CASE") / "setup" / "filenames"))
            return;

        if (fs::exists(fs::path("Z1_NULL") / "CASE" / "setup" / "filenames"))
        {
            fs::current_path("Z1_NULL");
        }
        else if (fs::exists(fs::path("Z4_Mercury") / "CASE" / "setup" / "filenames"))
        {
            fs::current_path("Z4_Mercury");
        }
        else if (fs::exists(fs::path("Z4_Mercury") / "9999setup" / "filenames"))
        {
            fs::create_directories(fs::path("Z1_NULL") / "CASE");
            if (!fs::exists(fs::path("Z1_NULL") / "CASE" / "setup"))
            {
                try
                {
                    fs::create_directory_symlink(fs::path("..") / ".." / "Z4_Mercury" / "9999setup",
                                                 fs::path("Z1_NULL") / "CASE" / "setup");
                }
                catch (...)
                {
                    fs::create_directories(fs::path("Z1_NULL") / "CASE" / "setup");
                    for (const auto &entry : fs::directory_iterator(fs::path("Z4_Mercury") / "9999setup"))
                        fs::copy_file(entry.path(),
                                      fs::path("Z1_NULL") / "CASE" / "setup" / entry.path().filename(),
                                      fs::copy_options::overwrite_existing);
                }
            }
            fs::current_path("Z1_NULL");
        }
        else if (fs::exists(fs::path("..") / "Z4_Mercury" / "CASE" / "setup" / "filenames"))
        {
            fs::current_path(fs::path("..") / "Z4_Mercury");
        }

        if (myid == 0)
            std::cout << "Z1_NULL CASE working directory: " << fs::current_path().string() << "\n";
    }
}
