#include "MarsCase.h"
#include "0_basic/MPI_WRAPPER.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace MARS
{
void PrepareCaseWorkdirIfNeeded(int rank)
{
    namespace fs = std::filesystem;
    if (fs::exists(fs::path("CASE") / "setup" / "filenames"))
        return;
    fs::path app;
    for (const fs::path candidate : {fs::path("Z5_Mars"), fs::path("../Z5_Mars"), fs::path(".")})
        if (fs::exists(candidate / "9999setup" / "filenames")) { app = candidate; break; }
    if (app.empty())
        throw std::runtime_error("MarsZ5: cannot locate Z5_Mars/9999setup/filenames");
    const fs::path setup = fs::absolute(app / "9999setup");
    if (rank == 0)
        fs::create_directories(app / "CASE");
    const fs::path link = app / "CASE" / "setup";
    if (rank == 0 && !fs::exists(link))
        fs::create_directory_symlink(setup, link);
    PARALLEL::mpi_barrier();
    fs::current_path(app);
    if (rank == 0)
        std::cout << "MarsZ5 CASE working directory: " << fs::current_path() << '\n';
}
}
