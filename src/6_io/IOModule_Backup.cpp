#include "6_io/IOModule.h"

#include <filesystem>
#include <cstdio>
#include <string>

void IOModule::CopyCurrentCheckpointToDirectory_(const std::string &dst_dir) const
{
    namespace fs = std::filesystem;

    if (myid_ == 0)
    {
        fs::create_directories(dst_dir);
    }

    PARALLEL::mpi_barrier();

    // 每个 rank 复制自己的 restart bin
    {
        fs::path src(restart_path_);
        fs::path dst = fs::path(dst_dir) / src.filename();

        if (fs::exists(src))
        {
            fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
        }
    }

    // rank0 复制 RunData.bin
    if (myid_ == 0)
    {
        fs::path src(run_.rundata_path_);
        fs::path dst = fs::path(dst_dir) / src.filename();

        if (fs::exists(src))
        {
            fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
        }
    }

    PARALLEL::mpi_barrier();
}

void IOModule::BackupDataDirectory(const std::string &backup_dir)
{
    namespace fs = std::filesystem;

    if (myid_ == 0)
    {
        fs::remove_all(backup_dir);
        fs::create_directories(backup_dir);
    }

    PARALLEL::mpi_barrier();

    CopyCurrentCheckpointToDirectory_(backup_dir);

    if (myid_ == 0)
    {
        std::printf("\n[IOModule] Backup current DATA to %s\n\n", backup_dir.c_str());
        std::fflush(stdout);
    }
}

std::string IOModule::MakeArchiveDirectoryName_(int step, double time) const
{
    char buf[256];

    std::snprintf(buf, sizeof(buf),
                  "./DATA_archive/Step_%09d_Time_%.6e",
                  step,
                  time);

    return std::string(buf);
}

void IOModule::MaybeArchiveDataDirectory(int step, double time)
{
    if (archive_output_interval_s_ <= 0.0)
        return;

    int do_archive = 0;

    const double now = MPI_Wtime();

    if (myid_ == 0)
    {
        if (now - last_archive_wall_s_ >= archive_output_interval_s_)
            do_archive = 1;
    }

    PARALLEL::mpi_bcast(&do_archive, 1, 0);

    if (!do_archive)
        return;

    const std::string archive_dir = MakeArchiveDirectoryName_(step, time);

    CopyCurrentCheckpointToDirectory_(archive_dir);

    last_archive_wall_s_ = MPI_Wtime();

    if (myid_ == 0)
    {
        std::printf("[IOModule] Archive current DATA to %s\n", archive_dir.c_str());
        std::fflush(stdout);
    }
}
