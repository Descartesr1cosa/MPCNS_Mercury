#include "6_io/IOModule.h"

#include <cmath>
#include <filesystem>
#include <cstdio>
#include <string>

void IOModule::CopyCurrentCheckpointToDirectory_(const std::string &dst_dir) const
{
    namespace fs = std::filesystem;

    // Let every rank make the idempotent call.  On a multi-node shared
    // filesystem, an MPI barrier does not by itself make directory metadata
    // created by rank 0 immediately visible on the other nodes.
    fs::create_directories(dst_dir);

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
    // Keep the previous backup directory and overwrite each rank's checkpoint
    // in place.  Removing and immediately recreating the same shared directory
    // can leave other compute nodes with stale directory metadata.
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

bool IOModule::IsArchiveDue(double time)
{
    if (archive_output_by_phy_time_)
    {
        if (archive_output_phy_interval_ <= 0.0 || time < next_archive_phy_time_)
            return false;

        return true;
    }

    if (archive_output_interval_s_ <= 0.0)
        return false;

    int do_archive = 0;

    if (myid_ == 0)
    {
        if (MPI_Wtime() - last_archive_wall_s_ >= archive_output_interval_s_)
            do_archive = 1;
    }

    PARALLEL::mpi_bcast(&do_archive, 1, 0);
    return do_archive != 0;
}

std::string IOModule::ArchiveDataDirectory(int step, double time)
{
    if (archive_output_by_phy_time_ && archive_output_phy_interval_ <= 0.0)
        return {};

    namespace fs = std::filesystem;
    const std::string archive_dir = MakeArchiveDirectoryName_(step, time);

    // Every rank creates the directory. On multi-node shared filesystems this
    // avoids relying on rank 0 directory metadata becoming visible
    // immediately. The barrier ensures the directory exists before any rank
    // opens its archive file.
    fs::create_directories(archive_dir);
    PARALLEL::mpi_barrier();

    const fs::path archive_restart =
        fs::path(archive_dir) / fs::path(restart_path_).filename();
    WriteRestartBinFileToPath_(step, time, archive_restart.string());

    if (myid_ == 0)
    {
        const fs::path archive_rundata =
            fs::path(archive_dir) / fs::path(run_.rundata_path_).filename();
        run_.WriteBinary(archive_rundata.string());
    }

    // Do not report completion or advance the schedule until all rank files
    // and RunData.bin are closed and visible.
    PARALLEL::mpi_barrier();

    if (archive_output_by_phy_time_)
    {
        // Keep a fixed nondimensional-time cadence even when a solver step or
        // checkpoint interval crosses one or more archive thresholds.
        const double intervals_crossed =
            std::floor((time - next_archive_phy_time_) /
                       archive_output_phy_interval_) + 1.0;
        next_archive_phy_time_ +=
            intervals_crossed * archive_output_phy_interval_;
    }
    else
    {
        last_archive_wall_s_ = MPI_Wtime();
    }

    if (myid_ == 0)
    {
        std::printf("[IOModule] Archive current DATA to %s\n", archive_dir.c_str());
        std::fflush(stdout);
    }

    return archive_dir;
}

void IOModule::CopyCheckpointToCurrentData(const std::string &src_dir)
{
    namespace fs = std::filesystem;

    const fs::path src_restart =
        fs::path(src_dir) / fs::path(restart_path_).filename();
    fs::copy_file(src_restart, restart_path_,
                  fs::copy_options::overwrite_existing);

    if (myid_ == 0)
    {
        const fs::path src_rundata =
            fs::path(src_dir) / fs::path(run_.rundata_path_).filename();
        fs::copy_file(src_rundata, run_.rundata_path_,
                      fs::copy_options::overwrite_existing);
    }

    PARALLEL::mpi_barrier();
}

void IOModule::MaybeArchiveDataDirectory(int step, double time)
{
    if (IsArchiveDue(time))
        (void)ArchiveDataDirectory(step, time);
}

void IOModule::SetArchivePhysicalTimeOrigin(double time)
{
    if (archive_output_by_phy_time_)
        next_archive_phy_time_ = time + archive_output_phy_interval_;
}
