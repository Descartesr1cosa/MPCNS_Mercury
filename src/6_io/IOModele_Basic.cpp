#include "6_io/IOModule.h"

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <filesystem>

void IOModule::Setup(Param *par, Grid *grd, Field *fld, int nvar)
{
    par_ = par;
    grd_ = grd;
    fld_ = fld;
    if (!par_ || !fld_ || !grd_)
        Fail_("[IOModule] Setup: null par/fld/grd");

    if (par_->GetInt("myid") == 0)
    {
        std::filesystem::create_directories("./DATA");
        // std::filesystem::create_directories("./DATA_backup");
        std::filesystem::create_directories("./DATA_archive");
    }
    PARALLEL::mpi_barrier();

    // 获取bin文件的输出路径
    const int myid = par_->GetInt("myid");

    archive_output_by_phy_time_ =
        par_->HasBoo("Archive_Output_By_PhyTime") &&
        par_->GetBoo("Archive_Output_By_PhyTime");

    const double archive_hours = par_->GetDou("Archive_Output_Time");
    if (archive_hours > 0.0)
        archive_output_interval_s_ = archive_hours * 3600.0;
    else
        archive_output_interval_s_ = 0.0;
    archive_output_phy_interval_ =
        par_->HasDou("Archive_Output_PhyTime")
            ? par_->GetDou("Archive_Output_PhyTime")
            : 0.0;
    if (archive_output_by_phy_time_ && archive_output_phy_interval_ <= 0.0)
        Fail_("[IOModule] Archive_Output_PhyTime must be positive when "
              "Archive_Output_By_PhyTime is enabled");
    next_archive_phy_time_ = archive_output_phy_interval_;

    last_archive_wall_s_ = MPI_Wtime();
    if (myid == 0)
    {
        if (archive_output_by_phy_time_)
            std::printf("[IOModule] Archive mode = physic_time, interval = %.6e nondimensional time\n",
                        archive_output_phy_interval_);
        else
            std::printf("[IOModule] Archive mode = CPU computing time, "
                        "interval = %.6e hours (%.6e seconds)\n",
                        archive_hours, archive_output_interval_s_);
        std::fflush(stdout);
    }

    auto rank4 = [](int id) -> std::string
    {
        char buf[8];
        std::snprintf(buf, sizeof(buf), "%04d", id);
        return std::string(buf);
    };
    restart_path_ = "./DATA/flow_field" + rank4(myid) + ".bin";
    tecplot_path_ = "./DATA/flow_field" + rank4(myid) + ".plt";

    myid_ = par->GetInt("myid");

    // residual_ref 长度由 solver 给定（nvar）
    run_.ResizeResidualRef(nvar, 1.0);
}

void IOModule::Fail_(const std::string &msg)
{
    std::fprintf(stderr, "%s\n", msg.c_str());
    std::abort();
}
