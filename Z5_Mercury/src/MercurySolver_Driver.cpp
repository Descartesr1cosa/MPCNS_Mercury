// Core

// Z4_Mercury
#include "MercurySolver.h"

void MercurySolver::Advance()
{
    // 初始先算一遍派生量（用于输出/诊断）
    // Boundary_Condition();
    calc_Bcell();
    // Calc_J_Edge();
    // calc_Jcell();
    calc_Jcell_from_Bcell_metric_();
    calc_PV();
    calc_Uplus();

    // step=0 也允许输出一次
    control_.UpdateSwitches(*run_data_);
    UpdateControlAndOutput();

    while (!control_.if_stop)
    {
        StepOnce();
    }
}

bool MercurySolver::StepOnce()
{
    // Calculate time step From CFL
    Compute_Timestep();

    // Time Advance
    Time_Advance();

    // Record and Update Runtime DATA
    {
        run_data_->dt = dt;
        run_data_->time += dt;
        run_data_->step += 1;
    }

    // Add Boundary Condition And Prepare for Next Step
    {
        mercury_bound_.Sync("Ucell");
        mercury_bound_.Sync("Bface");

        calc_Bcell();
        // Calc_J_Edge();
        // calc_Jcell();
        calc_Jcell_from_Bcell_metric_();
        calc_PV();
        calc_Uplus();
    }

    // When Stop/Output/Print
    control_.UpdateSwitches(*run_data_);
    return UpdateControlAndOutput();
}

bool MercurySolver::UpdateControlAndOutput()
{
    RunData &run = io_.Run();

    if (control_.if_outres)
    {
        calc_divB();
        runtime_data_->UpdateOnOutres(run);
        runtime_data_->PrintLineMinimal(run);
        PrintMinMaxDiagnostics_();
    }

    if (control_.if_outfile)
    {
        // 写新 DATA 前，先保护旧 DATA
        io_.BackupDataDirectory("./DATA_backup");

        // 写当前 checkpoint
        io_.WriteTecplotBinFile(run.step, run.time);
        io_.WriteRestartBinFile(run.step, run.time);
        io_.WriteRunDataFile();

        // 当前 checkpoint 写完之后，再判断是否需要长期归档
        io_.MaybeArchiveDataDirectory(run.step, run.time);
    }

    return control_.if_stop;
}