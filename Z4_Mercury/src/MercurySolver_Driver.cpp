// Core

// Z4_Mercury
#include "MercurySolver.h"

void MercurySolver::Advance()
{
    // 初始先算一遍派生量（用于输出/诊断）
    UpdateDerivedFields_();

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
    Compute_Timestep();

    ZeroRHS_();

    // Conservative Euler update: fluid RHS plus CT induction RHS.
    Scheme_U_();
    AddSourceToRHS_Fluid();
    AssembleRHS_Induction_CT_();
    ApplyUpdate_Euler_();

    // Mercury internal resistive correction acts only on B_face.
    if (resist_control.is_Mercury_resistance)
    {
        mercury_bound_.Sync("Bface");
        SolveImplicitResistiveDiffusion_(dt);
    }

    // Record and Update Runtime DATA
    {
        run_data_->dt = dt;
        run_data_->time += dt;
        run_data_->step += 1;
    }

    // Add Boundary Condition And Prepare for Next Step
    {
        mercury_bound_.Sync("Ucell");
        // The implicit update already finishes with a full Bface sync.  Avoid
        // immediately repeating the identical halo/owner exchange.
        if (!resist_control.is_Mercury_resistance)
            mercury_bound_.Sync("Bface");

        UpdateDerivedFields_();
    }

    // When Stop/Output/Print
    control_.UpdateSwitches(*run_data_);
    return UpdateControlAndOutput();
}

bool MercurySolver::UpdateControlAndOutput()
{
    RunData &run = io_.Run();
    const bool archive_due = io_.IsArchiveDue(run.time);

    if (control_.if_outres)
    {
        calc_divB();
        runtime_data_->UpdateOnOutres(run);
        runtime_data_->PrintLineMinimal(run);
        PrintMinMaxDiagnostics_();
    }

    // Archive channel: serialize directly into a new DATA_archive directory.
    // It never reads or writes DATA or DATA_backup and does not produce .plt.
    std::string archive_dir;
    if (archive_due)
        archive_dir = io_.ArchiveDataDirectory(run.step, run.time);

    // Normal output channel: write DATA, then mirror the newly completed
    // restart checkpoint to DATA_backup. If both channels are due, reuse the
    // archive checkpoint so the field data are serialized only once.
    if (control_.if_outfile)
    {
        io_.WriteTecplotBinFile(run.step, run.time);
        if (archive_due)
        {
            io_.CopyCheckpointToCurrentData(archive_dir);
        }
        else
        {
            io_.WriteRestartBinFile(run.step, run.time);
            io_.WriteRunDataFile();
        }
        io_.BackupDataDirectory("./DATA_backup");
    }

    return control_.if_stop;
}
