// Core

// Z3_Lunar
#include "LunarSolver.h"

void LunarSolver::Advance()
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

bool LunarSolver::StepOnce()
{
    Compute_Timestep();

    ZeroRHS_();

    // Conservative Euler update: fluid RHS plus CT induction RHS.
    Scheme_U_();
    AddSourceToRHS_Fluid();
    AssembleRHS_Induction_CT_();
    ApplyUpdate_Euler_();
    ApplyDensityFloor_();

    // Record and Update Runtime DATA
    {
        run_data_->dt = dt;
        run_data_->time += dt;
        run_data_->step += 1;
    }

    // Add Boundary Condition And Prepare for Next Step
    {
        lunar_bound_.Sync("Ucell");
        // Fluid is the only induction domain.  Reconcile the updated normal
        // flux to every Fluid/Solid face alias before the topology owner
        // synchronization, so a Solid owner cannot restore stale Bind.
        ReconcileFaceAliasesPreferFluid_(fid_.fid_B);
        lunar_bound_.Sync("Bface");

        UpdateDerivedFields_();
    }

    // When Stop/Output/Print
    control_.UpdateSwitches(*run_data_);
    return UpdateControlAndOutput();
}

bool LunarSolver::UpdateControlAndOutput()
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
        // Build application-owned nodal output from real quotient-mesh cell
        // incidence.  In particular, do not let Tecplot interpolation sample
        // degenerate block-local ghosts at a variable-valence singular edge.
        UpdateTecplotNodeFields_();

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
