// Core

// Z4_Mercury
#include "MercurySolver.h"

#include <algorithm>
#include <cmath>
#include <iostream>

void MercurySolver::Advance()
{
    // Initial derived fields for diagnostics/output.
    calc_Bcell();
    calc_Jcell_from_Bcell_metric_();
    calc_PV();
    calc_Uplus();

    control_.UpdateSwitches(*run_data_);
    UpdateControlAndOutput();

    while (!control_.if_stop)
    {
        Compute_Timestep();

#if HALL_IMPLICIT == 1
        double Emag0 = 0.0;
        if (control_.if_outres)
            Emag0 = ComputeMagEnergy_Cell_();
#endif

        ZeroRHS_();

        Scheme_U_();
        AddSourceToRHS_Fluid();

        AssembleRHS_Induction_CT_();
        ApplyUpdate_Euler_();

#if HALL_IMPLICIT == 1
        hall_implicit_.SolveOneStep(dt, control_.if_outres);
        calc_Bcell();

        if (control_.if_outres)
        {
            const double Emag1 = ComputeMagEnergy_Cell_();
            if (par_->GetInt("myid") == 0)
            {
                const double dE = Emag1 - Emag0;
                const double rel = dE / std::max(std::abs(Emag0), 1e-300);
                std::cout << "[HallOnlyEnergy] dt=" << dt
                          << " Emag0=" << Emag0
                          << " Emag1=" << Emag1
                          << " dE=" << dE
                          << " rel=" << rel
                          << std::endl
                          << std::endl
                          << std::endl;
            }
        }
#endif

        run_data_->dt = dt;
        run_data_->time += dt;
        run_data_->step += 1;

        mercury_bound_.Sync("Ucell");
        mercury_bound_.Sync("Bface");

        calc_Bcell();
        calc_Jcell_from_Bcell_metric_();
        calc_PV();
        calc_Uplus();

        control_.UpdateSwitches(*run_data_);
        UpdateControlAndOutput();
    }
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
