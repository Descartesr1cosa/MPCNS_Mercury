#include "5_io/RuntimeMonitor.h"
#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "5_io/RunData.h"
#include <time.h>
#include <cstdio>

static inline double cpu_time_seconds_()
{
#if defined(CLOCK_PROCESS_CPUTIME_ID)
    timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return double(ts.tv_sec) + 1e-9 * double(ts.tv_nsec);
#else
    return double(std::clock()) / double(CLOCKS_PER_SEC);
#endif
}

void RuntimeMonitor::Begin(RunData &run, Param *par, int64_t global_cells)
{
    par_ = par;
    global_cells_ = global_cells;
    cpu_last_s_ = cpu_time_seconds_();
    step_last_ = run.step;
    inited_ = true;

    // 只存通用 perf 信息（不带物理含义）
    run.SetMonitor("Perf/Cells", (int)global_cells_);
    run.SetMonitor("Perf/CPUInterval_s", 0.0);
    run.SetMonitor("Perf/s_per_MCells_step", 0.0);
}

void RuntimeMonitor::UpdateOnOutres(RunData &run)
{
    if (!inited_)
        return;

    const double cpu_now_s = cpu_time_seconds_();
    double cpu_interval_l = cpu_now_s - cpu_last_s_;
    if (cpu_interval_l < 0.0)
        cpu_interval_l = 0.0; // 防御
    cpu_last_s_ = cpu_now_s;

    const int32_t step_now = run.step;
    const int32_t step_delta = (step_now > step_last_) ? (step_now - step_last_) : 1;
    step_last_ = step_now;

    double cpu_interval_g = cpu_interval_l;
    PARALLEL::mpi_max(&cpu_interval_l, &cpu_interval_g, 1);

    const double denom = std::max(1e-30, double(step_delta) * (double(global_cells_) * 1e-6));
    const double s_per_mcell_step = cpu_interval_g / denom;

    run.SetMonitor("Perf/CPUInterval_s", cpu_interval_g);
    run.SetMonitor("Perf/s_per_MCells_step", s_per_mcell_step);
}

void RuntimeMonitor::PrintLineMinimal(const RunData &run) const
{
    if (!par_ || par_->GetInt("myid") != 0)
        return;

    const double cpu_dt = run.GetMonitorD("Perf/CPUInterval_s", 0.0);
    const double s_per = run.GetMonitorD("Perf/s_per_MCells_step", 0.0);
    const int32_t cells = run.GetMonitorI("Perf/Cells", 0);

    const double mdofs = double(cells) * 1e-6; // 暂用 Cells 作为 DOFs

    // ---- line 1 ----
    // Nstep: 6 columns, t/dt: width 10 with 5 decimals (adjust if you want wider)
    std::printf("[  Diag  ] Nstep   =%9d            t  =%10.5f        dt = %10.3e\n",
                run.step, run.time, run.dt);

    // ---- line 2 ----
    // indent: 11 spaces aligns under after "[Diag    ] "
    // Keys aligned by giving each key field a fixed width of 8 chars (right-aligned) before '=':
    // "CPU_Time", "DOFs", "speed" => '=' will align.
    std::printf("           %8s=%10.4f s    %8s=%8.3fM    %8s=%10.3f   s/MDOF/step\n",
                "CPU_Time", cpu_dt,
                "DOFs", mdofs,
                "speed", s_per);

    std::fflush(stdout);
}