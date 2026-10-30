#pragma once
// Optional trial-run profiler used to calibrate the Preprocess fluid/solid
// cost weights from measurements instead of empirical values.
//
// Disabled by default: every hook then reduces to one test of a global bool and
// no numerical state is read or written.  When enabled (setup key
// `profile_run`), each rank records, per solver phase, the wall time and the
// time spent blocked in MPI synchronisation (compute = wall - blocked), then
// rank 0 gathers the data, fits
//     compute_r = a_fluid * alloc_fluid_cells_r + a_solid * alloc_solid_cells_r + c * blocks_r
// where alloc cells include the ghost layers (per-block cost is dominated by the
// ghost-inclusive volume).
// and writes DATA/partition_profile.json.
#include <string>
#include <vector>

namespace PROFILE
{
extern bool g_enabled;
extern double g_blocked_seconds; // accumulated by the PMPI wrappers

struct BlockInfo
{
    std::string name;
    bool fluid = false;
    long long cells = 0;       // inner cells
    long long alloc_cells = 0; // cells including ghost layers
};

class PartitionProfiler
{
public:
    static constexpr int kMaxPhases = 16;

    void configure(bool enabled, int warmup_steps, int measured_steps);
    bool enabled() const { return g_enabled; }
    void set_blocks(std::vector<BlockInfo> blocks) { blocks_ = std::move(blocks); }
    void set_phase_name(int id, const char *name) { names_[id] = name; }
    void set_ghost_layers(int g) { ghost_layers_ = g; }

    // Called once per step after the phases; returns true when the measured
    // window is complete and the run should stop.
    bool end_step();

    void add_phase(int id, double wall, double blocked)
    {
        if (!measuring_)
            return;
        wall_[id] += wall;
        blocked_[id] += blocked;
    }

    // Collective: gather, fit and write the JSON file (rank 0).
    void write_json(const std::string &path, int nranks, int rank) const;

private:
    std::vector<BlockInfo> blocks_;
    std::string names_[kMaxPhases];
    double wall_[kMaxPhases] = {};
    double blocked_[kMaxPhases] = {};
    int warmup_ = 20, measured_ = 100, step_ = 0, ghost_layers_ = 0;
    bool measuring_ = false;
};

// RAII phase timer; a single branch on g_enabled when profiling is off.
class PhaseScope
{
public:
    PhaseScope(PartitionProfiler &p, int id) : p_(p), id_(id), on_(g_enabled)
    {
        if (on_)
            start_();
    }
    ~PhaseScope()
    {
        if (on_)
            stop_();
    }

private:
    void start_();
    void stop_();
    PartitionProfiler &p_;
    int id_;
    bool on_;
    double t0_ = 0.0, b0_ = 0.0;
};
} // namespace PROFILE
