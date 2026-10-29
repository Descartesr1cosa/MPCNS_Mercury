#pragma once
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "5_boundary/Boundary.h"
#include "6_io/IOModule.h"
#include "7_metric/Metric.h"
#include "8_dec/DecOps.h"
namespace MPCNS
{
enum class PhysicalCornerMode { MercuryCompatible, StageLocal, SolverManaged };
struct RuntimeIOConfig
{
    bool enabled = true;
    int residual_components = 0;
    std::vector<std::string> output_fields, restart_fields;
};
struct SimulationConfig
{
    int geometry_ghost = -1, dimension = -1;
    std::vector<FieldDescriptor> fields;
    std::vector<Field::PairKey> coupling_pairs;
    bool enable_metric_fields = false;
    PhysicalCornerMode physical_corner_mode = PhysicalCornerMode::MercuryCompatible;
    RuntimeIOConfig io;
    std::function<void(Grid &, TOPO::Topology &)> configure_topology;
};
struct SyncOptions
{
    bool physical = true, halo = true, coupling = true, owner_alias = true;
    std::optional<HaloLevel> through;
};
class MetricAccess
{
public:
    bool enabled() const noexcept { return enabled_; }
    void recompute();
    METRIC::MetricDiagnostics diagnostics();
private:
    friend class SimulationContext;
    MetricAccess(Field *f, Grid *g, bool e) : field_(f), grid_(g), enabled_(e) {}
    Field *field_; Grid *grid_; bool enabled_;
};
class DecAccess
{
public:
    void d0(const std::string &, const DEC::EdgeFormNames &, int = 0, int = 0);
    void d1(const DEC::EdgeFormNames &, const DEC::FaceFormNames &, int = 0, int = 0);
    void d2(const DEC::FaceFormNames &, const std::string &, int = 0, int = 0);
    void faraday(const DEC::EdgeFormNames &, const DEC::FaceFormNames &, double, int = 0, int = 0);
private:
    friend class SimulationContext;
    explicit DecAccess(Field *f) : field_(f) {}
    Field *field_;
};
class SimulationContext
{
public:
    SimulationContext(Param &, int, SimulationConfig);
    ~SimulationContext();
    SimulationContext(const SimulationContext &) = delete;
    SimulationContext &operator=(const SimulationContext &) = delete;
    void initialize();
    bool initialized() const noexcept { return initialized_; }
    std::uint64_t communication_generation() const noexcept { return communication_generation_; }
    void register_physical_handler(const std::string &, const std::string &, BOUND::PhysicalHandler);
    void register_coupling_handler(const std::string &, const std::string &, StaggerLocation,
                                   const std::string &, const std::string &, BOUND::CouplingHandler);
    void validate_physical_handlers() const;
    void sync(const std::string &);
    void sync(const std::string &, const SyncOptions &);
    void sync_group(const std::string &);
    void sync_group(const std::string &, const SyncOptions &);
    // Legacy physics adapters may still require their case parameters.  The
    // context only borrows this object; it neither owns nor mutates its
    // lifetime.  New infrastructure should prefer typed config instead.
    Param &param() noexcept { return param_; }
    const Param &param() const noexcept { return param_; }
    Grid &grid(); const Grid &grid() const;
    TOPO::Topology &topology(); const TOPO::Topology &topology() const;
    Field &field(); const Field &field() const;
    Halo &halo(); const Halo &halo() const;
    BoundaryCore &boundary(); const BoundaryCore &boundary() const;
    IOModule &io(); const IOModule &io() const;
    RunData &run_data(); const RunData &run_data() const;
    MetricAccess metric(); DecAccess dec();
    const std::vector<std::string> &output_fields() const noexcept { return config_.io.output_fields; }
    const std::vector<std::string> &restart_fields() const noexcept { return config_.io.restart_fields; }
private:
    struct SyncTarget { std::string group; std::vector<std::string> fields; bool by_field = false; };
    void require_initialized_(const char *) const;
    void validate_config_() const;
    void validate_collective_configuration_() const;
    void prepare_oriented_exchanges_();
    SyncTarget resolve_target_(const std::string &, bool) const;
    HaloLevel target_level_(const SyncTarget &, const SyncOptions &) const;
    std::vector<std::string> physical_fields_(const SyncTarget &) const;
    std::vector<int32_t> coupling_ids_(const SyncTarget &, const Field::PairKey &) const;
    void apply_physical_(const std::vector<std::string> &, HaloLevel);
    void sync_halo_(const SyncTarget &, HaloLevel);
    void sync_coupling_(const SyncTarget &, HaloLevel);
    void sync_owner_(const SyncTarget &);
    void sync_resolved_(const SyncTarget &, const SyncOptions &);
    Param &param_; int rank_; SimulationConfig config_; Grid grid_;
    std::unique_ptr<TOPO::Topology> topology_;
    std::unique_ptr<Field> field_; std::unique_ptr<Halo> halo_;
    std::unique_ptr<BoundaryCore> boundary_; std::unique_ptr<IOModule> io_;
    bool initialized_ = false;
    std::uint64_t communication_generation_ = 0;
};
} // namespace MPCNS
