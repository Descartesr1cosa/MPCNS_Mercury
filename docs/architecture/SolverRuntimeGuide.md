# Solver Runtime API

`MPCNS::SimulationContext` is the thin assembly layer for a new solver. It owns
`Grid`, `Topology`, `Field`, `Halo`, `BoundaryCore`, and (when enabled)
`IOModule`. It does not initialize or finalize MPI, and it does not own `Param`.
The caller must keep `Param` alive until the context is destroyed.

Include the public header:

```cpp
#include "9_runtime/SimulationContext.h"
```

## Declare fields

Fields use the existing `FieldDescriptor` and `FieldSyncContract`. A synchronized
field must have a non-empty group. An orientation-aware 1-form or 2-form group
must contain the complete xi/eta/zeta triplet with matching ghost width, level,
and value kind.

```cpp
FieldDescriptor state;
state.name = "Ucell";
state.location = StaggerLocation::Cell;
state.value_kind = FieldValueKind::ConservedVector;
state.ncomp = 5;
state.nghost = 2;
state.sync.group = "Ucell";
state.sync.do_halo = true;
state.sync.do_physical = true;
state.sync.do_coupling = true;
state.sync.halo_level = HaloLevel::Corner3D;

FieldDescriptor bxi;
bxi.name = "B_xi";
bxi.location = StaggerLocation::FaceXi;
bxi.value_kind = FieldValueKind::FaceContravariant2Form;
bxi.ncomp = 1;
bxi.nghost = 2;
bxi.sync.group = "Bface";
bxi.sync.do_halo = true;
bxi.sync.orientation_aware = true;
bxi.sync.owner_sync = OwnerSyncPolicy::FaceOwner;
bxi.sync.halo_level = HaloLevel::Corner3D;

FieldDescriptor beta = bxi;
beta.name = "B_eta";
beta.location = StaggerLocation::FaceEt;
FieldDescriptor bzeta = bxi;
bzeta.name = "B_zeta";
bzeta.location = StaggerLocation::FaceZe;
```

`Scalar`, `CartesianVector`, `ConservedVector`, `EdgeCovariant1Form`, and
`FaceContravariant2Form` retain their existing Halo and orientation semantics.
Owner policies are independent of `do_halo`; owner reconciliation is scheduled
after all selected Halo and Coupling stages.

## Initialize the context

MPI remains process-level infrastructure. Initialize it once before creating a
context and finalize it after the context has been destroyed.

```cpp
PARALLEL::mpi_initial(argc, argv);
{
    int rank = 0;
    PARALLEL::mpi_rank(&rank);

    Param param;
    param.ReadParam(rank);

    MPCNS::SimulationConfig config;
    config.fields = {state, bxi, beta, bzeta};
    config.coupling_pairs = {{"Solid", "Fluid"},
                             {"Fluid", "Solid"}};
    config.enable_metric_fields = true;
    config.io.residual_components = 5;
    config.io.output_fields = {"Ucell", "B_xi", "B_eta", "B_zeta"};
    config.io.restart_fields = config.io.output_fields;

    MPCNS::SimulationContext sim(param, rank, std::move(config));
    sim.initialize(); // collective

    // Register handlers and run the solver while sim is alive.
}
PARALLEL::mpi_finalize();
```

`initialize()` performs, once and in dependency order:

1. Grid preprocessing and generic topology construction.
2. Optional solver patch annotation through `config.configure_topology`.
3. Geometry, optional Metric fields, and declared Field allocation.
4. Declared Coupling channel and buffer construction.
5. Halo, Owner/Alias, Coupling, and oriented prepared-plan construction.
6. Shared workspace allocation and prepared oriented-path activation.
7. Physical boundary pattern construction and optional IO setup.

The configuration is fingerprinted collectively before construction. All ranks,
including ranks without a local block or local Coupling patch, must call
`initialize()` and use the same field and Coupling declarations. Field
registration is intentionally closed after initialization; create a new context
for a different registry instead of rebuilding metadata inside a timestep.

`configure_topology` may assign solver domain names or Coupling flags. It must
not change entity equivalence, connectivity geometry, or patch transforms.

## Register Physical and Coupling handlers

Handlers are registered after `initialize()`, when Boundary patterns exist:

```cpp
sim.register_physical_handler(
    "Ucell", "Farfield",
    [](FieldBlock &u, Field *fields,
       const BOUND::PhysicalRegion &region, int nghost) {
        BoundaryCore::DefaultPhysicalCopy(u, fields, region, nghost);
    });

sim.register_coupling_handler(
    "Solid", "Fluid", StaggerLocation::Cell,
    "Ucell", "Ucell",
    [](FieldBlock &destination, Field *fields,
       CouplingBufferBlock &buffer,
       const std::string &source_domain,
       const std::string &destination_domain,
       const std::string &channel) {
        BoundaryCore::DefaultCouplingCopy(
            destination, fields, buffer,
            source_domain, destination_domain, channel);
    });

sim.validate_physical_handlers();
```

The solver supplies physical behavior. The Runtime only stores callbacks and
schedules them. `validate_physical_handlers()` is explicit so a solver can
finish installing all handlers before requesting strict validation.

## Synchronize fields and groups

Both field names and declared group names are accepted:

```cpp
sim.sync("Ucell");
sim.sync_group("Bface");
// A member of an oriented triplet resolves to its complete group.
sim.sync("B_xi");
```

Each stage preserves Mercury ordering:

```text
Physical Boundary -> Halo -> Coupling transfer/application
```

After the requested Face/Edge/Vertex stages finish, Owner/Alias reconciliation
runs once. It therefore cannot overwrite a pending Coupling update.

Individual phases and the maximum stage can be selected without exposing MPI:

```cpp
MPCNS::SyncOptions options;
options.coupling = false;
options.owner_alias = false;
options.through = HaloLevel::Corner2D;
sim.sync("Ucell", options);
```

`PhysicalCornerMode::MercuryCompatible` is the default. `StageLocal` applies the
default physical corner fill immediately before each matching stage;
`SolverManaged` leaves corner treatment entirely to solver logic. Physical,
Halo, Coupling, and Owner phases may also be disabled independently in
`SyncOptions`.

## Access Field, Metric, DEC, IO, and RunData

```cpp
FieldBlock &u0 = sim.field().field("Ucell", 0);
const TOPO::Topology &topology = sim.topology();

if (sim.metric().enabled()) {
    sim.metric().recompute();
    METRIC::MetricDiagnostics diagnostics = sim.metric().diagnostics();
}

sim.dec().d1(DEC::EdgeFormNames{"E_xi", "E_eta", "E_zeta"},
             DEC::FaceFormNames{"B_xi", "B_eta", "B_zeta"});

sim.run_data().step += 1;
sim.io().WriteRestartBinFile(sim.run_data().step, sim.run_data().time);
```

`metric()` is available only when `enable_metric_fields` was set. `io()` and
`run_data()` require `config.io.enabled`; IO is enabled by default. The existing
`Field`, `Halo`, `BoundaryCore`, `IOModule`, Metric, and DEC APIs remain usable
through these accessors and are not reimplemented by the Runtime.
