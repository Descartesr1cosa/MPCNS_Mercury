# MPCNS Mercury Refactor Plan

Status: Accepted design baseline; Units 1-20 verified by the 3-rank Z0 suite. Unit 20 isolated the cell covariant-basis and regularized 3x3 pseudo-inverse kernels, and moved singular-edge legacy Hodge publication behind compat. All gates passed in all ten halo modes.
Owner: Ke Xu
Scope: bottom-layer architecture; preserve the validated Mercury physics

## 1. Objective and constraints

This document is the architecture contract for the structured multi-block, DEC-based Mercury multi-fluid MHD solver. The goal is to make geometry, discrete operators, communication, and physics independently evolvable without changing the validated production result accidentally.

- The present calculation is the scientific baseline. A refactor does not change stencils, boundary semantics, synchronization order, or results without an approved numerical or physical change.
- LumpedHodge remains the production default. ConsistentHodge is experimental until separate verification and Mercury regression pass.
- All local builds use -j4 or less and all local MPI runs use -np 4 or less.
- Geometry is immutable after startup. Simulation state and workspace are mutable.
- Every phase must build, test, and be reversible independently.

## 2. Naming standard

All new public code uses namespace mpcns.

~~~text
mpcns::core       basic types, ranges, errors, array views
mpcns::mesh       blocks, topology, and partition metadata
mpcns::geometry   coordinates, primal and dual geometry
mpcns::state      cochains, state layout, and storage
mpcns::discrete   incidence, Hodge, CT, and reconstruction
mpcns::parallel   exchange and owner/alias plans
mpcns::boundary   geometric, topological, and physical constraints
mpcns::linalg     distributed DoF maps and PETSc adapters
mpcns::physics    reusable physical closures
mpcns::solver     time integration and operator composition
mpcns::io         restart and visualization adapters
mpcns::apps::mercury  Mercury-specific configuration and wiring
~~~

| Kind | Rule | Example |
|---|---|---|
| Type | UpperCamelCase | GeometrySnapshot |
| Function or variable | lower_snake_case | build_exchange_plan() |
| Private member | lower_snake_case_ | geometry_ |
| Compile-time constant | kUpperCamelCase | kMaxGhostLayers |
| Boolean | is_, has_, needs_ prefix | is_owner |
| Source file | lower_snake_case.hpp/.cpp | lumped_hodge.hpp |
| Axis | full spelling | Axis::Xi, Axis::Eta, Axis::Zeta |

New public APIs do not expose Et, Ze, alpha, beta, or generic data. Names state mathematical meaning and representation.

| Legacy name | Canonical name |
|---|---|
| B_xi/B_eta/B_zeta | induced_magnetic_flux, primal 2-form |
| Badd_* | background_magnetic_flux, primal 2-form |
| E_* / Ehall_* / Eres_* | electric_emf / hall_electric_emf / resistive_electric_emf |
| J_* | electric_current, primal 1-form |
| B_cell / Bind_cell / J_cell | Cartesian derived field names |
| JDxi/JDet/JDze | face_area_vector |
| dr_*, dl_* | edge_vector, edge_length |
| dlstar_*, Sstar_*, Astar_* | dual geometry at an edge or face |

Legacy strings may remain in case files, restarts, and IO name maps. They do not control new numerical operators.

## 3. Source tree and dependency rule

~~~text
include/mpcns/
  core/  mesh/  geometry/  state/  discrete/  compat/  parallel/
  boundary/  linalg/  physics/  solver/  io/  apps/mercury/
src/ mirrors include/mpcns/
~~~

`compat/` is a temporary, explicit boundary for legacy Grid/Field layouts and
catalog names.  It may depend on legacy classes; the target geometry module
may not.  It is removed only after all applications are constructed through
the new geometry and state APIs.

Each area is a separate CMake target: mpcns_core, mpcns_mesh, mpcns_geometry, mpcns_state, mpcns_discrete, mpcns_parallel, mpcns_boundary, mpcns_linalg, mpcns_physics, mpcns_solver, mpcns_io, and mpcns_mercury_app.

~~~text
core -> mesh, state
geometry -> mesh
discrete -> mesh + geometry + state
parallel -> mesh + state
boundary -> mesh + geometry + state
linalg -> discrete + parallel
physics -> discrete + state
solver -> physics + boundary + linalg
apps/mercury -> solver + io
~~~

Geometry includes neither MPI, PETSc, Field, Halo, nor Mercury headers. Discrete includes neither MPI, PETSc, nor application headers. Physics makes no direct MPI call. Parallel implements no physical stencil or Hodge formula.


## 4. Mesh topology and MPI partition

Topology and partitioning are independent objects.

~~~cpp
namespace mpcns::mesh {
class MeshTopology;
class MeshPartition;
}
~~~

MeshTopology has no rank-dependent meaning. It owns block extents and transforms, global entity equivalence classes, orientation signs, d0/d1/d2 incidence, boundary/interface classification, and extraordinary-edge valence.

MeshPartition owns local and ghost entity maps, rank ownership, the neighbor graph, and distributed DoF numbering. Therefore the one-, two-, and four-rank forms of a case change MeshPartition but not MeshTopology.

Topology owns neither coordinates, geometric measures, Hodge values, field values, MPI requests, nor PETSc objects.

## 5. Unified geometry design

The legacy grid metric and newer DEC metric become one immutable data flow:

~~~text
NodeCoordinates
  -> PrimalGeometry
       cell measure and Jacobian
       covariant/contravariant bases
       face area vector and face area
       edge vector and edge length
       GCL face metrics and tangent bases
  -> DualGeometry
       dual edge length at faces
       dual face geometry at edges
       extraordinary-edge and axis closure
  -> Hodge caches and ReconstructionCache
~~~

~~~cpp
namespace mpcns::geometry {
class NodeCoordinates;
class PrimalGeometry;
class DualGeometry;
class ReconstructionCache;
class GeometrySnapshot {
public:
    const NodeCoordinates& coordinates() const;
    const PrimalGeometry& primal() const;
    const DualGeometry& dual() const;
    const ReconstructionCache& reconstruction() const;
};
}
~~~

Coordinate ghost exchange is a one-time parallel preparation step. GeometrySnapshot is then built from a closed local mesh view and remains read-only for the run; geometry code itself performs no MPI.

`GeometryBuilder` is the typed structural boundary between a geometry producer and geometry consumers. It accepts `GeometryBuildInput` (ordered `GeometryBlock` views) and returns `GeometryBuildResult` (a complete `GeometrySnapshot` plus its block count). During migration, `LegacyGeometryBuildPipeline` preserves the five legacy publisher calls and supplies those views to the builder; the geometry target still sees neither `Grid`, `Field`, MPI, nor solver types.

The dual polygon around an extraordinary edge is now closed by the pure `close_singular_edge_geometry` operation. `SingularEdgeRegistry` remains responsible only for quotient topology, ownership, MPI reductions, legacy field publication, and collecting the typed input. This keeps the physical/topological seam explicit without changing the current closure formula.

The reconstruction cache now delegates its cell covariant-basis average and regularized 3x3 pseudo-inverse to pure geometry kernels. The compatibility publisher retains only Grid/Field reads, the historical edge-local stencil selection, and field writes. The pseudo-inverse preserves the existing Gram-product order, `1e-14 * (trace + 1e-300)` regularization, eight retry attempts, `1e-300` determinant threshold, and zero fallback exactly.

PrimalGeometry is the only source for finite-volume and DEC primal metrics. A legacy JDxi-like face vector and a corresponding GCL face metric must use the same coordinate stencil or be checked by a geometry consistency diagnostic.

DualGeometry owns regular and variable-valence closure uniformly. A pole or extraordinary edge is a geometry/topology fact, not a branch hidden in Mercury physics.


## 6. State and cochain data model

Mathematical identity, storage layout, and communication requirements are separate concepts.

~~~cpp
namespace mpcns::state {
enum class ComplexKind { Primal, Dual };
enum class FormDegree { Zero, One, Two, Three };
template <ComplexKind, FormDegree> class CochainView;
using PrimalEdgeForm = CochainView<ComplexKind::Primal, FormDegree::One>;
using PrimalFaceForm = CochainView<ComplexKind::Primal, FormDegree::Two>;
class QuantityDescriptor;
class StateLayout;
class SimulationState;
}
~~~

QuantityDescriptor records mathematical identity, units, components, and external name. StateLayout records owned and ghost memory. Exchange requirements belong to an operation, not a field descriptor. Hot kernels use typed handles/views rather than string field names.

~~~cpp
struct MercuryState {
    FluidConservativeState hydrogen;
    FluidConservativeState sodium;
    state::PrimalFaceForm induced_magnetic_flux;
    state::PrimalFaceForm background_magnetic_flux;
    state::PrimalEdgeForm electric_emf;
    state::PrimalEdgeForm electric_current;
};
~~~

## 7. Discrete operators and Hodge

DEC owns the only incidence implementation. CTOperators is first converted to a compatibility wrapper and later removed.

~~~cpp
class IncidenceOperator {
public:
    void apply_d0(const NodeScalarField&, PrimalEdgeForm&) const;
    void apply_d1(const PrimalEdgeForm&, PrimalFaceForm&) const;
    void apply_d2(const PrimalFaceForm&, CellVolumeForm&) const;
    void apply_d1_transpose(const DualEdgeForm&, DualFaceForm&) const;
};
class HodgeOperator {
public:
    virtual void apply_star_2(const PrimalFaceForm&, DualEdgeForm&) const = 0;
    virtual void solve_star_1(const DualFaceForm&, PrimalEdgeForm&) const = 0;
};
~~~

LumpedHodge encapsulates the current diagonal scales formerly exposed as Face_beta and Edge_alpha. ConsistentHodge later implements the same interface with local Whitney/FEEC mass matrices and, if needed, a distributed solve.

~~~text
Hodge2ToDual1: primal 2-form -> dual 1-form
InverseHodge1: solve primal 1-form from a dual 2-form right-hand side
J = M1^{-1} D1^T M2 B
~~~

Reconstruction is separate from Hodge: face form to Cartesian cell vector and edge form to Cartesian cell vector. Hall, CT, output, and Mercury boundary code do not know which Hodge model is active.

The existing `METRIC::reconstruct_*` and `METRIC::diagnose_metric_fields` names
remain as backwards-compatible facades. Their legacy formulas and field
catalog reads live in `mpcns::compat::legacy_field_reconstruction` and
`mpcns::compat::legacy_metric_diagnostics`, respectively. This deliberately
separates result consumers from metric publication without changing either
public API or numerical arithmetic.

## 8. Parallel, boundaries, and implicit solve

ExchangePlan compiles an operator request into halo exchange, reduce-to-owner, and broadcast-from-owner operations. The executor begins exchange, computes the interior, finishes exchange, computes interfaces, and reconciles owners only if needed.

~~~text
GeometryClosure        axis, pole, dual closure, degenerate coordinates
TopologicalConstraint  panel equivalence, orientation, canonical values
PhysicalBoundary       far field, conducting wall, solid/fluid interface
~~~

PETSc is behind an adapter. ResistiveDiffusionOperator combines MeshPartition, incidence, Hodge, eta, and constraints. PetscLinearAdapter owns distributed Vec/Mat/KSP, options, and diagnostics. MercurySolver only decides solve timing and commits the magnetic-flux update.

An implicit operator application is observationally pure: input vector to output vector. Its workspace is not observable mutable MercurySolver state. This permits an assembled reference operator, useful preconditioners, and matrix-free execution to coexist.


## 9. Current-to-target mapping

| Current area | Target area |
|---|---|
| 0_basic | core |
| 1_grid | mesh/structured and geometry/primal |
| 2_topology | mesh/topology |
| 3_field | state |
| 4_halo | parallel |
| 5_boundary | boundary |
| 6_io | io |
| 7_metric/Metric.cpp | geometry, discrete/hodge, discrete/reconstruction |
| 7_metric/SingularEdgeRegistry | mesh/topology, geometry/dual, parallel |
| 8_dec | discrete/incidence |
| Z4_Mercury | apps/mercury |

## 10. Migration phases and gates

### Phase 0 — Baseline and guardrails

Freeze representative MESSENGER-validated cases, output epochs, diagnostics, and one-to-four-rank comparisons. Record compiler, MPI, PETSc, case revision, and acceptance tolerances. No numerical change begins before this regression command exists.

### Phase 1 — Typed compatibility facade

Introduce namespace, CMake target boundaries, typed cochain views, and GeometrySnapshot interfaces. Existing Grid, Field, and metric arrays remain the backend. Gate: no Mercury kernel change, common storage is used, and baseline output is unchanged.

### Phase 2 — Single incidence implementation

Make IncidenceOperator the only d0/d1/d2 implementation and provide d1 transpose. Convert CT to a compatibility caller. Gate: legacy CT and DEC d1 match exactly on random cochains; d1*d0 = 0 and d2*d1 = 0; orientation tests pass.

Current incremental scope: `d1` has been centralized behind the transitional
`mpcns::discrete::apply_d1` facade, with the former CT stencil retained only as
a Z0 regression oracle. `d0`, `d2`, and the transpose remain unchanged until
this unit has been reviewed and run in the Ubuntu/cluster build environment.

The next completed seam is documented in
`docs/architecture/GeometryInventory.md`: cochain storage now belongs to
`mpcns::state`, and all Mercury geometry/Hodge catalog-ID binding goes through
`mpcns::compat::LegacyGeometryFieldIds`.  This adds no geometry calculation or
solver behavior; it isolates the legacy backend before `GeometrySnapshot` is
introduced.

`mpcns::geometry::GeometrySnapshot` is now present as a geometry-only,
read-only API.  `mpcns::compat::borrow_legacy_geometry_snapshot` supplies its
temporary Grid/Field-backed implementation.  It is validated by exact Z0
readback of coordinates and every current primal/dual metric category; no
Mercury production path consumes it yet.

The reconstruction geometry cache now follows the same path: covariant bases
and the cell/edge pseudo-inverses are semantic members of GeometrySnapshot,
and Mercury's SolverFields binds their legacy catalog IDs only through
LegacyGeometryFieldIds.  This removes the remaining geometry-name lookups
from the Mercury solver facade without changing a metric calculation.

Legacy metric catalog spellings are now declared once in
`mpcns::compat::legacy_geometry_field_names`; both legacy metric-field
registration and the compatibility ID adapter consume this directory.  This
is a naming-boundary change only and requires the usual exact Z0 gate before
subsequent geometry-builder extraction.

The next seam is the geometry field-publishing contract.  The 43 legacy
geometry field registrations (catalog name, stagger location, component
count, ghost-width rule, and historical order) now live in
mpcns::compat::legacy_geometry_field_registration.  Metric.cpp is only a
backward-compatible caller.  The Z0 suite checks that the resulting catalog
matches that contract exactly; metric formulas and storage are unchanged.

The first pure metric kernel is now geometry_primitives::quadrilateral_area_vector.
It computes the oriented area vector of a non-planar quadrilateral using the
same two-triangle arithmetic as the legacy code.  Metric.cpp calls it only for
JDxi, JDet, and JDze.  Z0 independently recomputes every local JD value with
the former stencil and requires an exact zero difference before any subsequent
primal-geometry extraction proceeds.

The next pure kernel is geometry_primitives::cell_volume_from_oriented_faces.
It preserves the historical ordered six-face divergence-theorem accumulation:
xi-minus, xi-plus, eta-minus, eta-plus, zeta-minus, zeta-plus.  Metric.cpp
still owns the legacy face centroids and Field access; only the pure quadrature
is moved.  Z0 independently evaluates every local Jac entry using the prior
stencil and requires exact equality before the next primitive is extracted.

geometry_primitives::euclidean_norm is now used by the three primal-face
metric paths, the three primal-edge length paths, and the three dual-face
magnitude paths. Its three-rank Z0 gate has passed: every local edge length
and dual-face magnitude matches the former norm stencil exactly when combined
with the independent Sstar gate.

The axis-degeneracy treatment for inverse star-1 is represented by the pure
mpcns::discrete::AxisSafeInverseStar1Policy. Its default parameters preserve
the existing dual-collapse ratio (1e-3) and primal-edge absolute tolerance
(1e-12), while its DualFaceQuadrilateral input makes the four legacy
dual-center points explicit. The Z0 gate passed for constructed regular and
collapsed quadrilaterals and every field-derived inverse star-1 coefficient.

The existing face/edge metric publisher now owns the legacy catalog writes for
primal and dual edge/face measures plus inverse star-1. Unit 16 further splits
reconstruction-cache publication and the diagonal Hodge closure into two
compatibility entry points. Their call order and field backend are unchanged.
The combined 3-rank Z0 validation passed with exact field equality in all ten halo modes.

Unit 17 then moved primal face/cell and edge-vector publishing behind the same
compatibility boundary. The established JD/Jac gates and a new full-range dr_*
stencil gate also passed exactly in all ten modes.


The current unintegrated numerical foundation also exposes `d1^T` and a
`HodgeOperator` interface with `LumpedHodge` as its production implementation.
The existing CT adjoint curl is deliberately still active until a workspace and
an exact old-versus-new Hodge/adjoint regression are available.

That regression is now prepared in the Z0 DEC chain suite.  It compares the
current fused CT current stencil with `star_2 -> d1_transpose -> inverse_star_1`
on deterministic face-form data, including local ghost values.  The test must
pass exactly before the Mercury production call sites change.

### Phase 3 — GeometrySnapshot and LumpedHodge

Split the current metric implementation into primal geometry, dual geometry, lumped Hodge, reconstruction cache, and diagnostics. Gate: old/new geometry measures and Hodge action agree; extraordinary edges and poles remain stable; Mercury regression remains within tolerance.

### Phase 4 — Operator-facing Mercury physics

Replace direct field-id, Hodge, and stencil usage with electromagnetic and reconstruction operators plus typed state views. Gate: Mercury physics has no incidence stencil or direct Hodge-scale access and preserves baseline output.

### Phase 5 — ExchangePlan and constraints

Compile communication and ownership from MeshPartition and operator requests. Move physical boundary behavior out of generic halo execution. Gate: face, edge, vertex, orientation, and owner NaN-fill tests pass; one-to-four-rank results match; no unintended synchronization is added.

### Phase 6 — Implicit operator extraction

Move resistive implicit logic into ResistiveDiffusionOperator and PetscLinearAdapter. Add an assembled reference before optimizing matrix-free execution. Gate: linearity and manufactured-diffusion tests, PETSc residual checks, and baseline agreement pass.

### Phase 7 — Experimental ConsistentHodge

Add non-diagonal Hodge behind the established interface. It remains disabled by default until convergence, stability, and Mercury comparisons are complete.

## 11. First implementation unit

1. Create the mpcns::discrete compatibility header and source skeleton.
2. Expose typed 1-form and 2-form views over current IdTriplet storage.
3. Route the CT d1 operation through the single DEC implementation.
4. Add a legacy-versus-new exact equivalence test.
5. Build and test using no more than four cores/ranks.

No Hodge formula, MPI schedule, boundary treatment, solver parameter, or Mercury case setup changes in this first unit.

## 12. Unit 21: singular-edge data-source seam

The remaining singular-edge Grid/Field reads are now behind `mpcns::compat::LegacySingularEdgeDataSource` and `collect_legacy_singular_edge_measures`. The registry is reduced to topology-to-quotient mapping, its unchanged MPI reductions, valence filtering, pure dual-polygon closure, normalized contribution weights, and catalog publication. The shared alias and incident-entity records live in storage-independent `mpcns::mesh`.

This is intentionally a boundary-only migration: no Hodge formula, physical-edge classification rule, collective ordering, mesh partition behavior, or Mercury solve path changes. Both Z0 and MercuryZ4 compile with at most four jobs; the full three-rank, ten-mode Z0 gate passed. The next safe step is to make the collector consume a typed geometry/partition view rather than Field/Grid, only after a Mercury case-level comparison is available.

## 13. Unit 23: complete incidence and electromagnetic composition

`mpcns::discrete` is now the sole implementation of the three primal incidence maps: `apply_d0`, `apply_d1`, and `apply_d2`. The old `DEC` entry points remain compatibility facades only, including the former Faraday update. Thus a future physics module can depend on typed 0/1/2/3-form views without importing a Field catalog or duplicating an orientation stencil.

`mpcns::operators::ElectromagneticOperators` provides the allocation-free current composition `inverse_star_1(d1_transpose(star_2(B)))`. Scratch views are caller-owned and the caller remains responsible for halo/owner synchronization before invoking an operator. The operator owns no Grid, Field, topology construction, MPI, boundary condition, or time-step decision.

This is the intended MHD entry point: ideal, Hall, ambipolar, resistive, and implicit modules should first construct EMF/current cochain views, use this composition, and keep their physical closures outside the discrete layer.

## 14. Unit 24: Mercury current adoption

Mercury now owns a private face-located `Jstar_*` workspace with no output, coupling, physical-boundary, halo, or owner-sync contract. `ComputeRegularEdgeCurrent_` binds the unchanged legacy diagonal Hodge coefficients and invokes `ElectromagneticOperators` for both the production `B -> J` path and the implicit-resistive `dB -> dJ` path.

Alias reduction, singular-edge replacement, and the existing Jedge/dJ boundary synchronization remain outside the local operator in exactly their prior order. Consequently this migration removes the duplicated regular-edge Hodge/curl loops without changing the quotient-edge or boundary algorithms. MercuryZ4 compilation and the generated three-rank cube gate passed; a protected production Mercury case remains the next numerical acceptance gate.

## 15. Unit 25: topology-only Faraday interface

`FaradayOperator` isolates the metric-free induction topology from the
Hodge-dependent current composition. `compute_rhs(E, rhs)` writes
`rhs = -d1(E)`, while `advance(E, B, dt)` performs `B -= dt*d1(E)`.
Both operations consume typed primal cochains and own neither geometry nor
communication policy.

Mercury explicit induction, the implicit resistive matrix action, and the
accepted implicit resistive update use this interface. Clearing and boundary
synchronization remain in Mercury; the former CT facade stays for legacy clients.

## 16. Unit 26: contract-driven owner synchronization

Mercury boundary groups now derive their owner-alias policy from the
registered FieldDescriptor contracts rather than a hard-coded list of group
names. A group uses that policy only when its name matches the declared sync
group; this preserves historical alias groups such as Eres1form. A mixed
non-None policy in one declared group now aborts during setup instead of
silently choosing a reduction rule.

Physical handlers, coupling, staged halo transfer, and final owner-alias
ordering are unchanged. This unit creates no new MPI calls.

## 17. Unit 27: catalog-driven boundary groups

The standard Mercury boundary groups are now created from `SyncGroups()`.
This removes duplicated field lists and coupling flags from the boundary
implementation while retaining the runtime group names used by physics.

The Ucell three-stage boundary schedule remains an explicit compatibility
override because it is broader than its ordinary halo contract. Eres1form
also remains explicit as a historical alias group. No communication order
or field membership changed.

## 18. Unit 28: weighted cell reconstruction operator

`WeightedCellReconstruction` owns the local cached-weight multiplication
from primal face 2-forms or edge 1-forms to Cartesian cell vectors.

Mercury Bface-to-Bcell, BplusBadd-to-Bcell, and Jedge-to-Jcell now use it.
Pole closures and boundary synchronization remain in Mercury. The operator
has no Grid, Field catalog, geometry build, boundary, or MPI dependency.

## 19. Unit 29: Pole least-squares constraint kernel

Pole current-cell reconstruction now calls the pure
`solve_regularized_least_squares_3` constraint kernel. It preserves the
previous diagonal regularization, partial pivoting, singular threshold, and
in-place normal-matrix behavior. Patch traversal and cell overwrite policy
remain in Mercury pending the next extraction unit.

## 20. Unit 30: Pole current-cell constraint composition

`PoleCurrentCellConstraint` now owns the complete local Pole reconstruction:
patch-order traversal, the two normal-axis sample patterns, edge-length
filtering, normal-equation accumulation, least-squares solve, and ring-cell
overwrite.  The type is fed only a storage-level patch description containing
the legacy field blocks and node bounds; it has no Grid, topology, Field
catalog, boundary, or MPI dependency.

Mercury retains the deliberately small adapter that selects `Pole` topology
patches in their original order and binds their legacy blocks.  Thus boundary
classification and storage binding stay application-specific while the
numerical constraint can be reused and tested independently.  The sample
order, regularization, pivot rule, singular threshold, and overwrite order
are unchanged.

## 21. Unit 31: snapshot-backed cubic-sphere current path

`GeometrySnapshot` now carries the two directional diagonal-Hodge coefficient
sets alongside primal and dual measures. `SnapshotLumpedHodge` consumes these
read-only views directly. Mercury retains one validated snapshot for its
lifetime and uses it for every regular `B -> J` evaluation, including the
implicit resistive delta-current path.

This removes per-call legacy geometry catalogue binding from the principal
cubic-sphere electromagnetic path. Cochain storage, the diagonal coefficient
values, stencil order, singular-edge override, and all synchronization remain
unchanged. The old Field-backed `LumpedHodge` remains as the compatibility
implementation and regression oracle while state storage is migrated.

## 22. State Unit 1: immutable layout and legacy catalog adapter

`mpcns::state::StateLayout` is the immutable description of quantity identity,
mathematical kind, staggered storage layout, ghost extent, allocation domain,
and synchronization requirement. It owns neither values nor MPI plans.
`mpcns::compat::LegacyStateLayoutBinding` is the only bridge from this layout
to the current `FieldCatalog` integer IDs.

The adapter preserves the exact existing catalog order and descriptors; the Z0
gate checks every registered quantity and its reverse legacy-ID binding. This
is the first State migration unit only: cochain views still borrow `FieldBlock`
storage, while later units introduce typed state handles and workspace lifetime
management without changing that backend.

## 23. State Unit 2: zero-copy typed StateView

`StateView` and `ConstStateView` are non-owning, strided views over the legacy
continuous allocation. They expose logical and inner extents, components, and
inline element access without importing Grid, Field, Halo, MPI, or Mercury.
The compat adapter performs the only `FieldBlock` conversion and validates the
immutable `StateLayout` contract before returning a view.

The explicit Mercury Faraday RHS path now binds cached `QuantityId`s once at
solver construction and invokes the StateView overload. Its stencil remains in
the shared incidence operator; no communication or boundary ordering changed.

Unit 2 closure moved the legacy `FieldBlock` cochain structs into
`mpcns::compat::legacy_cochain_view`; the public State cochain header now
contains only StateView-based types. Both legacy and StateView `d1` entry
points instantiate one implementation, preserving the historical component
expressions and floating-point evaluation order. Z0 checks all three
directions, multiple components, logical inner and ghost extents, const access,
invalid IDs, and a deliberately allocation-domain-inactive `FieldBlock`.

## 24. State Unit 3: minimal SimulationState and Workspace lifetime seam

`SimulationState` is a non-owning lifetime anchor for one immutable
`StateLayout`. It produces checked `StateHandle`s only; it neither allocates a
second state array nor owns a Field catalog. `Workspace` accepts only handles
whose `StorageRole` is `Workspace` and explicitly binds/release borrowed
`StateView`s. It copies only view metadata, never field values, and callers
must release/clear bindings before their storage backend is reconfigured or
destroyed.

The first Unit 3 slice intentionally has no new production-operator migration:
`LegacyStateViewAdapter` remains the sole FieldStorage adapter and existing
Mercury Faraday use remains unchanged. Z0 verifies direct workspace writes are
visible through the original FieldBlock, rejects non-workspace and unbound
handles, and confirms release semantics. The next State unit should migrate one
additional explicitly-scoped scratch consumer to this lifecycle API rather than
introducing an owning SimulationState.

## 25. State Unit 4A: role-separated, multi-block Workspace bindings

`QuantityKind` now contains only mathematical representation. `StorageRole`
independently distinguishes state, derived, workspace, and coupling-buffer
storage. The legacy descriptor has matching metadata; the catalog adapter
preserves its existing IDs, allocation, location, components, ghosts, and
synchronization contract while mapping the two axes separately. Mercury marks
only the existing local `Jstar_xi/eta/zeta` face-2-form scratch fields as
`Workspace` role in this unit. They remain `PrimalFaceForm`, and no current
operator has been migrated to consume them through `Workspace`.

`StateView` carries its immutable quantity ID and stagger location in addition
to its non-owning layout metadata. A `Workspace` binding is now keyed by
`(StateHandle, local_block_id)`. A replacement bind updates only that pair;
release likewise leaves other blocks intact, while clear removes all borrowed
bindings without touching FieldStorage. Handles are tied to a non-copyable,
non-movable `SimulationState`; bind validates owning state, workspace role,
quantity identity, stagger, and component count. These are setup/lifecycle
operations, not numerical hot-path calls.

Z0 covers simultaneous two-block binding, independent writes into the legacy
blocks, single-block release, clear, role rejection, mismatched-view rejection,
and cross-`SimulationState` handle rejection. Unit 4B should migrate only the
regular `Jstar` binding/setup portion of the current path to this API, leaving
the Hodge, `d1` transpose, owner/alias reconciliation, and all numerical loops
unchanged for its separate regression.

## 26. State Unit 4B: production Jstar Workspace binding

`StateView` now records the local process block index provided by the legacy
adapter. `Workspace::bind` verifies that this identity matches its binding key,
so a correct quantity view borrowed from another local block is rejected before
it can replace a binding.

MercurySolver owns one `SimulationState`, one `Workspace`, and three Jstar
handles. Construction creates them after the legacy layout/adapter and binds
every complete local Jstar triplet once. The setup checks that each borrowed
view aliases the exact `FieldBlock` allocation. Destruction order releases the
Workspace before the State/layout adapter; application-owned FieldStorage must
outlive the solver. Inactive complete triplets are skipped, while partial
triplets fail setup explicitly.

This is a production storage-lifetime migration, not yet a StateView numerical
operator migration. `SnapshotLumpedHodge`, `d1_transpose`, and
`ElectromagneticOperators` still consume legacy FieldBlock cochain views in
`ComputeRegularEdgeCurrent_`, for both ordinary B-to-J and implicit delta-B to
delta-J calls. Their formulas, field addresses, stencil order, alias handling,
and synchronization order are unchanged. The next State unit should add the
smallest compatible StateView Jstar bridge for this current operator only.

## 27. State Unit 4C-1: StateView snapshot star-2 production bridge

`SnapshotLumpedHodge::apply_star_2` now has a StateView cochain overload. The
legacy and StateView entries instantiate one local directional kernel: it uses
the same complete target extent, GeometrySnapshot diagonal coefficient view,
loop nesting, and `coefficient * source` expression. Jstar remains a dual
edge 1-form mathematically, despite its primal-face storage location.

Mercury caches the three validated Workspace views by local block at setup.
During both B-to-J and delta-B-to-delta-J current calls, source B views are
borrowed by legacy field ID without string lookup, snapshot star-2 writes those
cached Jstar StateViews, and the existing FieldBlock d1-transpose/inverse-star
continuation consumes the same allocation. `compute_current_from_star_2`
provides the sole shared downstream composition; star-2 is not recomputed.
No Workspace scan, bind, catalog build, or cochain allocation occurs in the
current hot loop.

StateView star-2 and continuation gates compare all three directions over the
full local extent and the final B-to-J output against the existing field-backed
path. A later Unit 4C-2 may migrate only d1-transpose's Jstar input to
StateView; inverse-star-1 and all communication remain outside that unit.

## 28. State Unit 4C-2: StateView d1-transpose input bridge

`apply_d1_transpose` now accepts either the legacy `DualEdgeFormView` or a
`DualEdgeFormStateView`, while retaining the FieldBlock-backed `DualFaceFormView`
output. Both public entries instantiate one directional implementation: the
xi/eta/zeta signs, negative-offset input reads, output-inner iteration ranges,
parentheses, and multiplication order are identical. Thus a dual edge 1-form
such as Jstar keeps its mathematical identity even though its borrowed storage
is located on the legacy primal-face allocation.

The StateView overload of `ElectromagneticOperators::compute_current_from_star_2`
passes that Jstar directly to the shared incidence entry and then invokes the
unchanged legacy inverse-star-1 stage. `MercurySolver::ComputeRegularEdgeCurrent_`
uses its setup-time cached Workspace views for both ordinary B-to-J and
delta-B-to-delta-J calls; it no longer obtains Jstar FieldBlocks in the hot
loop. Jedge remains the existing FieldBlock-backed dual-face intermediate and
primal-edge final output. No bind, lookup by name, allocation, MPI call, halo
operation, or owner/alias action was added to the numerical loop.

The Z0 DEC gate independently snapshots the legacy transpose output before
inverse-star-1, reruns it from the same zero-copy StateViews, and compares all
three directions across the full local allocation (including unchanged ghost
storage); maximum error is exactly zero. The next minimal unit, 4C-3, may add
a StateView output to inverse-star-1 while retaining its one existing
axis-safe numerical kernel and keeping all reconstruction, singular-edge,
and communication steps unchanged.

## 29. State Unit 4C-3: complete StateView local current chain

The production regular current composition is now entirely typed StateView
data access: induced (or delta) primal-face flux, Workspace-backed dual-edge
Jstar, StateView dual-face curl storage, and StateView primal-edge Jedge
output. Jedge deliberately aliases one existing FieldStorage allocation across
its intermediate dual 2-form and final primal 1-form identities. The legacy
adapter is the only FieldStorage boundary; Mercury no longer passes a
FieldBlock cochain into star-2, d1-transpose, or inverse-star-1.

`d1_transpose` has a fully StateView output overload. Its legacy, mixed, and
fully StateView entries instantiate one stencil with the same inner output
range and ghost-read offsets. Snapshot inverse-star-1 likewise has legacy and
StateView entries that call one component kernel; it retains `source *
coefficient` evaluation order and writes only the historical target inner
extent. A narrow virtual StateView bridge on `HodgeOperator` rejects unsupported
implementations explicitly; `SnapshotLumpedHodge` supplies the production
override, without a cast or a hierarchy rewrite.

The Z0 DEC gate seeds non-zero magnetic flux through the entire ghost extent,
then exactly compares legacy and StateView Jstar, curl, inverse-star-1 result,
and complete continuation for every direction and local block. It also repeats
the intentional in-place Jedge reuse. No temporary cochain, extra State array,
or hot-loop bind/allocation was introduced. The remaining FieldBlock uses in
the current workflow are storage-allocation checks, the explicit legacy/compat
operator entries, and the later singular-edge/owner-alias/boundary stages.

State Unit 5 should be an acceptance-only audit before any broader migration:
exercise a fixed Mercury case for ordinary B-to-J and implicit delta-B-to-delta-J,
compare B/J/divB and fluid monitors against the 4C-2 baseline, and inventory
the remaining FieldBlock consumers. It must not start Parallel or Physics
extraction by itself.

## 30. State Unit 5: acceptance audit and phase freeze

Acceptance was performed at git revision `45a9bec967cd9f616f646ea708397963388713fb`
with the Unix Makefiles generator, `/usr/bin/cc` and `/usr/bin/c++`, and no
configured build type. The exact successful build commands were
`cmake --build build --parallel 4`, `cmake -S Z0_NULL -B build/Z0_NULL`, and
`cmake --build build/Z0_NULL --target z0_null --parallel 4` (all exit 0).
The ten bundled Z0 modes were run independently and serially with three MPI
ranks: cell/node/edgeforms/faceforms/forms at vertex level; all at face, edge,
and vertex level; all with owner synchronization; and owner-only. Every command
exited 0. All report `D1*D0 == 0`, `D2*D1 == 0`, owner/alias consistency, the
Faraday RHS reference, StateView star-2, StateView d1 transpose, StateView
inverse-star-1, and full current continuation at exact maximum error 0.

The requested periodic-cube command,
`NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh`, was also
executed. Preprocess built only in `/tmp` and its source checkout was not
modified. The first all-vertex case printed every Z0/DEC assertion as PASS,
including the State checks, but one non-root MPI process returned exit 1 after
the root reported normal completion. The harness therefore stopped before its
owner modes. This is **not a passed periodic-cube validation**; it remains an
integration-harness limitation to reproduce and diagnose before treating cube
coverage as a Parallel-phase gate.

No runnable Mercury benchmark is stored in this checkout: `Z4_Mercury` has
reference controls but no `CASE/geometry`. Consequently no physical B, J,
divB, fluid-monitor, Pole, singular-edge, or implicit-resistive dB-to-dJ
baseline comparison was fabricated. Code inspection confirms both ordinary
and delta current calls use `ComputeRegularEdgeCurrent_`; the latter is called
from the implicit resistive path with `fid_dB` and `fid_dJ`. A supplied fixed
case must compare that shared path to `d86cc1d643e34b5bf39938838df9506444b09627`
at equal rank count, time, and controls.

The frozen public State seam is: immutable `StateLayout`/`QuantityDescriptor`,
non-owning typed `StateView`/`ConstStateView`, typed cochain views,
non-copyable `SimulationState` handles, and non-owning multi-block `Workspace`
bindings. `FieldStorage` remains the sole production allocation backend.
`LegacyStateLayoutBinding` is the one catalog mapping and
`LegacyStateViewAdapter` is the zero-copy borrowing boundary; State does not
duplicate the catalog or an array. Views retain the borrowed quantity, location,
and local block identity. Workspace validates all three and requires its
borrowed storage to outlive it; Mercury's setup binds/caches Jstar once before
the current hot loop.

There is one mathematical implementation each for StateView/legacy d1,
d1-transpose, snapshot star-2, and snapshot inverse-star-1. The StateView
inverse-star-1 virtual bridge fails explicitly for a non-supporting Hodge
implementation; only `SnapshotLumpedHodge` is used by the migrated current
path. No string lookup, Workspace bind, StateLayout construction, temporary
cochain allocation, MPI action, or physical-boundary action was added to that
path. Jstar remains a dual-edge 1-form despite primal-face storage; Jedge's
dual-face intermediate and primal-edge final interpretations remain explicit.

Remaining compatible technical debt is intentionally isolated: legacy
cochain adapter entries, allocation checks, reconstruction, singular-edge and
Pole handling, owner/alias reconciliation, physical boundaries, and the
unmigrated physics operators still use FieldBlock. The periodic-cube non-root
exit belongs to integration/Parallel follow-up, not a State arithmetic defect.

**State phase decision: A — State is frozen and may proceed to Parallel.**
This decision rests on exact local numerical regression and all ten standard
MPI modes. A real Mercury case and a passing generated-cube harness remain
mandatory acceptance evidence for future end-to-end or Parallel changes.



## 31. Parallel Unit 0: MPI acceptance baseline and periodic-cube exit diagnosis

**Root cause and minimal repair.** The retained periodic-cube fixture reproduced the original failure in all_vertex: rank 0 printed every assertion as PASS and normal completion, while ranks 1 and 2 returned exit 1. The new rank/module summary identified the first false module as DEC / State, not Topology, Physical IO, or Halo. Its failing assertion was SimulationState Workspace is multi-block, zero-copy, and safely keyed with max=inf on ranks 1 and 2. The cause was test-only: its fixture required field.num_blocks() >= 2 on every rank. The generated 3-rank periodic cube validly assigns only one local block to ranks 1 and 2. No owner/alias inconsistency, topology failure, halo failure, MPI-finalize failure, or post-normal-completion error was found.

The test now exercises one local block when that is all a rank owns, including zero-copy binding, local-block identity rejection, role/component/foreign-state rejection, release and clear. Where two local blocks exist it additionally retains the original independent multi-block/rebind assertions. This is not a threshold change or a skipped test; it makes the test storage-lifetime contract valid for both legal local block cardinalities.

**MPI conclusion rule.** Each rank retains four local module outcomes: Topology, Physical IO, Halo Communication, and DEC / State. After every module has run, every normal participant calls one MPI_Allreduce(MAX) through the existing wrapper, followed by the same-order MPI_Gather of the four failure bits. Root prints normal completion only after the all-reduced pass state is true. On failure root lists every failed rank and module; non-root assertion failures emit rank-specific diagnostics.

**Acceptance results (three MPI ranks; build parallelism at most four).**

| Command | Result |
| --- | --- |
| NP=3 BUILD_JOBS=4 KEEP_PREPROCESS_CUBE_FIXTURE=1 bash Z0_NULL/run_preprocess_cube_validation.sh | exit 0; all_vertex, all_owner_vertex, and owner_only all passed; retained fixture /tmp/mpcns-preprocess-cube.SEIvOE |
| NP=3 bash Z0_NULL/run_halo_validation.sh | serial per mode and consumes the MPI program exit status; cell, node, edgeforms and an independent faceforms/vertex re-run passed with the global summary. The remaining modes require one uninterrupted full-harness rerun before claiming the ten-mode gate complete. |

The retained failing fixture before repair was /tmp/mpcns-preprocess-cube.mujwpZ; its all_vertex log records the original non-root exit. The post-repair retained fixture above is the minimum reproducer and passing evidence. Preprocess source was not changed.

**Read-only communication map for Parallel Unit 1.**

| Area | Current responsibility | Unit 1 placement |
| --- | --- | --- |
| Topology, TopologyEquivDetail, TopologyBuilder | Builds local/parallel patches; reconciles node, edge and face equivalence; selects owner and orientation sign | mpcns::mesh ownership/equivalence data |
| FieldDescriptor, FieldSyncContract, FieldHaloRequest, FieldCatalog | Declares location, representation, group, halo level, owner policy and orientation awareness | durable contract types; keep compatibility adapters |
| Halo::build_registered_patterns plus HaloBuildFace/Edge/Vertex | Converts topology patches to face/edge/vertex send/receive boxes and per-neighbor tags | candidate mpcns::parallel planning layer |
| HaloExchangeParallel, HaloTriplet, HaloCouplingParallel | Pack, nonblocking send/receive, wait, unpack component and oriented form transfers | candidate mpcns::parallel execution layer |
| HaloEdgeOwnerSync | Owner-to-alias local and MPI synchronization with orientation sign | split plan/execution later; preserve compat entry points now |
| Z0_Boundary and legacy Halo public API | Existing production call order and boundary/coupling orchestration | must remain compat during Unit 1 |

The numerical compatibility constraints are: stage order is face (1D), edge (2D), vertex (3D); registered groups/channels are processed serially; form triplets preserve transform permutation and sign; sender/receiver tags agree per neighbor and stage; owner-alias synchronization occurs only after the relevant halo data is present. Same-rank multi-block paths and cross-rank patches both remain required coverage. Unit 1 should extract read-only mesh and planning representations first, without moving these production calls or altering tag/order semantics.


## 32. Parallel Unit 1: read-only MeshPartition acceptance and freeze

MeshPartition is the immutable, MPI-free description of one rank local blocks
and face/edge/vertex adjacency. Its public data uses only LogicalIndex and
LogicalBox. BlockAddress is always (rank, rank-local block index); no global
block number is invented. Neighbor indices are validated in their own rank
namespace and never against the current rank block count.

A patch belongs to this_block; neighbor_block is a topological neighbor, not a
communication direction. Inner remains same-rank, Parallel remains cross-rank,
and Physical has an invalid neighbor and nullopt transform. Coupling retains
the legacy direction neighbor to this_block. Patch order is exactly legacy
Inner, Parallel, then Physical order within each codimension, without sorting
or deduplication. The adapter copies boxes and stored permutation/sign/offset
verbatim and never recomputes topology, owner, quotient IDs, or orientation.
PhysicalPatch::raw, MPI tags, buffers, and requests do not enter MeshPartition.
The scoped const_cast in the adapter is read-only compatibility for legacy
Array1D, which has no const index operator.

Z0 independently compares every local block and every face, edge, and vertex
field: identities, names, six box bounds, transform permutation/sign/offset,
directions, coupling, and physical metadata. Diagnostics include rank,
codimension, kind, index, both blocks, field, expected, and actual. Constructor
tests cover invalid rank/size, block ordering, local patch identity, incomplete
neighbors, Physical contradictions, missing transforms, and Inner/Parallel
rank relations. Legal single/multi-block, repeated periodic, same-rank,
cross-rank, and physical inputs remain accepted.

Acceptance evidence on three ranks:
- cmake -S . -B build and cmake --build build --parallel 4: exit 0; MercuryZ4 links separate mpcns_mesh and mpcns_legacy_mesh targets.
- Z0 configure and z0_null build with parallel 4: exit 0.
- Ten standard modes each exited 0 in serial order; every log reports checked=9229 failed=0, global MPI pass, and normal completion.
- NP=3 BUILD_JOBS=4 periodic-cube validation: exit 0 for all three modes; each reports checked=4079 failed=0.

Topology equivalence, owner/alias maps and GIDs, collapsed/singular entities,
Halo patterns, tags, pack/send/recv/unpack, synchronization order, and all
numerical operators remain authoritative in legacy Topology/Halo. Parallel
Unit 1 is frozen. Unit 2A may consume MeshPartition, field sync contracts, and
legacy Halo pattern observations for a read-only ExchangePlan description, but
must not yet replace execution.

## 33. Parallel Unit 2A: read-only ExchangePlan acceptance and freeze

`mpcns::parallel::ExchangePlan` is an immutable, MPI-free description of one
component-copy pattern key. The plan key is `FieldLocation` plus ghost depth;
each ordered `ExchangeRegion` records stage, one of the six explicit
inner/parallel roles, local and peer `BlockAddress`, role-valid logical boxes,
the required permutation/sign/offset transform, role-valid optional MPI tags,
and its index in the source legacy pattern. The public header depends only on
core logical indices and MeshPartition. It contains no Grid, Topology, Halo,
Field storage, MPI request/status, communicator, or Mercury parameter type.

MeshPartition remains the authority for block identity and topological
adjacency. ExchangePlan describes only communication work that the legacy Halo
builders actually retained. `LegacyHaloPlanAdapter` is the compatibility
boundary: it reads a copy returned by the const
`Halo::debug_component_copy_patterns` diagnostic and never derives boxes,
transforms, tags, or filtered regions from MeshPartition. The diagnostic copies
already-built containers in production traversal order, performs no MPI,
does not rebuild a pattern, and exposes no mutable internal reference.

The exact legacy mapping is:

| Legacy container | Exchange kind | Valid boxes | Valid tags |
| --- | --- | --- | --- |
| inner face/edge/vertex | InnerCopy at its source stage | send and recv | none |
| parallel face | ParallelFace | send and recv | send and recv |
| parallel edge send / recv | ParallelEdgeSend / ParallelEdgeRecv | send / recv respectively | send / recv respectively |
| parallel vertex send / recv | ParallelVertexSend / ParallelVertexRecv | send / recv respectively | send / recv respectively |

This distinction preserves the historical face pattern's bidirectional local
view and the independently negotiated edge/vertex send and receive views.
Edge/vertex transforms are copied from their respective actual containers, so
the inverse sender view is not reconstructed from the receiver. Empty
role-inapplicable boxes and uninitialized inner tags are never read. Pattern
keys, special-interface filtering, vector order, and assigned tags remain
legacy-authoritative.

Z0 constructs the plan only after `build_registered_patterns` has completed.
Its expected side reads each `LegacyHaloPatternSnapshot` directly and checks
location, ghost depth, stage, kind/role, both addresses, every valid box bound,
all transform permutation/sign/offset entries, every valid tag, per-pattern
region index, total size, and flattening order. It also verifies explicit
rejection of an absent key, invalid block address, missing transform, and an
orientation-aware request. Failures join the existing Halo Communication
module and existing four-module MPI global summary; no new collective summary
or runtime buffer was introduced.

Acceptance used three MPI ranks and build parallelism four or less:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built and linked |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; all ten modes passed serially |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; all_vertex, all_owner_vertex, owner_only passed serially |

Across the ten standard modes, 72 rank-local plans and 6,032 actual regions
produced 198,488 exact checks with zero differences. The generated periodic
cube added 36 plans, 1,704 regions, and 55,821 checks with zero differences.
Every cube log reports the MPI-global pass conclusion and normal completion.
Existing DEC, StateView, Hodge, Faraday, owner/alias, and MeshPartition gates
also remained passing.

This unit supports only scalar and component-wise Cartesian/conserved copying.
Orientation-aware edge 1-form and face 2-form triplets reuse some geometric
patterns but still require component grouping and sign/permutation application.
Directed coupling needs source/destination channel identity and its separate
pattern family. Node/edge/face owner-alias synchronization needs equivalence
owner/alias identities, reduce/broadcast role, and orientation sign. None is
represented as supported ExchangePlan work here. Legacy Halo remains solely
responsible for registration, pattern construction, buffers, tags, MPI calls,
pack/unpack, waits, face-to-edge-to-vertex ordering, forms, coupling, and
owner/alias execution.

Parallel Unit 2A is frozen. A later Unit 2B may define read-only descriptions
and parity gates for one currently unsupported family at a time, beginning only
after its extra grouping, direction, and orientation inputs are explicit. It
must not make ExchangePlan execute production communication implicitly.


## 34. Parallel Unit 2B-1: Face-level oriented-form transfer-plan acceptance

`mpcns::parallel::OrientedFormTransferPlan` is the immutable, MPI-free
description of one orientation-aware three-member form group at the Face
stage. A plan names the sync group and its Xi/Eta/Zeta quantities, form kind,
location, axis, component count, and ghost depth. Each ordered transfer records
Inner or Parallel kind, actual source and destination `BlockAddress`, source
and destination axes and quantities, source and destination legacy region
indices, send and receive boxes, the stored transform, orientation sign, and
role-valid tags. Packing order is explicitly transformed destination
i/j/k/component order. The record does not reinterpret the triplet as a
Cartesian vector.

The legacy mathematical convention is frozen rather than reformulated. For
destination axis `d`, source axis `s` is the unique axis satisfying
`transform.permutation[s] == d`. Edge 1-form sign is
`transform.sign[s]`. Face 2-form sign is
`sign[0] * sign[1] * sign[2] * sign[s]`, exactly matching
`HALO_TOOLS::face_2form_orientation_sign`; every stored and computed sign is
validated as plus or minus one. Inner records describe the directed local copy
from the source pattern's `this_block` to the destination pattern's
`neighbor_block`. Parallel records describe the local outbound source and
peer destination; their receive box is the matching local inbound box in the
same bidirectional Face pattern. Topological this/neighbor identity therefore
is never used as a generic send/receive inference.

`LegacyOrientedFormPlanAdapter` is the only conversion boundary. It observes
the already registered triplet group and already built Face `HaloPattern`
snapshots through const diagnostics. It copies field order, region order,
boxes, transforms, and assigned tags; it performs no pattern construction, tag
allocation, MPI call, synchronization, numeric field read, or runtime-buffer
allocation. Plan construction is initialization/test-only. The request-based
entry point accepts only `do_halo && orientation_aware` EdgeCovariant1Form or
FaceContravariant2Form requests. A non-orientation-aware request is rejected
from this plan and remains eligible for the Unit 2A component-copy adapter.
The group-only conversion is private, so callers cannot bypass that request
classification. An additional three-rank `Z0_FORM_ORIENTATION=0` probe exited
zero with 18 Unit 2A component-copy plans and zero Unit 2B-1 plans.
Missing members, mixed kinds, inconsistent component/ghost metadata, invalid
axes/permutations/signs, rank-role errors, absent tags, and inconsistent
source/destination pattern identity fail explicitly.

Z0 independently reads the registered group and the actual legacy Face
patterns rather than comparing two adapter products. It checks group and
member identity, location, both axes and quantities, both block addresses,
all box bounds, all permutation/sign/offset entries, both tags, shape,
per-region indices, and the exact legacy traversal order: all Inner axes and
regions, then all Parallel axes and regions. Distinct rank/block/axis/component
markers are synchronized by the unchanged production Halo. A shadow oracle
replays directed same-rank writes in legacy order and then applies peer
receives, so overlapping regions are judged by their real final writer.
Separate pure mapping cases guarantee non-identity permutation and both
one-form and two-form negative-sign detection even when a mesh contains only
identity transforms.

Acceptance used three MPI ranks and build parallelism four or less:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built and linked |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; all ten modes passed serially |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; all_vertex, all_owner_vertex, and owner_only passed serially |

Across the ten standard modes, 42 rank-local plans and 2,856 transfers produced
4,515,552 exact checks with zero failures. They include 1,008 Inner and 1,848
Parallel transfer records, 4,620 non-identity permutation-axis observations,
1,596 negative orientation signs, and 4,393,536 final numeric value checks:
1,957,900 same-rank and 2,435,636 cross-rank. The periodic cube added 18 plans,
540 transfers, and 52,191 exact checks with zero failures, including 7,716
same-rank and 21,156 cross-rank final-value checks. Its 180 Inner and 360
Parallel records have identity transforms and positive signs, so nontrivial
orientation coverage comes from the standard topology and 18 synthetic
mapping/sign assertions. Every mode also retained passing ExchangePlan,
MeshPartition, topology, owner/alias, physical I/O, DEC, Faraday, Hodge, and
StateView results and the MPI-global normal-completion conclusion.

No production pattern, tag, pack/unpack, send/receive, wait, buffer, stage
order, or numerical operator was changed. Legacy Halo remains the only
executor. Edge- and Vertex-stage form paths remain deliberately unfrozen:
their parallel functions construct `TripletCornerMeta` during synchronization,
exchange it with `MPI_Alltoallv`, and create temporary data buffers and MPI
requests before pack/send/receive/unpack. Moving those paths requires a later
initialization-time protocol that preserves directed metadata pairing,
send/receive transforms and boxes, peer ordering, assigned tags, buffer
lifetimes, and rebuild invalidation. Directed coupling and owner/alias
reduce/broadcast are also outside this unit.

Parallel Unit 2B-1 is complete and frozen. The next smallest unit should
describe and compare one Edge- or Vertex-stage triplet metadata family at a
time; it must not introduce an ExchangeExecutor or replace production Halo.


## 35. Parallel Unit 2B-2A: Edge-stage Edge 1-form metadata snapshot acceptance

`mpcns::parallel::OrientedCornerTransferPlan` is the immutable, MPI-free
description of the cross-rank Edge-stage portion of one orientation-aware
edge 1-form triplet. Each receive or send record owns its source and
destination `BlockAddress`, source and destination axes and quantity names,
destination receive box, source node box, two source-side directions, stored
receive-to-send transform, orientation sign, peer rank, message tag, shape,
message length, local plan order, per-peer order, and the frozen packing order
`DestinationRecvIJKComponent`. Send records additionally retain the actual
source send box computed by legacy Halo; receive records expose no synthetic
send box. Valid empty half-open receive/send boxes and zero-length messages
are retained because they are present in the legacy protocol.

The Edge-stage convention is intentionally distinct from the frozen Face
plan. For destination/receive axis `d`, the source/send axis is directly
`transform.permutation[d]`, and the edge 1-form sign is directly
`transform.sign[d]`. Packing traverses the destination receive box in
i/j/k/component order, maps each destination index through the stored
receive-to-send transform, and applies that sign. The computed source send box
is diagnostic metadata; the unchanged legacy packer is driven by the mapped
destination traversal. Pure contract cases include a non-self-inverse
permutation, so accidentally applying the Face-stage inverse-axis rule is
detected independently of mesh coverage.

The compatibility snapshot is opt-in and disabled by default. When enabled,
`Halo::exchange_parallel_edge_edge_1form_triplet_` copies the actual receive
metadata built from parallel edge patches and the actual send metadata
delivered by the existing `MPI_Alltoallv`, after the existing source box and
message length have been computed. It retains values only: no field pointer,
buffer, request, or mutable pattern reference escapes. Capturing adds no MPI
call, tag negotiation, buffer, or communication-stage change. Snapshots are
cleared whenever registered patterns are rebuilt, and disabling capture does
not affect synchronization. `LegacyOrientedCornerPlanAdapter` only validates
and converts this captured evidence; it never reconstructs the cross-rank
pairing or reads numerical field storage.

Z0 derives receive expectations independently from raw parallel
`TOPO::EdgePatch` traversal. A test-only all-gather of those expected records
constructs the remote sender view in the same rank-grouped, stable order as
the legacy all-to-all metadata exchange. It then compares every identity,
axis, quantity, box bound, direction, transform entry, sign, peer, tag,
length, and both ordering indices. Distinct rank/block/axis/index/component
markers are packed by an independent oracle and compared with the unchanged
production synchronization. The oracle first replays the legacy same-rank
Inner Edge writes so that cross-rank sends observe the same pre-send field
state, and resolves overlapping receives by actual receive-plan final-writer
order. Capture-default-off and rebuild-invalidation behavior are also tested.

Acceptance used three MPI ranks and build parallelism four or less:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built and linked |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; all ten modes passed serially |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; all_vertex, all_owner_vertex, and owner_only passed serially |

Across the ten standard modes, 18 rank-local captured plans contained 792
receive and 792 send records. All 792 global sender/receiver pairs matched,
195,846 exact metadata/contract checks and 55,776 final numeric value checks
had zero failures. Coverage included 1,440 non-identity axis mappings, 504
negative signs, 48 retained zero-length messages, all 18 active rank/mode
instances with multiple local blocks, and 30 synthetic direction cases. The
periodic cube added nine plans, 180 receive and 180 send records, 180 matched
pairs, 36,201 exact checks, and 4,224 final-value checks with zero failures.
Its transforms were identity and positive; the standard topology and pure
mapping cases supply the required nontrivial orientation coverage. Every log
retained the existing MPI-global pass and normal-completion conclusion, and
the ExchangePlan, MeshPartition, DEC, Faraday, Hodge, StateView, and
owner/alias gates remained passing.

No production topology builder, pattern builder, tag, pack/unpack formula,
send/receive/wait order, buffer ownership, stage order, numerical operator,
coupling path, or owner/alias path was changed. Legacy Halo remains the only
executor and still creates this metadata during synchronization. Vertex-stage
edge 1-form metadata, Edge/Vertex face 2-form metadata, coupling,
owner/alias, and an ExchangeExecutor remain out of scope. A later Unit 2B-2B
may freeze the Vertex-stage family, but first must preserve its three-direction
metadata construction, all-to-all pairing, mapped destination traversal,
zero-length behavior, tag/order contract, and rebuild invalidation without
moving production communication.

Parallel Unit 2B-2A is complete and frozen.


## 36. Parallel Unit 2B-2B: Vertex-stage Edge 1-form metadata snapshot acceptance

`mpcns::parallel::OrientedVertexTransferPlan` is the immutable, MPI-free
description of the cross-rank Vertex-stage portion of one orientation-aware
edge 1-form triplet. It reuses the frozen receive/send role and destination
i/j/k/component packing-order identities from the Edge-stage plan, but keeps a
separate transfer and plan type. Each record owns the actual source and
destination block addresses and axes, quantity identities, orientation sign,
3D destination receive box, source node box, all three source directions,
receive-to-send transform, peer, tag, shape, message length, and stable local
and per-peer order. Only Send records expose the actual source send box
computed by legacy Halo. The public types retain no MPI, Field, Halo,
Topology, buffer, request, or legacy metadata pointer.

Edge and Vertex share the direct legacy axis/sign convention:
`send_axis = transform.permutation[recv_axis]` and
`sign = transform.sign[recv_axis]`. Vertex is not an alias for Edge:
it consumes `TOPO::VertexPatch`, constructs 3D corner boxes, and preserves
three directions. The first direction uses the legacy normal mapping while the
second and third use the legacy tangential mapping. If the mapped normal is
not initially in the first slot, the unchanged rule swaps it with the second
or third slot whose absolute axis equals
`permutation[abs(dir1)-1] + 1`. All three resulting absolute axes must be
valid and distinct. No inverse-axis Face-stage convention is substituted.

The legacy receiver traverses filtered non-coupling parallel Vertex patches
in source order and appends recv axes 0, 1, 2. Tags start at zero and increase
independently per peer. These receiver records are grouped by destination
peer; the existing `MPI_Alltoall` plus `MPI_Alltoallv` delivers sender
records grouped by origin rank while retaining within-rank order. The snapshot
copies the actual local receive records and the actually delivered send
records after the 3D send box and message length are known. It adds no second
metadata exchange and does not alter a tag, buffer, request, send, receive,
wait, pack, unpack, or stage-order operation.

Packing remains driven by the destination receive box in
i/j/k/component order. Every destination index is mapped through the stored
receive-to-send transform; an out-of-bounds source contributes zero, and a
non-finite in-bounds source also contributes zero. The source send box is
retained as diagnostic metadata but does not replace this mapped destination
traversal. Receive records are unpacked in their actual vector order, so the
last real receive remains authoritative for overlaps.

Vertex capture is opt-in and disabled by default. Its enable flag and snapshot
map are independent of the frozen Edge-stage capture. Values, not temporary
pointers, are retained. Enabling Vertex capture clears only prior Vertex
snapshots; disabling it leaves the captured result readable. Rebuilding
registered patterns invalidates both stage-specific snapshot families because
both describe topology-derived work that is no longer current.
`LegacyOrientedVertexPlanAdapter` only validates and converts the captured
evidence. It does not reconstruct sender/receiver pairing from MeshPartition
or inspect numerical field values.

Z0 derives expected receiver metadata independently from raw
`TOPO::VertexPatch` traversal and performs a test-only all-gather to form the
expected sender view in the same rank-grouped order as legacy metadata
exchange. Every identity, axis, field, direction, box bound, transform entry,
sign, peer, tag, length, and ordering index is compared exactly. The numeric
oracle seeds rank/block/axis/index/component-distinct values, runs the real
Face and Edge stages first, snapshots that resulting state, independently
replays same-rank Vertex writes, and predicts the cross-rank payload and final
writer for every destination. A separate all-NaN source probe verifies the
production non-finite-to-zero rule. Pure contract cases cover a
non-self-inverse permutation, negative sign, three-direction reordering,
duplicate-direction rejection, invalid sign/permutation/rank/role rejection,
missing Send box, and acceptance of a valid empty box with zero message
length.

Acceptance used three MPI ranks and build parallelism four or less:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built and linked |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; all ten modes passed serially |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; all_vertex, all_owner_vertex, and owner_only passed serially |

Across the ten standard modes, 18 rank-local plans contained 144 receive and
144 send records; all 144 global pairs matched. The suite made 29,118 exact
metadata/contract checks, 1,152 finite mapped ghost checks, and 1,152
non-finite-source zero checks with no failures. Actual topology coverage
included 252 non-identity permutation-axis entries, 120 negative signs, 108
three-direction reorder records, 96 out-of-bounds source components, and all
18 active rank/mode instances with multiple local blocks. The standard
topology produced no zero-length Vertex messages, so legal zero-length
acceptance is covered by the pure contract case rather than claimed as real
mesh coverage.

The periodic cube added nine plans, 72 receive and 72 send records, 72 matched
pairs, 14,481 exact checks, 576 finite mapped ghost checks, and 576
non-finite-source zero checks, all without failure. Its Vertex transforms were
identity and positive, and it produced no direction reorder, out-of-bounds, or
zero-length record; those cases are supplied by the standard topology and the
13 pure contract cases. Every mode retained MPI-global pass and normal
completion, and the frozen MeshPartition, ExchangePlan, Face oriented-form,
Edge metadata, DEC, Faraday, Hodge, StateView, and owner/alias regressions
remained passing.

No production MPI timing or metadata exchange was moved. Initialization-time
prebuilding would still require finalized topology and triplet registration,
a stable communicator/rank assignment, preservation of rank-grouped
all-to-all sender ordering and per-peer tags, durable 3D boxes/transforms and
buffer lengths, and rebuild-safe invalidation tied to every topology or
registration change. Those prerequisites are recorded only; legacy Halo
remains the sole executor.

Face 2-form Edge/Vertex metadata, coupling, owner/alias communication, and an
ExchangeExecutor remain out of scope. Parallel Unit 2B-2B is complete and
frozen; no production metadata negotiation migration was started.


## 37. Parallel Unit 2B-3A: Edge-stage Face 2-form metadata snapshot acceptance

`mpcns::parallel::Face2FormEdgeTransferPlan` is the immutable, MPI-free
description of the cross-rank Edge-stage portion of one orientation-aware face
2-form triplet. It reuses the frozen receive/send role and destination
i/j/k/component packing-order identities, together with BlockAddress,
LogicalBox, IndexTransform, HaloStage, and FormQuantityIdentity. Each record
owns the actual source and destination blocks, axes and FaceXi/FaceEta/FaceZeta
quantity identities, orientation sign, destination receive box, source node
box, two source directions, receive-to-send transform, peer, tag, shape,
message length, and stable local and per-peer order. Only Send records expose
the actual source send box. The public types retain no MPI, Field, Halo,
Topology, buffer, request, or legacy metadata pointer.

The axis rule is deliberately different from the frozen Edge1Form plans. For
destination axis `d`, Face2Form uses the unique source axis `s` satisfying
`transform.permutation[s] == d`; it never uses
`transform.permutation[d]`. The orientation sign is the unchanged legacy
rule `sign[0] * sign[1] * sign[2] * sign[s]`, normalized to plus or minus
one. A non-self-inverse synthetic permutation makes the direct and inverse
answers differ and verifies a negative 2-form sign. The orientation-unaware
legacy axis-triplet path retains sign one and is excluded from this plan.

Receiver metadata is built from filtered non-coupling parallel Edge patches in
legacy order, with recv axes 0, 1, 2 appended for each patch. Tags begin at
zero and increase independently per peer. Existing `MPI_Alltoall` and
`MPI_Alltoallv` operations group records by peer rank and deliver sender
records in origin-rank order while preserving within-rank order. The optional
snapshot copies the actual receive plan and the actually delivered send plan
after the 2D source send box and message length are known. It performs no
second metadata exchange and does not alter tag assignment, communicator,
buffer ownership, or point-to-point ordering.

Packing remains driven by the destination receive box in
i/j/k/component order. Each destination coordinate is mapped through the
stored receive-to-send transform. Out-of-bounds source coordinates contribute
zero, as do non-finite in-bounds source values; otherwise the exact 2-form sign
is applied. The source send box is diagnostic metadata and does not replace
the mapped destination traversal. Receive records are unpacked in their actual
vector order, retaining the real final writer for overlapping destinations.

Capture is independent, opt-in, and disabled by default. It engages only when
both its flag and the target function's `oriented` argument are true.
Enabling it clears only prior Face2Form Edge snapshots; disabling it leaves a
captured value snapshot readable. Registered-pattern rebuild invalidates it
along with the previously frozen stage snapshots. A direct execution of the
legacy `oriented=false` axis-triplet path on all three ranks produced no
snapshot. `LegacyFace2FormEdgePlanAdapter` only validates and converts
captured evidence and never reconstructs MPI pairing from MeshPartition or
reads field values.

Z0 independently traverses raw `TOPO::EdgePatch` data, computes inverse axes
and legacy 2-form signs, and uses a test-only all-gather to construct expected
sender records in the actual rank-grouped order. It compares every rank,
block, axis, quantity, direction, box bound, transform entry, sign, peer, tag,
length, and ordering index. The numeric oracle seeds
rank/block/axis/index/component-distinct values, executes the real Face stage,
snapshots that state, independently replays same-rank Inner Edge writes, and
predicts every cross-rank payload and final writer. A separate all-NaN probe
verifies non-finite-to-zero behavior. Pure contract tests cover inverse versus
direct permutation, negative and incorrect 2-form signs, invalid
permutation/rank/direction/role, missing Send box, and legal empty half-open
boxes with zero length.

Acceptance used three MPI ranks and build parallelism four or less:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built and linked |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; all ten modes passed serially |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; all_vertex, all_owner_vertex, and owner_only passed serially |

Across the ten standard modes, 18 rank-local plans contained 792 receive and
792 send records; all 792 global pairs matched. The suite made 250,488 exact
metadata/contract checks, 55,104 finite mapped ghost checks, and 55,104
non-finite-source zero checks with no failures. Actual topology coverage
included 1,440 non-identity permutation-axis entries, 384 negative 2-form
signs, 96 zero-length messages, 5,544 out-of-bounds source components, and all
18 active rank/mode instances with multiple local blocks. The unchanged
orientation-unaware path was executed and rejected from capture on all 18
active rank/mode instances.

The periodic cube added nine plans, 180 receive and 180 send records, 180
matched pairs, 40,026 exact checks, 3,984 finite mapped ghost checks, and 3,984
non-finite-source zero checks, all without failure. Its transforms were
identity and positive and produced no out-of-bounds or zero-length message;
the standard topology and 11 pure contract cases provide those nontrivial
conditions. Every log retained MPI-global pass and normal completion, and the
frozen MeshPartition, ExchangePlan, Face-stage oriented form, Edge1Form
Edge/Vertex metadata, DEC, Faraday, Hodge, StateView, and owner/alias
regressions remained passing.

No production MPI timing or metadata negotiation was moved. A later
Vertex-stage Face2Form unit may reuse the value-snapshot lifecycle,
rank-grouped pairing validation, inverse-axis/sign checks, and numeric oracle,
but must add the third direction and its existing reorder rule plus 3D corner
boxes; it must observe the separate
`exchange_parallel_vertex_face_2form_triplet_` path rather than extending
this two-direction contract. Initialization-time prebuilding still requires
finalized topology and triplet registration, stable communicator/rank
assignment, preserved peer grouping and per-peer tags, durable boxes and
transforms, planned buffer lifetimes, and invalidation on every relevant
topology or registration rebuild.

Coupling, owner/alias communication, an ExchangeExecutor, and production
metadata migration remain out of scope. Parallel Unit 2B-3A is complete and
frozen.

## 38. Parallel Unit 2B-3B: Vertex-stage Face 2-form metadata snapshot acceptance

`mpcns::parallel::Face2FormVertexTransferPlan` is the immutable, MPI-free
description of the cross-rank Vertex-stage portion of one orientation-aware
face 2-form triplet. It is deliberately separate from the frozen Edge-stage
Face2Form type. Each record owns source and destination `BlockAddress`
values, inverse-permutation source and destination axes, stable
FaceXi/FaceEta/FaceZeta quantity identities, the legacy 2-form orientation
sign, the actual 3D destination receive box, optional actual source send box,
source node box, all three mapped source directions, the receive-to-send
transform, peer, tag, shape, message length, local order, and per-peer
message order. Only Send records expose a source send box. Empty half-open
boxes and zero message lengths remain legal. The public type contains no MPI,
Field, Halo, Topology, request, buffer, or legacy metadata pointer.

The Face2Form rule remains
`send_axis = inverse_permuted_axis(transform, recv_axis)`, followed by the
unchanged sign
`sign[0] * sign[1] * sign[2] * sign[send_axis]`. Vertex does not inherit the
direct-axis Edge1Form convention. The first Vertex direction uses the legacy
normal mapping, while the second and third use the tangential mapping. If the
mapped normal is not initially first, the unchanged legacy rule swaps it with
the direction whose absolute axis is
`permutation[abs(dir1)-1] + 1`. The public contract requires all three
resulting absolute axes to be valid and distinct. A pure non-self-inverse
mapping test independently distinguishes inverse from direct permutation,
checks a negative 2-form sign, and exercises three-direction reordering.

Receiver metadata is observed from the actual non-coupling parallel
`TOPO::VertexPatch` traversal, with destination axes 0, 1, 2 appended in
legacy order and tags increasing independently per peer. The existing
metadata `MPI_Alltoall`/`MPI_Alltoallv` remains the only source of the
sender view. The optional snapshot copies the actual local receive plan and
the actually delivered sender plan only after legacy Halo has computed each
3D source send box and message length. It adds no metadata exchange and
changes no communicator, tag, allocation, pack, send, receive, wait, unpack,
or stage-order operation.

Packing remains destination-driven in receive-box i/j/k/component order.
Each destination coordinate is mapped through the stored transform. An
out-of-bounds source or a non-finite in-bounds source contributes zero;
otherwise the exact 2-form sign is applied. The source send box is retained
for diagnostics rather than substituted for that traversal. Receive records
are unpacked in their actual vector order, so the last real writer remains
authoritative for overlaps.

Capture is optional and disabled by default. It is both stage-specific and
form-specific and engages only when its flag and the target function's
`oriented` argument are true. Enabling it clears only prior Vertex Face2Form
snapshots; disabling it leaves the copied values readable. Rebuilding
registered patterns invalidates the snapshot. Z0 explicitly verifies that
Vertex Face2Form capture creates no Edge Face2Form, Edge Edge1Form, or Vertex
Edge1Form snapshot, and that executing the real `oriented=false`
Vertex-level axis-triplet path on every active rank creates no snapshot.
`LegacyFace2FormVertexPlanAdapter` only validates and converts captured
evidence; it does not reconstruct topology or MPI pairing and does not read
field values.

The independent Z0 oracle traverses raw `TOPO::VertexPatch` data, derives
all three directions, inverse axes, 2-form signs, receive boxes, per-peer tags,
and receiver order without consulting the adapter. A test-only all-gather
forms the expected sender view in the same rank-grouped order as legacy
metadata exchange. Every rank, block, axis, field, direction, box bound,
transform entry, sign, peer, tag, length, and ordering index is compared
exactly. The numeric oracle seeds rank/block/axis/index/component-distinct
values, runs the real Face and Edge prestages, snapshots that state,
independently applies same-rank Vertex writes, then predicts every cross-rank
payload and final destination writer. A separate probe runs the same two
prestages, replaces the sources with NaNs, and verifies production
non-finite-to-zero behavior at Vertex destinations.

Acceptance used three MPI ranks, build parallelism four, and serial harness
execution:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built and linked |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| focused `Z0_TEST_DRIVEN_SYNC=1 Z0_HALO_MODE=faceforms Z0_HALO_LEVEL=vertex mpirun -np 3 ./z0_null` | 0; metadata and numeric oracle passed |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; all ten modes passed serially |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; all_vertex, all_owner_vertex, and owner_only passed serially |

Across the ten standard modes, 15 rank-local plans contained 120 receive and
120 send records, and all 120 sender/receiver pairs matched. The suite made
24,420 exact metadata/contract checks, 960 finite mapped ghost checks, and
960 non-finite-source zero checks with no failures. Actual standard-topology
coverage included 210 non-identity permutation-axis entries, 40 negative
2-form signs, 90 reordered three-direction records, 180 out-of-bounds source
components, and all 15 active rank/mode instances with multiple local blocks.
The unchanged orientation-unaware path was executed and excluded from capture
on all 15 active instances.

The periodic cube added nine plans, 72 receive and 72 send records, 72
matched pairs, 14,526 exact checks, 576 finite mapped ghost checks, and 576
non-finite-source zero checks, all without failure. Its Vertex transforms were
identity and positive and produced no direction reorder, out-of-bounds
source, or zero-length record. The standard topology supplies the nontrivial
transform, sign, reorder, and out-of-bounds coverage. Neither real topology
produced a zero-length Vertex Face2Form message, so zero-length acceptance is
covered by the pure contract case and is not claimed as observed mesh
coverage. Fourteen pure contract checks cover inverse/direct distinction,
negative and incorrect 2-form signs, direction reordering and duplication,
invalid axes/permutation/rank/sign/role, a missing Send box, and a legal empty
message.

Every log retained the MPI-global pass and `Z0_NULL finished normally`
conclusion. MeshPartition, component-copy ExchangePlan, Face-stage oriented
forms, Edge1Form Edge/Vertex snapshots, Face2Form Edge snapshots,
owner/alias, DEC, Faraday, Hodge, and StateView regressions remained passing.
No production communication path, numerical formula, topology builder, tag,
buffer, synchronization order, coupling path, or owner/alias operation was
changed.

Legacy Halo remains the sole executor and still constructs Vertex metadata
during synchronization. Moving this work to initialization would require
final topology and triplet registration, stable communicator/rank assignment,
preservation of rank-grouped all-to-all sender order and per-peer tags,
durable 3D boxes/transforms/lengths, planned buffer lifetimes, and rebuild
invalidation for every relevant topology or registration change. Coupling,
owner/alias reduce/broadcast, an ExchangeExecutor, and production metadata
migration remain out of scope. Parallel Unit 2B-3B is complete and frozen;
no production Halo executor migration was started.


## 39. Parallel Unit 2B-4A: initialization-time prepared Edge1Form Edge metadata acceptance

This unit adds an explicit, initialization-only preparation path for the
existing non-coupling Parallel Edge-stage Edge1Form protocol.  The public
entry point is `Halo::prepare_edge_1form_corner_metadata()`.  It is never
called by MercuryZ4 initialization or by a production time step; Z0 invokes it
explicitly as acceptance coverage.  It returns false when registered patterns
are not built.  The diagnostic accessors return value snapshots only while the
prepared cache is valid, including a generation and an accumulated metadata
exchange counter.

Preparation is eligible only for a complete, orientation-aware
`EdgeCovariant1Form` triplet at Corner2D or Corner3D level.  The final
eligibility decision also reads the active `FieldDescriptor` contract.  This
prevents a stale legacy triplet view from producing an oriented plan after
`Z0_FORM_ORIENTATION=0` has selected the ordinary component-copy path.  The
private cache owns copied legacy-shaped receive/send metadata; the only public
conversion is `LegacyOrientedCornerPlanAdapter::build_prepared`, which produces
the already frozen MPI-free `mpcns::parallel::OrientedCornerTransferPlan`.
It never calls prepare, performs MPI, reads field values, or reconstructs a
pattern.  No new public transfer type was introduced.

The prepared receiver sequence exactly follows filtered non-coupling legacy
`TOPO::EdgePatch` order, appending destination axes 0, 1, and 2 for each patch.
It preserves the direct Edge1Form source-axis rule
`send_axis = transform.perm[recv_axis]`, the exact `transform.sign[recv_axis]`,
the existing two-direction normal-first reorder rule, 2D receive and source
inner-ghost boxes, tags starting at zero independently for each peer, message
length, and local/per-peer order.  Its private metadata record has the same
field layout and byte-level sender-discovery semantics as legacy
`TripletCornerMeta`; it is not a production MPI buffer or a Field view.

Before any per-group metadata collective, all ranks build a deterministic
schedule sorted by group name and compare a descriptor containing the prepared
stage, generation state, group/field identities, levels, ghost depth,
orientation flag, locations, value kinds, and component counts.  The agreement
check uses `MPI_Allgather` followed by `MPI_Allgatherv`; a disagreement is
reported with rank and exits before the metadata collectives.  Every rank then
participates for every agreed group, including ranks with no local Parallel
Edge patch.  Per group the unchanged rank-grouped sender-discovery order is
reproduced with `MPI_Alltoall` counts and `MPI_Alltoallv` metadata.  This is
explicit initialization/test communication only: it neither changes the
runtime communicator/tag allocation nor invokes field pack, send, receive,
wait, unpack, or buffer allocation.

The cache is value-owned and immutable to callers.  It is invalidated on every
`register_halo_field`, `register_halo_fields`, and
`build_registered_patterns` call.  Registration increments a registry
generation and makes the pattern state unbuilt; successful pattern construction
increments a pattern generation.  A prepared entry is usable only when both
stored generations still match.  Direct external mutation of topology remains
outside Halo's observability boundary and must be followed by
`build_registered_patterns` before another explicit prepare.

`Z0_PreparedOrientedCornerPlanTest.cpp` independently derives expected
receiver data from raw `TOPO::EdgePatch` records, uses a test-only all-gather to
form the expected sender order, and compares every identity, axis, sign, box
bound, transform entry, peer, tag, size, local order, and peer order.  On a
mismatch it reports rank, group, Edge stage, role, peer, block, axis, tag,
local order, peer order, plan generation, field, expected, and actual values.
The numeric oracle still invokes the untouched legacy Edge-stage runtime sync
and checks distinct source values after component remapping and sign
application.  It verifies: absent-before-prepare rejection; first prepare;
byte-for-byte prepared/runtime snapshot parity; adapter conversion; repeated
prepare cache reuse without a new metadata exchange; rebuild invalidation and
regeneration; and atomic three-member triplet registration followed by the
required rebuild.  Ten pure contract checks retain invalid axis/transform/rank,
length/tag/role, and missing-send-box coverage.

Acceptance used exactly three MPI ranks and build parallelism no greater than
four; the two harnesses were run serially:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 compiled and linked `HaloPreparedMetadata.cpp` |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| `Z0_TEST_DRIVEN_SYNC=1 Z0_FORM_ORIENTATION=0 Z0_HALO_MODE=edgeforms Z0_HALO_LEVEL=edge mpirun -np 3 ./z0_null` from `build/Z0_NULL` | 0; component-copy path retained, no prepared plan or metadata exchange |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; explicit marker `validation_logs/unit2b4a_halo.exit` is `0`; all ten modes passed serially |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; explicit marker `preprocess_cube_validation_logs/unit2b4a_cube.exit` is `0`; all_vertex, all_owner_vertex, and owner_only passed serially |

The six active standard oriented Edge1Form mode instances produced 18
rank-local prepared plans, 792 receive records, 792 send records, and 792
matched pairs.  They made 770,034 exact metadata/lifecycle checks and 223,104
numeric value checks with zero failures; the four inactive modes added 180
negative-path checks.  Per active MPI invocation the lifecycle materialized
three group preparations per rank (the logged aggregate counter is 9), reused
the cache once per rank, rebuilt twice per rank, and invalidated once through
atomic triplet registration.  Observed topology coverage includes a rank with
zero local Parallel Edge patches, a rank with multiple peer ranks, and ranks
with multiple messages to the same peer.

The periodic cube's three modes produced nine rank-local plans, 180 receive
records, 180 send records, 180 matched pairs, 141,753 exact checks, and 16,896
numeric value checks, all with zero failures.  Cube transforms are identity
and positive; the standard topology supplies the non-identity permutations and
negative signs.  All logs retain the Z0 MPI-global pass summary and
`Z0_NULL finished normally`; MeshPartition, component-copy ExchangePlan,
Face-level oriented forms, all prior Edge/Vertex metadata snapshots, owner/
alias, DEC, Faraday, Hodge, and StateView regressions remain passing.

`Halo::exchange_parallel_edge_edge_1form_triplet_` remains the sole runtime
executor.  It still independently builds sender metadata during each sync and
continues to own all pack/send/recv/wait/unpack and buffers.  This unit adds no
runtime cache consumption, no ExchangeExecutor, no Edge/Vertex Face2Form
preparation, no coupling or owner/alias migration, and no production tag or
stage-order change.  The next smallest unit may prepare one remaining oriented
Corner family only after defining its own final-topology/registration barrier,
rank-global schedule agreement, immutable data lifetime, and invalidation
contract; it must not conflate that work with runtime executor migration.

Parallel Unit 2B-4A is complete and frozen.


## 40. Parallel Unit 2B-4B: initialization-time prepared Vertex Edge1Form metadata acceptance

This unit adds `Halo::prepare_vertex_edge_1form_metadata()` and a separate
value-owned cache for the non-coupling Parallel Vertex-stage Edge1Form
protocol.  It is an explicit initialization/test operation only: no
MercuryZ4 initialization or time-step path calls it, and legacy
HaloTriplet.cpp remains the sole runtime executor. debug_prepared_vertex_edge_1form_plan() and its generation/counter accessors expose copies only while the cache agrees with the current registration and pattern generations.
`LegacyOrientedVertexPlanAdapter::build_prepared()` is the only public
conversion boundary; it reads that immutable snapshot and constructs the
already-frozen `mpcns::parallel::OrientedVertexTransferPlan`.  It never
prepares, runs MPI, reads field values, or rebuilds a Halo pattern.

Eligibility is restricted to complete orientation-aware EdgeCovariant1Form
triplets registered at Corner3D and to non-coupling Parallel `VertexPatch`
records.  As with the legacy implementation, a receiver axis maps directly to
`transform.perm[recv_axis]` and uses `transform.sign[recv_axis]`.  The three
mapped directions retain the legacy normal-first correction, while receive
and source boxes retain the 3D ghost/inner-ghost rules.  Receiver patch order,
per-peer tag numbering, rank-grouped sender discovery, source/destination
block identities, message lengths, and local/per-peer orders are unchanged.
The schedule is rank-global and deterministic: each rank agrees on the sorted
triplet descriptors through `MPI_Allgather`/`MPI_Allgatherv`, then participates
in each agreed group’s unchanged `MPI_Alltoall`/`MPI_Alltoallv` metadata
exchange, including ranks with no local Parallel Vertex patch.

Vertex and Edge prepared caches are deliberately independent maps, readiness
states, counters, and snapshots.  They share only the registration/pattern
generations and the required invalidation boundary: every
`register_halo_field`, `register_halo_fields`, and
`build_registered_patterns` call invalidates both caches.  A prepare for one
stage does not evict the other.  Pattern rebuild or re-registration makes a
prepare fail until patterns are rebuilt.  This preserves the final-topology
and final-registration barrier without introducing a shared mutable plan or a
runtime buffer.

`Z0_PreparedOrientedVertexPlanTest.cpp` derives expected receive records from
raw legacy `TOPO::VertexPatch` data and derives sender order with a test-only
rank-grouped all-gather; it does not use the prepared adapter as an oracle.
It compares every group member, rank/block, direct axis map, sign, all three
mapped directions, box bound, transform component, tag, message length, and
local/peer order exactly.  Its numeric oracle invokes the untouched legacy
Vertex sync after the real Face/Edge prestages using distinct source values;
it proves the prepared snapshot and runtime capture have identical component
mapping and sign behavior.  Failures report rank, group, stage, role, peer,
block, axis, tag, region, and the exact field/expected/actual pair.

The test additionally verifies absent-plan adapter rejection, first prepare,
repeated cache reuse without another metadata exchange, Edge/Vertex cache
independence, rebuild invalidation and generation removal, complete-triplet
re-registration invalidation of both caches, rejection before rebuild, and
successful regeneration.  An `all_edge` regression showed why test activation
must follow the actual Halo registration rather than only a mode-limited Field
request: an earlier lifecycle test can legally register a Corner3D triplet
while normal sync coverage remains Edge-limited.  The acceptance test now
uses the registration-backed prepared snapshot as that decision point; no
production registration or stage behavior changed.

Acceptance used exactly three MPI ranks, build parallelism at most four, and
serial harness execution:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built and linked |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| focused `Z0_TEST_DRIVEN_SYNC=1 Z0_HALO_MODE=edgeforms Z0_HALO_LEVEL=vertex mpirun -np 3 ./z0_null` | 0; 3 prepared plans, 24 receive records, 24 send records, and 192 numeric checks |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; marker `validation_logs/unit2b4b_halo.exit` is `0`; all ten serial modes retain the MPI-global pass summary |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; marker `preprocess_cube_validation_logs/unit2b4b_cube.exit` is `0`; all_vertex, all_owner_vertex, and owner_only pass |

The focused orientation-disabled component-copy regression (`Z0_FORM_ORIENTATION=0`, `edgeforms`, `vertex`, three ranks) also exited 0: no prepared Vertex plan or metadata exchange was created.

Each active standard or periodic invocation reports three rank-local prepared
plans, 24 receive records, 24 send records, and 24 matched pairs with zero
metadata or numeric failures.  The standard topology exercises non-identity
permutations, negative Edge1Form signs, normal-first direction reordering,
and out-of-bounds source positions; the periodic cube supplies the
identity/positive periodic cases and retains MPI-global success.  Real Vertex
messages in these meshes have no zero-length packet; the existing pure
contract checks cover legal empty records.  Coupling vertices, Face2Form
preparation, ownership/alias paths, an ExchangeExecutor, and moving any
runtime pack/send/receive/wait/unpack or metadata exchange out of legacy Halo
remain deliberately out of scope.

Parallel Unit 2B-4B is complete and frozen.  A later unit may consider a
separate prepared Face2Form family only with its own exact legacy snapshot,
rank-global schedule, independent cache lifetime, and regression evidence;
it must not begin executor migration or merge Edge/Vertex metadata lifetimes.

## 41. Parallel Unit 2B-4C: initialization-time prepared Edge Face2Form metadata acceptance

This unit adds the explicit, initialization-only
`Halo::prepare_edge_face_2form_metadata()` path and a third, independent
value-owned prepared cache for the existing non-coupling Parallel Edge-stage
Face2Form protocol. The cache stores copied receive/send metadata only while
the current registration and pattern generations match. Its debug accessors
return no stale snapshot; `LegacyFace2FormEdgePlanAdapter::build_prepared()`
is the sole conversion boundary to the frozen
`mpcns::parallel::Face2FormEdgeTransferPlan`. The adapter neither starts
preparation nor performs MPI, accesses field values, rebuilds patterns, or
reconstructs topology pairing.

Eligibility requires a complete FaceXi/FaceEta/FaceZeta
`FaceContravariant2Form` triplet with `do_halo` and
`orientation_aware` true on all three final Field descriptors, a level of at
least Corner2D, and a non-coupling Parallel EdgePatch. The final catalog
check preserves component-copy behavior when `Z0_FORM_ORIENTATION=0`, even
though the legacy registry may normalize oriented form metadata internally.
With orientation disabled, the prepared Bface plan is absent and no Face2Form
prepared metadata exchange occurs.

Face2Form deliberately differs from Edge1Form: for each receiver axis,
`send_axis` is the inverse permutation satisfying
`transform.perm[send_axis] == recv_axis`, and the frozen sign is
`HALO_TOOLS::face_2form_orientation_sign(transform, send_axis)`
(the product of all transform signs times the source-axis sign). The
preparation preserves raw Parallel EdgePatch order, receiver axes 0/1/2,
normal-first two-direction reorder, receive ghost box, source inner-ghost
diagnostic box, transform, independently incremented per-peer tag,
message length, local order, and peer order. It does not reinterpret source
box ordering as a pack order.

All ranks first agree on a deterministic group schedule sorted by group name,
then each agreed group runs the legacy-shaped `MPI_Alltoall` count exchange
and rank-grouped `MPI_Alltoallv` sender discovery, including ranks without a
local Parallel EdgePatch. This explicit preparation has no production
time-step caller. `HaloTriplet.cpp` remains byte-for-byte unchanged and is
still the sole runtime executor, including its original per-sync metadata
collective and pack/send/receive/wait/unpack buffers.

The Face2Form cache is independent of prepared Edge1Form Edge and Vertex
caches. Preparing any one family does not evict another. Every
`register_halo_field`, `register_halo_fields`, and
`build_registered_patterns` invalidates all three, and a post-registration
prepare is rejected until patterns are rebuilt. External in-place topology
mutation remains caller-owned and requires an explicit rebuild.

`Z0_PreparedFace2FormEdgePlanTest.cpp` independently derives raw
TOPO::EdgePatch receive records and sender order, compares Raw/Prepared/
runtime snapshots field-for-field without sorting, and checks the public
adapter conversion. It covers first preparation, repeated cache reuse, two
legacy runtime sync captures, rebuild and registration invalidation, cache
independence, zero local Parallel patches, multi-block and multi-peer cases,
same-peer multiple messages, legal zero-length messages, negative signs,
out-of-bounds sources, non-finite source zeroing, and a non-self-inverse
permutation contract that rejects direct-axis mapping and the wrong 2-form
sign. Failures carry rank, group, Edge stage, role, peer, block, axis, tag,
order, field, expected, and actual values.

Acceptance used three MPI ranks, build parallelism at most four, and serial
harness execution:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| focused `Z0_TEST_DRIVEN_SYNC=1 Z0_HALO_MODE=faceforms Z0_HALO_LEVEL=edge mpirun -np 3 ./z0_null` from `build/Z0_NULL` | 0; three plans, 132 receive + 132 send records, 132 matched pairs, 87,296 exact/lifecycle checks, 9,184 finite mapped-value checks, 9,184 non-finite-source zero checks, 240 non-identity axes, 64 negative signs, and 16 zero-length messages |
| focused `Z0_TEST_DRIVEN_SYNC=1 Z0_FORM_ORIENTATION=0 Z0_HALO_MODE=faceforms Z0_HALO_LEVEL=edge mpirun -np 3 ./z0_null` from `build/Z0_NULL` | 0; no Face2Form prepared plan and zero prepared metadata exchanges |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; ten serial modes passed and all logs report the MPI-global normal conclusion |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; all_vertex, all_owner_vertex, and owner_only passed serially |

Across the six active standard modes, the new test produced 18 rank-local
plans, 792 receive and 792 send records, 792 matched pairs, 523,776 exact/
lifecycle checks, 55,104 finite mapped-value checks, and 55,104
non-finite-source zero checks, all with zero failures. Each active invocation
reported three cache hits, six rebuilds, three registration invalidations, and
nine preparation metadata exchanges across its three ranks. The periodic cube
added nine rank-local plans, 180 receive and 180 send records, 180 matched
pairs, 102,432 exact/lifecycle checks, 3,984 finite mapped-value checks, and
3,984 non-finite-source zero checks, also with zero failures. Its transforms
are identity and positive; the standard topology supplies non-identity,
negative-sign, zero-length, multi-block, and out-of-bounds coverage.

MeshPartition, component-copy ExchangePlan, both existing Edge1Form prepared
caches, all earlier Face2Form snapshots, DEC, Faraday, Hodge, StateView, and
owner/alias regression remain covered by the same global MPI conclusions.
No production tag, communicator, communication order, numerical formula,
coupling, or owner/alias path changed. Parallel Unit 2B-4C is complete and
frozen. The next smallest scoped unit is Vertex-stage Face2Form preparation,
which may reuse the schedule/cache architecture but must retain its distinct
three-direction metadata and must not migrate the runtime executor.

## 42. Parallel Unit 2B-4D: initialization-time prepared Vertex Face2Form metadata acceptance

Parallel Unit 2B-4D adds the fourth independent prepared-metadata cache:
`Halo::prepare_vertex_face_2form_metadata()`. It is explicit and
initialization/test-only. `debug_prepared_vertex_face_2form_plan()` returns a
value copy only while the registration and pattern generations still match;
the matching generation and metadata-exchange counter accessors expose no
stale state. `LegacyFace2FormVertexPlanAdapter::build_prepared()` is the
single read-only conversion boundary to the frozen
`mpcns::parallel::Face2FormVertexTransferPlan`. It neither initiates
preparation, calls MPI, accesses FieldBlock values, rebuilds patterns, nor
reconstructs topology.

Eligibility uses the final Field descriptors and requires a complete,
orientation-aware FaceXi/FaceEta/FaceZeta `FaceContravariant2Form` triplet
whose registered level includes Corner3D, plus a non-coupling Parallel
`VertexPatch`. This intentionally follows the actual Halo registry rather
than an initial mode-limited request. With `Z0_FORM_ORIENTATION=0`, Bface
remains on the existing component-copy path: no prepared Vertex Face2Form
plan or metadata collective is created.

The preparation is a literal read-only snapshot of the legacy Vertex Face2
protocol. For each receive axis it finds the inverse-permutation source axis
and uses `HALO_TOOLS::face_2form_orientation_sign(transform, send_axis)`.
It retains the three legacy directions with normal-first correction, 3D
face ghost and source inner-ghost boxes, patch traversal order, receiver and
sender local order, rank-grouped sender discovery, per-peer tags, block
addresses, transform, signs, message lengths, and component identities.
The rank-global sorted schedule is checked before the unchanged Alltoall/
Alltoallv metadata discovery; ranks with no local Parallel Vertex patch still
participate. This is preparation metadata only, not a new field-data path.

The Vertex Face2 cache has its own map, ready state, saved generations, and
counter. It is independent from prepared Edge1Form Edge, Edge1Form Vertex,
and Face2Form Edge caches. Every `register_halo_field`,
`register_halo_fields`, and `build_registered_patterns` invalidates all four
prepared families. The Z0 lifecycle test checks absence and generation
removal for all four after rebuild, rejects an adapter request for the stale
Vertex Face2 plan, and proves re-preparation does not evict the other three.

`Z0_PreparedFace2FormVertexPlanTest.cpp` derives raw receiver records from
`TOPO::VertexPatch` and sender order through a test-only rank-grouped
exchange. It compares raw, prepared, public-adapter, and two captured legacy
runtime snapshots in original order, field by field. The check includes the
three quantities and locations, rank/block identities, all transform
components, source/target axes, sign, three directions, 3D boxes, tag,
message length, and local/peer sequence. Distinct numeric component markers
then exercise the untouched legacy Face/Edge/Vertex runtime sequence,
including sign and component remapping. Diagnostics identify rank, group,
stage, role, peer, block, axis, tag, order, field, expected, and actual.

`HaloTriplet.cpp` is unchanged (SHA-256
`e03291c68741662289cd7fb9a04c105224e3f41522530c65e68662a8508feb30`), and no
production tag, send/receive/pack/unpack, wait order, buffer, owner/alias,
or numerical path changed.

Acceptance used at most four build jobs and exactly three MPI ranks; the two
harnesses were run serially:

| Command | Exit/result |
| --- | --- |
| `cmake --build build --parallel 4` | 0; MercuryZ4 built and linked |
| `cmake -S Z0_NULL -B build/Z0_NULL` | 0 |
| `cmake --build build/Z0_NULL --target z0_null --parallel 4` | 0 |
| focused `Z0_TEST_DRIVEN_SYNC=1 Z0_HALO_MODE=faceforms Z0_HALO_LEVEL=vertex mpirun -np 3 ./z0_null` from `build/Z0_NULL` | 0; 3 plans, 24 receive + 24 send records, 24 matched pairs, 13,485 exact/lifecycle checks, and 192 finite plus 192 non-finite mapped-value checks, all zero failures |
| focused orientation-disabled equivalent with `Z0_FORM_ORIENTATION=0` | 0; 0 plans, 0 metadata exchanges, 51 contract checks, zero failures |
| `NP=3 bash Z0_NULL/run_halo_validation.sh` | 0; all ten serial modes passed; every log reports the MPI-global normal conclusion |
| `NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh` | 0; all_vertex, all_owner_vertex, and owner_only passed serially |

The five active standard modes (faceforms_vertex, forms_vertex, all_vertex,
all_owner_vertex, owner_only) produced 15 rank-local plans, 120 receive and
120 send records, 120 matched pairs, 67,425 exact/lifecycle checks, 960
finite mapped-value checks, and 960 non-finite-source-zero checks: all
zero failures. Each active invocation recorded three cache reuses, six
rebuilds, three registry invalidations, and nine prepared metadata exchange
groups. The standard case provides 42 non-identity axis mappings, eight
negative signs, 18 normal-first direction reorders, three multi-local-block
ranks, and 36 out-of-bounds source positions per active run.

The three periodic-cube modes produced nine plans, 72 receive and 72 send
records, 72 matched pairs, 40,455 exact/lifecycle checks, and 576 each of
finite and non-finite-source-zero mapped-value checks, all zero failures.
Those generated periodic cases are identity/positive as expected and still
exercise the MPI-global conclusion and cache lifecycle.

Parallel Unit 2B-4D is complete and frozen. Coupling Vertex Face2Form
semantics, any new ExchangeExecutor, and every production Halo execution
migration remain out of scope. A later executor unit must consume these
frozen snapshots without changing legacy runtime ordering or metadata
collective semantics.
