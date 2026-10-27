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
