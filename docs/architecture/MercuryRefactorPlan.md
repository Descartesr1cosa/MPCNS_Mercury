# MPCNS Mercury Refactor Plan

Status: Accepted design baseline; Units 1-11 implementation in review (3-rank Z0 Hodge, GeometrySnapshot, field-name, field-publishing, quadrilateral-area, and cell-volume gates passed; face-metric gate awaiting validation)
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

geometry_primitives::euclidean_norm is now used only by the three primal-face
metric paths: face area, dual-edge length, and the diagonal star-2 scale.
The existing edge/dual-face and inverse-Hodge calculations still use their
legacy local implementation.  Z0 independently checks every local area,
dual-edge length, and star-2 value against the former stencil before this
boundary expands.

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
