# Mercury Geometry Inventory and Migration Boundary

Status: audited; typed GeometryBuilder hand-off, pure reconstruction kernels, borrowed GeometrySnapshot compatibility boundary, legacy field-publishing contract, and Metric consumer facades implemented

## 1. The two existing geometry layers

The current implementation has two valid, but entangled, representations of
the same mesh geometry.

| Layer | Owner today | Contents | Correct responsibility after migration |
|---|---|---|---|
| Coordinate and grid-preprocessing layer | `Grid::Block` | node coordinates `x/y/z`, legacy Jacobian and metric arrays, dual-cell centers `dual_x/y/z`, GCL face metrics and tangent bases | `NodeCoordinates`, `PrimalGeometry`, and reconstruction input |
| Field-derived DEC layer | `METRIC::compute_metric_fields` | face-area vectors, primal/dual measures, diagonal Hodge scales, and reconstruction pseudo-inverses | `PrimalGeometry`, `DualGeometry`, `LumpedHodge`, and `ReconstructionCache` |
| Quotient-edge closure | `METRIC::SingularEdgeRegistry` | pole/extraordinary-edge aliases, canonical edge, variable-valence dual area, and corrected inverse Hodge value | topology equivalence plus `DualGeometry` closure |

`Grid::Grid_Preprocess()` owns the one-time coordinate exchange and grid
preparation.  `METRIC::compute_metric_fields()` subsequently reads that closed
coordinate view and writes the Field-backed derived data.  It must remain free
of MPI.  `SingularEdgeRegistry` is a separate topological closure step, not a
generic Hodge stencil.

## 2. Canonical meaning of the production fields

The compatibility adapter `mpcns::compat::LegacyGeometryFieldIds` is now the
only place where the following legacy catalog strings are resolved for the
Mercury solver.

| Canonical identity | Legacy field family | Location |
|---|---|---|
| `primal_cell_volume` | `Jac` | cell |
| `primal_face_area_vector` | `JDxi/JDet/JDze` | primal face |
| `primal_edge_vector` | `dr_xi/dr_eta/dr_zeta` | primal edge |
| `cell_covariant_basis` | `a_xi/a_eta/a_zeta` | cell |
| `primal_face_area` | `Area_*` | primal face |
| `dual_edge_length_at_face` | `dlstar_*` | primal face / dual edge |
| `primal_edge_length` | `dl_*` | primal edge |
| `dual_face_area_vector_at_edge` | `Sstar_*` | primal edge / dual face |
| `dual_face_area_at_edge` | `Astar_*` | primal edge / dual face |
| `lumped_star_primal_2_to_dual_1` | `Hodge_star_2form_to_1form_face_*_lumped` | primal face |
| `lumped_inverse_star_primal_1` | `Hodge_star_inverse_2form_to_1form_edge_*_lumped` | primal edge |
| `cell_reconstruction_gradient_pseudoinverse` | `pinvGT_cell` | cell |
| `edge_reconstruction_gradient_pseudoinverse` | `pinvGT_*` | primal edge |
| `edge_reconstruction_area_pseudoinverse` | `pinvAT_*` | primal edge |

The compatibility aliases `Face_beta` and `Edge_alpha` remain inside
`SolverFields` only so existing physics is unchanged.  New code must use the
canonical Hodge meaning rather than either alias.

## 3. Dependency boundary now established

~~~text
legacy Grid + Field catalog
        |
        v
mpcns::compat::LegacyGeometryFieldIds  <-- only legacy field-name mapping
        |
        v
Mercury SolverFields                   <-- existing numerical storage unchanged
~~~

`mpcns::state::PrimalEdgeFormView` and `PrimalFaceFormView` now belong to the
state layer.  Incidence code consumes those typed cochains and no longer owns
their storage definition.  `mpcns_discrete` and `mpcns_legacy_geometry` are
explicit CMake targets; the former is not compiled through an application
source glob.  The compatibility layer can also bind a block's exact existing
Hodge fields as `LumpedHodgeCoefficients`; no copied geometric state or new
formula is introduced.

`mpcns::compat::legacy_geometry_field_names` is the single catalog-spelling
directory for this migration.  `METRIC::register_metric_fields` and
`LegacyGeometryFieldIds` use the same constants, so a future catalog rename
cannot silently register one name and bind another.

The registration metadata itself is now a separate compatibility contract:
legacy_geometry_field_registration owns all 43 legacy entries and preserves
their historical registration order, stagger location, component count, and
ghost-width rule.  It is the publishing adapter from a future geometry
builder to the current Field backend; it performs no metric calculation.
Z0 verifies that the live catalog exactly matches this contract.

## 4. GeometrySnapshot compatibility boundary

`mpcns::geometry::GeometrySnapshot` is now a read-only, type-erased view of
all local geometry blocks.  The geometry target itself includes neither legacy
`Grid`/`Field`, MPI, nor Mercury headers.  Its vocabulary is deliberately
small and stable:

~~~text
NodeCoordinates
PrimalGeometry: cell volume, face-area vector, edge vector, face area, edge length
DualGeometry: dual-edge length, dual-face-area vector, dual-face area
ReconstructionCache: covariant bases and cell/edge pseudo-inverses
GeometrySnapshot: ordered local GeometryBlock views
~~~

The first calculation extracted from the legacy metric routine is the pure
geometry primitive geometry_primitives::quadrilateral_area_vector.  It has no
Grid, Field, MPI, or Mercury dependency and preserves the historical
two-triangle orientation and operation order.  The legacy metric publisher
uses it only for the three primal face-area-vector families.  A Z0 regression
recomputes all local JDxi, JDet, and JDze entries using the former stencil;
only an exact comparison permits further primal-geometry movement.

The cell-volume quadrature is now isolated as
geometry_primitives::cell_volume_from_oriented_faces.  Its input is six
already-oriented face area vectors and centroids plus the historical dual-cell
reference point.  The primitive preserves the original ordered accumulation;
the compatibility metric code retains all Grid/Field reads.  Z0 compares every
local Jac entry with an independent copy of the legacy six-face stencil before
the migration can continue.

Dual-face area vectors at primal edges (`Sstar_xi`, `Sstar_eta`, and
`Sstar_zeta`) now use the same pure quadrilateral-area primitive. Their
topology-dependent choice of the four surrounding dual-cell centers remains
in the legacy publisher, together with the field writes and inverse-Hodge
calculation. The Z0 dual-face gate has passed: it independently evaluated the
former two-triangle stencil for every vector component and local ghost entry.

geometry_primitives::euclidean_norm now also supplies primal-edge lengths and
dual-face magnitudes; its Z0 norm-stencil gate has passed. Axis-degeneracy
treatment is in the pure discrete-layer AxisSafeInverseStar1Policy, with
constructed-branch and field-derived Z0 checks passed.

The legacy metric writes are now organized as four publisher responsibilities:
primal face/cell measures (JDxi, JDet, JDze, Jac), face/edge measures plus
inverse star-1, reconstruction caches, and diagonal Hodge closure. The primal
edge vectors remain a separate ordered publisher because their historical
execution location is retained intentionally. Unit 17 moved the primal two
responsibilities from Metric.cpp without changing their formulas. Existing
primal-face and cell-volume gates plus the new full-range dr_* stencil gate
passed with exact equality in all ten three-rank Z0 halo modes.

The remaining Field-backed metric consumers no longer share the publisher
implementation file. `legacy_field_reconstruction` owns the four historical
face/edge form reconstructions at cell/node samples, while
`legacy_metric_diagnostics` owns the metric-health scans. `METRIC` preserves
the former public functions as direct facades, and the Z0 suite checks their
results and all diagnostic counters against the extracted modules in every
halo mode. This is an interface-preservation gate; it introduces no new
geometry storage or formula.

Unit 19 adds the first complete typed build hand-off. `LegacyGeometryBuildPipeline` owns legacy registration and the contractual five publisher calls; it then converts the borrowed views into `GeometryBuildInput` and asks `GeometryBuilder` to return a validated `GeometryBuildResult`. The result is still a zero-copy view, but no geometry consumer now needs to know the legacy construction path. Z0 checks that the pipeline result reads exactly the same local coordinates, metrics, and reconstruction cache as the Field backend.

The extraordinary-edge dual polygon is also now a pure `mpcns::geometry` operation. `SingularEdgeRegistry` gathers aliases, centers, lengths, and fallback measures through the existing topology/MPI route, while `close_singular_edge_geometry` computes center ordering, polygon area, fallback selection, inverse Hodge, and canonical vector. A test-only copy of the pre-extraction polygon/fallback arithmetic agrees exactly for nondegenerate, insufficient, and collinear center sets.

Unit 20 moves the remaining pure reconstruction arithmetic out of the Field publisher: `cell_covariant_basis_from_corners` preserves the historical four-edge averages, while `regularized_pseudoinverse_transpose` preserves the exact Gram product, regularization, retry, and zero-fallback sequence. The publisher now owns only legacy stencil gathering and catalog writes. The pre-extraction Z0 reference still compares every reconstruction-cache field exactly. The final singular-edge `Astar_*` and inverse-Hodge writes are likewise isolated in `legacy_singular_edge_metric_publisher`; registry ownership, alias orientation, and MPI behavior remain unchanged.


`mpcns::compat::borrow_legacy_geometry_snapshot(Grid&, Field&)` is the only
legacy-aware constructor.  It borrows the current coordinates and derived
metric fields exactly; it allocates nothing and changes no numerical result.
Its lifetime is explicitly tied to the supplied legacy `Grid` and `Field`.
The compatibility snapshot therefore establishes the target interface before
the storage migration, rather than copying a second, potentially divergent
geometry cache.

Z0 now checks every snapshot category against the original storage at a local
sample point, including all directions and vector components.  This is exact
readback validation, not a tolerance-based geometry comparison.  The complete
ten-mode, three-rank Z0 halo suite passed this gate with `max=0.000000`.

## 5. Next safe geometry migrations

1. Replace the borrowed compatibility backing with owned immutable geometry
   storage only after a full field-by-field comparison against each protected
   Mercury case is available.
2. Replace the borrowed `GeometryBuildResult` with owned `PrimalGeometry`,
   `DualGeometry`, `LumpedHodge`, and `ReconstructionCache` storage while
   retaining the old field catalog as a publishing adapter.
3. Replace the completed singular-edge collector boundary with typed geometry and partition views only after a protected Mercury case-level comparison; retain its MPI schedule and current closure as the reference.

4. Replace `Face_beta` and `Edge_alpha` callers with a `LumpedHodge` API only
   after field-by-field agreement diagnostics exist.

No consistent/non-diagonal Hodge, MPI schedule, boundary action, or solver
formula is changed by the current adapter work.

## 6. Regression gate before production migration

Z0 registers a local dual-edge workspace with the same ghost reach as the
current metric fields.  Its DEC-chain test computes electric current through
`LumpedHodge::apply_star_2`, `d1^T`, and
`LumpedHodge::apply_inverse_star_1`, then compares every local inner edge
value to the former fused CT formula.  It uses no workspace halo exchange.
Only an exact zero difference permits migration of Mercury's current and
implicit-resistive paths to the new composition.

## 7. Unit 21: singular-edge collection split

The variable-valence registry no longer reads legacy metric fields or Grid coordinates directly. `legacy_singular_edge_collector` is the compatibility adapter for local lengths, vectors, cell/face measures, and dual centers; the registry consumes its local result and remains responsible for the unmodified topology/MPI/closure sequence. This leaves a precise future replacement seam for typed GeometrySnapshot and partition views, without exposing those storage choices to physics callers.

Z0 and MercuryZ4 compiled with no more than four build jobs, and the complete three-rank Z0 halo suite passed in all ten modes.
