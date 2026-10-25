# Mercury Geometry Inventory and Migration Boundary

Status: audited; borrowed GeometrySnapshot compatibility boundary implemented

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
2. Move the field-derived metric construction behind a builder that produces
   `PrimalGeometry`, `DualGeometry`, `LumpedHodge`, and
   `ReconstructionCache` separately while retaining the old field catalog as
   a publishing adapter.
3. Give `SingularEdgeRegistry` a topology/dual-geometry input instead of
   `Field` and `Grid`; retain its current numerical closure as the reference.
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
