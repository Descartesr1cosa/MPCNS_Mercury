# Legacy Metric Publisher Seams

Status: Units 16-21 accepted; the three-rank Z0 halo suite passed with exact equality in all ten modes.
Verification: the pre-Unit-16 reference and new publishers matched every reconstruction-cache and diagonal-Hodge value exactly; the direct dr_* reference stencil also matched exactly (max=0.000000) in each halo mode. Unit 18 additionally verified that every legacy reconstruction and diagnostic facade has exactly the same result as its extracted compatibility implementation. Unit 19 verified the typed GeometryBuild readback and the pre-extraction singular-edge polygon/fallback closure reference in every mode. Unit 20 verified every reconstruction-cache field against the pre-extraction implementation after moving its pure kernels.

## Purpose

This document defines the temporary compatibility boundary that keeps the current `Grid` and `Field` backend stable while metric construction is separated into independently reviewable publishers. A publisher may read legacy storage and write the existing catalog, but it must not own MPI synchronization, a physics stencil, or a solver decision.

## Ordered production path

```text
Metric::compute_metric_fields
  1. populate_legacy_primal_face_and_cell_metrics  (primal face vectors and cell volume)
  2. populate_legacy_lumped_face_and_edge_metrics (face / edge / dual metrics and inverse star-1)
  3. populate_legacy_reconstruction_cache         (covariant bases and pseudo-inverse caches)
  4. populate_legacy_primal_edge_vectors            (three primal edge-vector fields)
  5. populate_legacy_diagonal_hodge_closure        (star 0<->3 and forward star 1->2)
```

The order is contractual. In particular, step 1 publishes JDxi, JDet, JDze, and Jac; step 3 consumes `JDxi`, `JDet`, and `JDze`; step 5 consumes `Jac` and the inverse edge Hodge scales published in step 2. Step 4 is intentionally kept after reconstruction to preserve historical execution order. No caller may reorder these steps or repeat a subset without a dedicated regression case.

## Publisher contracts

| Publisher | Reads | Writes | Dependency rule |
| --- | --- | --- | --- |
| populate_legacy_primal_face_and_cell_metrics | node coordinates | JDxi, JDet, JDze, Jac | Compatibility layer; delegates pure quadrilateral and cell-volume kernels |
| populate_legacy_primal_edge_vectors | node coordinates | dr_xi, dr_eta, dr_zeta | Compatibility layer only; no Hodge, topology, or MPI |
| populate_legacy_lumped_face_and_edge_metrics | node coordinates, primal face vectors, Jac | dl_*, dlstar_*, Sstar_*, Astar_*, inverse star-1 | Compatibility layer; delegates pure areas/norms and pure axis policy |
| populate_legacy_reconstruction_cache | node coordinates, JDxi, JDet, JDze | a_xi, a_eta, a_zeta, pinvGT_cell, pinvGT_*, pinvAT_* | Compatibility layer only; no Hodge field and no MPI |
| populate_legacy_diagonal_hodge_closure | Jac, inverse star-1 | star 0->3, star 3->0, star 1->2 | Compatibility layer only; no coordinates, topology, or MPI |

The primal publisher preserves the historical field loops and calls only the already-isolated pure geometry primitives. The reconstruction publisher retains the legacy Tikhonov regularization and retry sequence in the same arithmetic order. The diagonal closure retains its zero guards (`1e-300`) and has no axis or singular-edge decision logic.


## Field-backed consumers

`Metric.cpp` is now a coordinator and compatibility facade. The Field-backed
consumer algorithms that used to be interleaved with construction are owned by
two explicit adapters:

| Adapter | Owns | Public compatibility path |
| --- | --- | --- |
| `legacy_field_reconstruction` | face 2-form and edge 1-form reconstruction at cell and node samples | `METRIC::reconstruct_*` forwards directly |
| `legacy_metric_diagnostics` | finite/positive metric scans and axis counters | `METRIC::diagnose_metric_fields` maps the result structure directly |

Neither adapter performs MPI, creates geometric fields, selects a Hodge
model, or makes a solver decision. The Z0 suite compares each `METRIC`
facade call and every diagnostic counter against the extracted adapter at an
inner sample of every local block; this is a delegation gate, not an
independent numerical oracle.

## Typed build hand-off

`LegacyGeometryBuildPipeline` now owns the compatibility sequence: register the 43 legacy fields, publish the five ordered metric groups, borrow the resulting field views, and pass an ordered `GeometryBuildInput` to `GeometryBuilder::finalize`. The result is a structurally complete, immutable-in-use `GeometrySnapshot`; it still borrows the current storage and makes no numerical or MPI decision.

The variable-valence edge registry now supplies `SingularEdgeGeometryInput` to the pure geometry closure. It retains alias equivalence, field collection/publication, ownership, and MPI reductions. The closure module owns only center ordering, polygon area, the existing volume-over-length fallback, and the corresponding inverse Hodge value.

## Pure reconstruction kernels and singular publication

`reconstruction_kernels` now owns the four-edge cell covariant-basis average and the regularized 3x3 pseudo-inverse. `legacy_reconstruction_metric_publisher` owns only legacy coordinate/Field access and edge-local stencil assembly. The Z0 pre-Unit-16 reference compares all `a_*`, `pinvGT_*`, and `pinvAT_*` field entries exactly in every mode.

`legacy_singular_edge_metric_publisher` owns only the final `Astar_*` and inverse-Hodge catalog writes for a closed singular edge. `SingularEdgeRegistry` retains topology aliases, MPI reductions, geometry input collection, closure, and weights. No owner or synchronization schedule changed.

## Allowed next changes

1. Replace the borrowed GeometryBuild result with owned immutable storage only after field-by-field Mercury comparison is available.
2. Replace the borrowed GeometryBuild result with owned immutable storage only after field-by-field Mercury comparison is available.
3. Replace the singular-edge collector Field/Grid source with typed geometry and partition views only after a protected Mercury case-level comparison, while retaining its MPI reduction schedule.
4. Keep non-diagonal/consistent Hodge work behind `HodgeOperator`; it must not modify this compatibility path.

## Unit 21: singular-edge collection boundary

`mpcns::mesh` now owns the small, storage-free record vocabulary for quotient-edge aliases and weighted incident entities. `mpcns::compat::LegacySingularEdgeDataSource` is the sole Grid/Field boundary for singular-edge construction, and `collect_legacy_singular_edge_measures` owns the historical reads of `dl_*`, `dr_*`, `Jac`, `Area_*`, and dual-cell centers.

`SingularEdgeRegistry::build` now consumes topology plus that explicit data source. It retains only candidate selection, quotient ownership, the existing MPI all-reduce/all-gather schedule, physical-face valence classification, pure geometry closure, weight normalization, and the existing singular-Hodge publisher. Thus the registry public header no longer exposes Grid, and the legacy data reads have one reviewable replacement point. No formula, alias orientation rule, MPI collective order, or solver behavior changed.

Verification for this unit: Z0 and MercuryZ4 both compiled with at most four build jobs; the three-rank Z0 suite completed all ten halo modes successfully.
