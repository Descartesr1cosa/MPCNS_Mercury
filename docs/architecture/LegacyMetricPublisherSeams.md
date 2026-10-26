# Legacy Metric Publisher Seams

Status: Units 16-18 accepted; the three-rank Z0 halo suite passed with exact equality in all ten modes.
Verification: the pre-Unit-16 reference and new publishers matched every reconstruction-cache and diagonal-Hodge value exactly; the direct dr_* reference stencil also matched exactly (max=0.000000) in each halo mode. Unit 18 additionally verified that every legacy reconstruction and diagnostic facade has exactly the same result as its extracted compatibility implementation.

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

## Allowed next changes

1. Give a future `GeometryBuilder` a typed input and result model, then publish the result through these adapters.
2. Extract a pure 3x3 regularized pseudo-inverse only after its legacy equivalence test passes.
3. Give `SingularEdgeRegistry` topology/dual-geometry inputs and retain its current closure as the reference.
4. Keep non-diagonal/consistent Hodge work behind `HodgeOperator`; it must not modify this compatibility path.
