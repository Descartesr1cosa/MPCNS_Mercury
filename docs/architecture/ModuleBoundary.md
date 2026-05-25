# MPCNS Mercury Module Boundary

Status: Draft  
Owner: Ke Xu  

## 1. Purpose

This document defines the responsibility boundary of the refactored MPCNS Mercury code.

The main design principle is:

- `grid` handles geometry;
- `topology` handles metric-free connectivity and oriented complexes;
- `field` handles data storage and mathematical field identity;
- `halo` handles synchronization across blocks, MPI ranks, and owner-alias entities;
- `metric` handles Hodge stars, inner products, and geometry-induced reconstruction;
- `dec` handles exterior derivative, codifferential, and structure-preserving diagnostics;
- `fv` handles cell-centered compressible finite-volume discretization;
- `physics` combines generic operators into concrete physical models;
- `solver` handles time integration and nonlinear/linear solver flow;
- `app` handles concrete case configuration.

The purpose of this document is to prevent geometry, topology, synchronization, DEC operators, physics, and solver logic from being mixed together.

---

## 2. Module responsibility table

| Module | Responsibility | Must not do |
|---|---|---|
| `include/1_grid` | coordinate mapping, Jacobian, volume, area, edge length, basis vectors | block connectivity, incidence, curl/div, physics |
| `include/2_topology` | block connection, index transform, entity equivalence, owner/alias/sign, incidence d0/d1/d2 | metric, Hodge star, halo communication, physical operators |
| `include/3_field` | field descriptor, field storage, mathematical identity of fields | geometry, topology construction, physical operators |
| `include/4_halo` | ghost/interface/owner-alias synchronization | curl/div/Hall/Lorentz force |
| `include/5_metric` | Hodge star, inner product, vector-form projection, reconstruction | topology incidence, MPI communication |
| `include/6_dec` | d0, d1, d2, codifferential, structural diagnostics | physical model decisions |
| `include/7_fv` | reconstruction, limiter, Riemann flux, conservative update | DEC topology, Hodge star |
| `include/8_physics` | ideal MHD, Hall EMF, resistive induction, multi-ion source terms, planetary induction model | reimplementing curl/div/Hodge/halo |
| `include/9_solver` | time stepping, operator splitting, IMEX, PETSc SNES/KSP, residual assembly | geometry/topology ownership |
| `apps/mercury` | mesh choice, physical parameters, initial condition, boundary model, output | generic numerical algorithms |

---

## 3. Golden rule

Topology tells how discrete objects are connected.

Geometry and metric tell how these objects are measured.

DEC combines topology and metric into structure-preserving differential operators.

Physics only combines these generic operators into concrete physical models.

Solver only controls time advancement and nonlinear/linear solution procedures.

Application code only configures a specific physical case.

---

## 4. Layer dependency rule

The intended dependency direction is:

```text
app
  -> solver
      -> physics
          -> fv
          -> dec
              -> topology
              -> metric
                  -> grid
          -> halo
              -> topology
              -> field
          -> field
```

Lower-level modules should not depend on higher-level modules.

For example:

- `topology` must not depend on `physics`;
- `metric` must not depend on `solver`;
- `halo` must not depend on Hall-MHD details;
- `field` must not know how to compute curl or flux;
- `grid` must not know how blocks are connected.

---

## 5. What belongs to topology?

The topology module is responsible for metric-free discrete structure.

It owns:

- block-to-block connectivity;
- local-to-remote index transform;
- node/edge/face/cell entity keys;
- entity orientation convention;
- topological equivalence classes;
- owner and alias relations;
- orientation signs;
- incidence relations d0, d1, d2.

It does not own:

- length;
- area;
- volume;
- Jacobian;
- metric tensor;
- Hodge star;
- numerical curl/div/grad;
- MPI send/recv;
- Hall EMF;
- Lorentz force;
- PETSc vectors or matrices.

---

## 6. What belongs to grid?

The grid module is responsible for geometry.

It owns:

- coordinate mapping;
- physical coordinates;
- Jacobian;
- covariant basis vectors;
- contravariant basis vectors;
- cell volume;
- face area;
- edge length;
- local geometric quality information.

It does not own:

- block connectivity;
- entity equivalence class;
- owner/alias relation;
- topological incidence;
- halo communication;
- physical source terms.

---

## 7. What belongs to field?

The field module is responsible for data and mathematical field identity.

It owns:

- field name;
- field location;
- number of components;
- value kind;
- form degree;
- orientation awareness;
- owner-sync requirement;
- storage layout.

Example field identities:

| Field | Location | Mathematical identity | Orientation-aware |
|---|---|---|---|
| density | cell | scalar/conserved variable | no |
| momentum | cell | Cartesian vector/conserved variable | no |
| magnetic flux B | face | 2-form | yes |
| electric EMF E | edge | 1-form | yes |
| potential phi | node | 0-form | no |
| current J | edge | 1-form or edge-collocated vector, depending on formulation | depends on descriptor |

The field module must not implement physical or differential operators.

---

## 8. What belongs to halo?

The halo module is responsible for synchronization.

It owns:

- cell ghost synchronization;
- face ghost synchronization;
- edge ghost synchronization;
- block interface copy;
- MPI send/recv;
- owner-to-alias synchronization;
- alias-to-owner reduction if needed;
- applying topology orientation sign when the field descriptor requires it.

The halo module uses:

- topology relation;
- field descriptor;
- communication pattern.

The halo module must not compute:

- curl;
- divergence;
- Hodge star;
- Hall EMF;
- Lorentz force;
- fluid flux.

---

## 9. What belongs to metric?

The metric module converts geometry into discrete metric structures.

It owns:

- Hodge star operators;
- primal/dual volume ratios;
- edge/face/cell inner products;
- vector-to-form projection;
- form-to-vector reconstruction;
- B_face to B_cell reconstruction;
- J_edge to J_cell reconstruction;
- near-axis metric regularization;
- degenerate geometry treatment at the metric level.

The metric module uses geometry from `grid` and may use topological entity information from `topology`.

It must not define incidence signs or owner-alias relations.

---

## 10. What belongs to DEC?

The DEC module provides structure-preserving discrete differential operators.

It owns:

- d0: node 0-form to edge 1-form;
- d1: edge 1-form to face 2-form;
- d2: face 2-form to cell 3-form;
- codifferential operators;
- discrete curl;
- discrete divergence;
- discrete gradient skeleton;
- div-curl diagnostic;
- curl-grad diagnostic;
- energy inner product diagnostic.

DEC combines:

- incidence from topology;
- Hodge star and inner products from metric.

DEC must not implement Hall-MHD or fluid Riemann solvers.

---

## 11. What belongs to FV?

The finite-volume module handles cell-centered compressible flow discretization.

It owns:

- reconstruction;
- limiter;
- Riemann solver;
- numerical flux;
- conservative update;
- cell-centered source term infrastructure.

It must not own:

- DEC incidence;
- Hodge star;
- owner-alias synchronization;
- block topology.

---

## 12. What belongs to physics?

The physics module combines generic numerical operators into physical models.

It owns:

- ideal MHD EMF construction;
- Hall EMF construction;
- resistive electric field;
- Lorentz force;
- multi-ion source terms;
- electron-pressure closure;
- planetary induction model;
- Mercury-specific physical closures.

It must call generic operators from DEC, metric, FV, field, and halo instead of reimplementing them.

---

## 13. What belongs to solver?

The solver module controls solution procedures.

It owns:

- explicit time stepping;
- implicit time stepping;
- operator splitting;
- IMEX schemes;
- nonlinear residual assembly;
- PETSc SNES/KSP setup;
- convergence control;
- timestep selection.

The solver may use topology owner information to build global unknown numbering, but topology itself must not depend on PETSc.

---

## 14. What belongs to app?

Application code configures a concrete simulation.

For `apps/mercury`, it owns:

- mesh selection;
- Mercury physical parameters;
- solar wind parameters;
- planetary magnetic field model;
- boundary condition selection;
- initial condition;
- output variables;
- restart configuration.

It must not implement generic topology, halo, DEC, metric, or solver algorithms.

---

## 15. Refactor rule

When adding new code, ask:

1. Is this about measuring geometry?  
   Then it belongs to `grid` or `metric`.

2. Is this about which discrete object touches which object?  
   Then it belongs to `topology`.

3. Is this about how data is stored and identified?  
   Then it belongs to `field`.

4. Is this about copying/synchronizing data?  
   Then it belongs to `halo`.

5. Is this about d, curl, div, codifferential, or structural identities?  
   Then it belongs to `dec`.

6. Is this about finite-volume fluid flux?  
   Then it belongs to `fv`.

7. Is this about MHD/Hall/multi-ion/Mercury physics?  
   Then it belongs to `physics`.

8. Is this about timestep or nonlinear/linear solve?  
   Then it belongs to `solver`.

9. Is this about choosing one specific Mercury simulation setup?  
   Then it belongs to `apps/mercury`.

---

## 16. Current priority

The first refactor target is `include/2_topology`.

The initial topology refactor should establish:

1. entity key system;
2. local orientation convention;
3. incidence stencils d0, d1, d2;
4. block connectivity;
5. entity equivalence classes;
6. owner/alias/sign relation;
7. tests for d1*d0 = 0 and d2*d1 = 0;
8. tests for interface owner/sign consistency.
