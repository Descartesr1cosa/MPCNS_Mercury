# Topology Module

This directory contains the public C++ interface of the topology module.

The topology module defines the metric-free discrete structure of the multi-block mesh.

It provides:

- entity keys for node, edge, face, and cell;
- local orientation convention;
- block connectivity;
- index transformation across interfaces;
- entity equivalence classes;
- owner / alias / orientation sign relation;
- incidence stencils d0, d1, d2.

The topology module does not compute:

- geometry;
- Jacobian;
- edge length;
- face area;
- cell volume;
- Hodge star;
- curl, div, grad as metric-scaled numerical operators;
- halo communication;
- physical source terms;
- PETSc solvers.

---

## 1. Directory role

This directory is the C++ public interface of the topology layer.

It should be included by:

- halo;
- metric;
- DEC;
- solver-side global DOF builders;
- topology tests.

It should not depend on:

- physics;
- FV;
- solver;
- PETSc;
- MPI communication implementation;
- Mercury-specific application code.

---

## 2. Main files

| File | Responsibility |
|---|---|
| `TopologyTypes.h` | basic enums and structs |
| `EntityIndex.h` | helper functions for creating node/edge/face/cell keys |
| `EntityOrientation.h` | orientation utility functions, if needed |
| `BlockConnectivity.h` | block-to-block connection graph |
| `TopologyEquiv.h` | equivalence class, owner, alias, and sign relation |
| `Incidence.h` | metric-free d0, d1, d2 incidence stencils |
| `Topology.h` | public facade of the topology module |

External modules should preferably include only:

```cpp
#include "2_topology/Topology.h"
```

Internal implementation files may include more specific headers.

---

## 3. Conceptual responsibilities

### 3.1 Entity identity

Topology defines how to identify a discrete object:

```text
Node  : block, i, j, k
Edge  : block, i, j, k, axis
Face  : block, i, j, k, axis
Cell  : block, i, j, k
```

The `axis` field distinguishes:

```text
EdgeXi, EdgeEta, EdgeZeta
FaceXi, FaceEta, FaceZeta
```

---

### 3.2 Orientation

Topology defines the reference orientation of edges, faces, and cells.

Example:

```text
EdgeXi(i,j,k) goes from Node(i,j,k) to Node(i+1,j,k)
```

Face and cell orientations are defined in `docs/topology/TopologySpec.md`.

---

### 3.3 Incidence

Topology provides metric-free incidence stencils:

```text
d0: node -> edge
d1: edge -> face
d2: face -> cell
```

These stencils only contain integer signs.

They do not contain:

```text
length
area
volume
Jacobian
Hodge star
```

---

### 3.4 Equivalence class

Topology defines when several local entities represent the same topological entity.

Examples:

```text
same interface node stored by two blocks
same interface edge stored by two blocks
same interface face stored by two blocks
repeated axis node
periodic duplicate entity
```

Each equivalence class has exactly one owner.

---

### 3.5 Orientation sign

Topology provides the sign between an alias and its owner.

The sign is used by orientation-aware fields such as:

```text
edge 1-form
face 2-form
```

The sign is ignored by non-orientation-aware fields such as:

```text
cell scalar
cell conserved vector
Cartesian vector field
```

The decision to apply the sign belongs to halo and field descriptor logic.

---

## 4. Dependency rule

Allowed dependencies:

```text
topology -> standard library
```

Possible future dependency:

```text
topology -> basic utility module
```

Forbidden dependencies:

```text
topology -> grid metric implementation
topology -> halo communication implementation
topology -> physics
topology -> solver
topology -> PETSc
topology -> MPI send/recv
topology -> Mercury application
```

Topology may store MPI rank id as metadata if needed for deterministic owner selection, but it must not perform MPI communication.

---

## 5. First-stage implementation target

The first implementation stage should provide:

1. `EntityKey`;
2. `EntityDim`;
3. `EdgeAxis`;
4. `FaceAxis`;
5. entity construction helpers;
6. local incidence stencils;
7. minimal `TopologyEquiv`;
8. minimal `BlockConnectivity`;
9. minimal `Topology` facade;
10. tests for `d1*d0 = 0` and `d2*d1 = 0`.

This first stage should not be connected to the Mercury production solver yet.

---

## 6. Testing requirements

The topology module must eventually pass:

| Test | Purpose |
|---|---|
| single-block incidence | verify local orientation |
| two-block interface | verify owner/alias relation |
| reversed orientation interface | verify sign relation |
| periodic interface | verify periodic equivalence |
| rotational axis | verify collapsed/repeated entity handling |

The most fundamental identities are:

```text
d1 * d0 = 0
d2 * d1 = 0
```

If these fail, the topology orientation convention is inconsistent.

---

## 7. Relation to other modules

| Question | Responsible module |
|---|---|
| Which block touches which block? | topology |
| How does local index map across an interface? | topology |
| Is an edge reversed across an interface? | topology |
| Should halo multiply by a sign for this field? | halo + field descriptor |
| What is the physical edge length? | grid / metric |
| What is the Hodge star? | metric |
| How to compute curl E? | DEC |
| How to compute Hall EMF? | physics |
| How to solve implicit Hall equation? | solver |

---

## 8. Recommended usage

Most external modules should depend on the public facade:

```cpp
#include "2_topology/Topology.h"
```

Specialized topology tests may include:

```cpp
#include "2_topology/Incidence.h"
#include "2_topology/TopologyEquiv.h"
#include "2_topology/BlockConnectivity.h"
```

Avoid including internal implementation headers outside topology and tests.

---

## 9. Current status

This module is under refactor.

The current goal is to establish a clean topology contract before upgrading halo, metric, DEC, and physics.
