# Topology Module Specification

Status: Draft  
Owner: Ke Xu  

## 1. Purpose

The topology module defines the metric-free discrete complex of a multi-block structured mesh.

It answers the following questions:

1. Which local node, edge, face, or cell belongs to which topological entity?
2. How are entities connected across block interfaces?
3. Which entity is the owner of an equivalence class?
4. What is the orientation sign between an alias entity and its owner?
5. What are the incidence relations d0, d1, and d2?

The topology module is independent of geometry, metric, physics, halo communication, and time integration.

---

## 2. Non-goals

The topology module must not compute:

- Jacobian;
- cell volume;
- face area;
- edge length;
- metric tensor;
- Hodge star;
- numerical curl, div, or grad;
- B-face to B-cell reconstruction;
- J-edge to J-cell reconstruction;
- electric field reconstruction;
- Hall EMF;
- Lorentz force;
- MPI send/recv;
- PETSc vectors, matrices, SNES, or KSP objects.

Those tasks belong to grid, metric, DEC, halo, physics, or solver modules.

---

## 3. Core principle

Topology gives the combinatorial structure.

It defines:

```text
who is connected to whom
who is the boundary of whom
which duplicate entities are actually the same topological entity
which duplicate entity is the owner
whether an alias has the same or opposite orientation as the owner
```

Topology does not define:

```text
how long an edge is
how large a face is
how large a cell is
how to compute metric inner products
how to compute Hodge star
how to compute physical fluxes
```

---

## 4. Entity model

The topology module uses four entity dimensions.

| Entity | Dimension | Meaning |
|---|---:|---|
| Node | 0 | vertex |
| Edge | 1 | oriented line segment |
| Face | 2 | oriented surface |
| Cell | 3 | oriented volume |

Each entity is represented by an `EntityKey`.

Conceptually:

```cpp
enum class EntityDim {
    Node = 0,
    Edge = 1,
    Face = 2,
    Cell = 3
};

struct EntityKey {
    EntityDim dim;
    int block;
    int i;
    int j;
    int k;
    int axis;
};
```

For nodes and cells:

```text
axis = -1
```

For edges:

| axis | edge type |
|---:|---|
| 0 | EdgeXi |
| 1 | EdgeEta |
| 2 | EdgeZeta |

For faces:

| axis | face type |
|---:|---|
| 0 | FaceXi |
| 1 | FaceEta |
| 2 | FaceZeta |

---

## 5. Local edge orientation

The local positive edge directions are:

| Edge type | Positive direction |
|---|---|
| EdgeXi | from Node(i,j,k) to Node(i+1,j,k) |
| EdgeEta | from Node(i,j,k) to Node(i,j+1,k) |
| EdgeZeta | from Node(i,j,k) to Node(i,j,k+1) |

Therefore the boundary of one edge is:

```text
boundary EdgeXi(i,j,k)
  = + Node(i+1,j,k)
    - Node(i,j,k)

boundary EdgeEta(i,j,k)
  = + Node(i,j+1,k)
    - Node(i,j,k)

boundary EdgeZeta(i,j,k)
  = + Node(i,j,k+1)
    - Node(i,j,k)
```

This defines the topological d0 operator.

---

## 6. Local face orientation

The positive face normal directions are:

| Face type | Positive normal |
|---|---|
| FaceXi | +xi |
| FaceEta | +eta |
| FaceZeta | +zeta |

The boundary edge orientation must satisfy the right-hand rule with respect to the positive face normal.

The current convention is:

```text
boundary FaceXi(i,j,k)
  = + EdgeEta (i, j,   k)
    + EdgeZeta(i, j+1, k)
    - EdgeEta (i, j,   k+1)
    - EdgeZeta(i, j,   k)
```

```text
boundary FaceEta(i,j,k)
  = + EdgeZeta(i,   j, k)
    + EdgeXi  (i,   j, k+1)
    - EdgeZeta(i+1, j, k)
    - EdgeXi  (i,   j, k)
```

```text
boundary FaceZeta(i,j,k)
  = + EdgeXi (i,   j,   k)
    + EdgeEta(i+1, j,   k)
    - EdgeXi (i,   j+1, k)
    - EdgeEta(i,   j,   k)
```

This defines the topological d1 operator.

The exact convention must be used consistently by DEC, CT, and halo synchronization.

---

## 7. Local cell orientation

The positive cell orientation is `(xi, eta, zeta)`.

The boundary of a cell is:

```text
boundary Cell(i,j,k)
  = + FaceXi  (i+1, j,   k)
    - FaceXi  (i,   j,   k)
    + FaceEta (i,   j+1, k)
    - FaceEta (i,   j,   k)
    + FaceZeta(i,   j,   k+1)
    - FaceZeta(i,   j,   k)
```

This defines the topological d2 operator.

---

## 8. Incidence operators

Topology provides only metric-free incidence stencils.

The mathematical operators are:

```text
d0: C0 -> C1
d1: C1 -> C2
d2: C2 -> C3
```

Where:

| Operator | Meaning | Input | Output |
|---|---|---|---|
| d0 | edge boundary / gradient skeleton | node 0-form | edge 1-form |
| d1 | face boundary / curl skeleton | edge 1-form | face 2-form |
| d2 | cell boundary / divergence skeleton | face 2-form | cell 3-form |

Topology does not divide by edge length, face area, or cell volume.

Topology only returns signed integer stencils.

Example conceptual interface:

```cpp
struct IncidenceEntry {
    EntityKey entity;
    int sign;
};

class Incidence {
public:
    std::vector<IncidenceEntry> boundary_of_edge(const EntityKey& edge) const;
    std::vector<IncidenceEntry> boundary_of_face(const EntityKey& face) const;
    std::vector<IncidenceEntry> boundary_of_cell(const EntityKey& cell) const;
};
```

---

## 9. Required topological identities

The following identities must hold exactly at the topological level:

```text
d1 * d0 = 0
d2 * d1 = 0
```

They correspond to:

```text
curl grad = 0
div curl = 0
```

These identities must hold before metric, Hodge star, halo communication, and physics are applied.

Any violation of these identities is a topology bug.

---

## 10. Block connectivity

Block connectivity describes how block boundary patches are connected.

A connection contains:

```cpp
struct BlockFaceConnection {
    int block_l;
    int face_l;

    int block_r;
    int face_r;

    IndexTransform index_transform;
    OrientationTransform orientation_transform;
};
```

The block connectivity layer describes:

- neighboring block id;
- local face id;
- remote face id;
- local-to-remote index transformation;
- axis permutation;
- orientation reversal.

The block connectivity layer must not contain metric information.

It must not contain:

- physical coordinates;
- Jacobian;
- face area;
- edge length;
- normal vector length;
- flux;
- electric field;
- magnetic field.

---

## 11. Equivalence class

A topological equivalence class represents the same physical topological entity stored in multiple local block coordinate systems.

Examples:

- the same interface node stored by two neighboring blocks;
- the same interface edge stored by two neighboring blocks;
- the same interface face stored by two neighboring blocks;
- repeated nodes on a rotational axis;
- repeated edges on a periodic boundary.

Each equivalence class must have exactly one owner.

Conceptual interface:

```cpp
class TopologyEquiv {
public:
    EntityKey owner_of(const EntityKey& e) const;

    bool is_owner(const EntityKey& e) const;

    int sign_to_owner(const EntityKey& e) const;

    const std::vector<EntityKey>& aliases_of(const EntityKey& owner) const;
};
```

---

## 12. Owner rule

Every topological equivalence class must have one and only one owner.

The owner selection must be deterministic.

Recommended owner priority:

1. lower MPI rank;
2. lower block id;
3. lower entity dimension if needed;
4. lexicographically smaller local index;
5. lower axis id.

The exact rule may be changed, but it must be deterministic and documented.

The topology module only chooses owner/alias relationships.

The topology module does not perform MPI communication.

---

## 13. Orientation sign

The orientation sign describes whether an alias entity has the same orientation as the owner entity.

| Entity type | Possible sign |
|---|---:|
| Node | +1 |
| Edge | +1 or -1 |
| Face | +1 or -1 |
| Cell | usually +1 |

For orientation-aware fields:

```text
alias_value = sign_to_owner(alias) * owner_value
```

For non-orientation-aware fields, the sign is ignored.

Examples:

| Field type | Uses topology sign? |
|---|---|
| node scalar | no |
| cell conserved vector | no |
| edge 1-form | yes |
| face 2-form | yes |
| Cartesian vector stored on edge | usually no |
| covariant 1-form | yes |
| flux 2-form | yes |

The topology module provides the sign.

The halo module decides whether to use the sign based on the field descriptor.

---

## 14. Interaction with field

The field module defines the mathematical identity of each field.

Examples:

```text
E_edge:
  location = edge
  form_degree = 1
  orientation_aware = true
  owner_sync = true

B_face:
  location = face
  form_degree = 2
  orientation_aware = true
  owner_sync = true

U_cell:
  location = cell
  value_kind = conserved_vector
  orientation_aware = false
  owner_sync = false or normal ghost sync only
```

Topology itself does not store physical field values.

---

## 15. Interaction with halo

Halo consumes topology information.

Topology provides:

- which entities are aliases;
- which entity is owner;
- sign from alias to owner;
- block interface index mapping.

Halo performs:

- ghost copy;
- MPI send/recv;
- owner-to-alias synchronization;
- optional alias-to-owner reduction;
- applying sign when field descriptor requires it.

Topology must not allocate MPI buffers or call MPI send/recv.

---

## 16. Interaction with metric

Metric consumes topology and grid information.

Topology provides:

- entity identity;
- incidence relation;
- orientation convention.

Grid provides:

- coordinates;
- Jacobian;
- lengths;
- areas;
- volumes;
- basis vectors.

Metric constructs:

- Hodge star;
- primal/dual metric ratios;
- form-vector projection;
- reconstruction weights;
- inner products.

Topology must not compute any of these metric quantities.

---

## 17. Interaction with DEC

DEC consumes topology incidence and metric Hodge star.

Topology provides:

```text
d0 skeleton
d1 skeleton
d2 skeleton
```

Metric provides:

```text
star0
star1
star2
star3
inner products
```

DEC provides:

```text
grad
curl
div
codifferential
J = delta B
div curl diagnostic
curl grad diagnostic
energy inner product
```

Topology alone does not provide DEC operators with physical scaling.

---

## 18. Axis and collapsed entities

For rotationally generated grids, some entities may be geometrically collapsed.

Topology may mark an entity as collapsed, but it must not regularize metric values.

Allowed topology responsibilities:

- identify repeated axis nodes;
- identify collapsed edges;
- identify collapsed faces;
- build equivalence classes for repeated axis entities;
- exclude invalid entities from owner lists when necessary;
- provide stable owner selection for repeated entities.

Forbidden topology responsibilities:

- modify Hodge star;
- cap metric factors;
- regularize edge length;
- regularize dual area;
- reconstruct B;
- reconstruct J;
- add numerical dissipation.

Those operations belong to metric, DEC, FV, or physics.

---

## 19. Minimal public API

The first-stage topology API should provide:

```cpp
class Topology {
public:
    const BlockConnectivity& connectivity() const;
    const TopologyEquiv& equiv() const;
    const Incidence& incidence() const;

    EntityKey owner_of(const EntityKey& e) const;
    int sign_to_owner(const EntityKey& e) const;
    bool is_owner(const EntityKey& e) const;
};
```

This API is intentionally small.

Future extensions may add:

- entity range iterators;
- local-owned entity lists;
- collapsed entity query;
- periodic entity query;
- boundary entity query;
- interface patch query.

---

## 20. Minimal test plan

The topology module must pass the following tests.

### 20.1 Single-block incidence test

Build a small single block.

Check:

```text
d1 * d0 = 0
d2 * d1 = 0
```

This test verifies local orientation consistency.

### 20.2 Two-block interface test

Build two connected blocks.

Check:

- interface node owner uniqueness;
- interface edge owner uniqueness;
- interface face owner uniqueness;
- alias-to-owner relation is deterministic.

### 20.3 Orientation sign test

Build a block interface with reversed orientation.

Check:

- edge 1-form sign;
- face 2-form sign;
- sign consistency under owner-alias synchronization.

### 20.4 Periodic boundary test

Build periodic equivalence classes.

Check:

- periodic owner uniqueness;
- orientation sign correctness;
- no duplicated owner in one equivalence class.

### 20.5 Rotational axis test

Build repeated axis nodes.

Check:

- repeated axis nodes belong to the same equivalence class;
- collapsed edges are identified;
- collapsed faces are identified if needed;
- owner selection is stable;
- invalid entities are excluded from owner lists if required.

---

## 21. Initial refactor scope

The initial topology refactor should implement only:

1. `EntityKey`;
2. entity construction helpers;
3. local orientation convention;
4. local incidence d0, d1, d2;
5. topology identity tests;
6. placeholder block connectivity;
7. placeholder owner/alias/sign interface.

It should not yet modify the Mercury Hall-MHD solver.

The first goal is to build a reliable topology foundation before connecting it to halo, DEC, metric, and physics.
