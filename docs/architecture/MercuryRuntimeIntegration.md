# Mercury Runtime Integration and Special Numerical Cases

`MercuryZ4` is a two-ion, constrained-transport (CT) Mercury magnetosphere
solver on a multi-block cubic-sphere grid.  This note explains where the new
MPCNS runtime stops and where Mercury-specific numerical rules begin.  Read it
before changing field synchronization, a magnetic boundary condition, or the
implicit resistive operator.

## Ownership and startup

The application entry point uses `MPCNS::SimulationContext` as the sole owner
of the reusable numerical substrate.  Its initialization is collective across
MPI ranks and freezes the common field/coupling registry before the solver is
created.

```text
MPI + PETSc
    -> Param
    -> SimulationContext
       Grid -> Topology -> metric Field storage -> Coupling/Halo/owner plans
    -> SingularEdgeRegistry + edge-owner pattern
    -> MercurySolver -> MercuryBoundary + Mercury physics + Mercury I/O
```

`MercurySolver` borrows `Grid`, `Topology`, `Field`, `Halo`, and `Param` from
the context.  Therefore the context must outlive the solver.  `main.cpp`
creates them in that order so destruction occurs before `PetscFinalize()` and
MPI finalization.

`MERCURY_FIELD::FieldSpecs()` is the single source of truth for the Mercury
field catalogue.  `RuntimeFieldDescriptors()` expands those specs for the new
runtime; the older `RegisterFields()` helpers remain only for compatible test
or legacy bootstrap code.  Do not duplicate the list in `main.cpp`.

The runtime intentionally has `io.enabled = false` in this application.
Mercury's `IOModule` chooses its own restart fields, normalized metadata,
archive cadence, and Tecplot/VTK output.  Disabling the optional generic I/O
does not disable runtime ownership of the grid, topology, geometry, fields, or
communication plans.

## Discrete quantities: values are integrals, not Cartesian samples

| Quantity | Location | Meaning used by the solver |
| --- | --- | --- |
| `B_xi`, `B_eta`, `B_zeta` | face | induced magnetic-flux 2-form |
| `Badd_*` | face | prescribed/background magnetic-flux 2-form |
| `E_*`, `Ehall_*`, `Eres_*` | edge | electric circulation 1-form |
| `J_*` | edge | current circulation 1-form |
| `B_cell`, `Bind_cell`, `J_cell` | cell | reconstructed Cartesian vectors |

Consequences of this distinction:

- `Calc_J_Edge()` differentiates the **induced** flux `B_*`.  `Badd_*` is
  deliberately excluded from that current calculation.
- Hall physics uses the reconstructed total `B_cell = Bind_cell + Badd_cell`,
  while the CT update advances only the induced face flux.
- A block rotation changes component signs and locations for 1-forms and
  2-forms.  Treating these arrays as Cartesian vectors when copying across a
  panel interface silently breaks CT orientation.

The typed StateView, GeometrySnapshot, LumpedHodge, ElectromagneticOperators,
FaradayOperator, and reconstruction operators carry the regular local part of
this chain.  The legacy storage adapters are deliberately zero-copy bridges;
they do not create a second magnetic state.

## Why Mercury does not use generic `SimulationContext::sync()` for every step

The runtime provides the normal staged order
`Physical -> Halo -> Coupling -> Owner/Alias`.  Mercury preserves the same
basic order, but adds equations-specific work between those stages.  In
particular it has a historical three-stage corner policy, an `Eface` special
transfer path, and an `Eres1form` compatibility group.  More importantly, it
must build or overwrite some values only after all aliases have contributed.

`MercuryBoundary::Sync_()` is thus the physics-level scheduler.  It is not a
second owner of mesh metadata: it uses the `Halo` already created by the
runtime.  Do not replace a call such as `mercury_bound_.Sync("Jedge")` with a
generic sync unless the special operations before it have been proved
equivalent on a production Mercury case.

## Cubic-sphere seams, poles, and singular physical edges

A regular edge shared by panels has several oriented aliases.  Each block may
form a slightly different floating-point candidate in its local basis.
`ReduceEdgeAliasCandidatesToOwners_()` first applies each alias orientation,
reduces the candidates, writes the canonical owner value, and then the normal
owner-to-alias synchronization publishes that value.  Skipping the reduction
would choose a panel arbitrarily and can produce a seam forcing.

Poles and extraordinary cubic-sphere edges are different.  They have a
variable number of incident faces and cannot use the regular two-face local
stencil.  `SingularEdgeRegistry` is constructed after `SimulationContext`
finishes building topology and metrics.  It provides canonical orientation,
incident entities, and dual geometry to:

1. overwrite regular-edge current with a variable-valence current assembly;
2. overwrite ideal/Hall edge EMF with the corresponding physical-edge
   assembly; and
3. keep aliases and MPI owners consistent.

Pole face-flux collapse and pole cell-current constraints are likewise
intentional constraints, not generic smoothing.  The `divB` diagnostic masks
the first off-axis shell next to a collapsed pole because that quotient value
is not a meaningful local divergence measurement.

## Induction and conducting Mercury surface

For each explicit step, the induction path is:

```text
clear E
  -> ideal + ambipolar + artificial-resistivity terms
  -> reduce regular edge aliases; assemble singular non-Hall edges
  -> assemble/sync/add Hall E
  -> impose stationary conducting-wall non-resistive condition
  -> sync Eedge
  -> FaradayOperator: RHS_B = -d1(E)
```

At a solid/fluid interface, tangential non-resistive electric circulation is
zeroed only **after** all ideal, ambipolar, artificial, and Hall terms are
assembled.  Applying the wall condition earlier lets a later contribution
leak back onto the common edge.  When Mercury resistivity is enabled, the
separate resistive substep supplies the allowed tangential `eta J` term; do
not zero that term with the non-resistive wall operation.

## Hall and implicit resistivity

Hall EMF is assembled into `Ehall_*` separately, synchronized, and then added
to `E_*`.  This separation makes it possible to apply the Hall limiter,
radial taper, Rusanov diffusion, and singular-edge replacement before it
enters CT.  A change to the Hall expression must preserve this staging and
must use the floor returned by `Hall_Num_Limiter()`; do not divide by a raw
near-vacuum electron density.

The interior Mercury-resistivity correction is an implicit, matrix-free PETSc
solve.  It advances only `Bface` after the explicit Euler update.  Its shell
operator depends on synchronized face increments and edge-current increments,
so the B-face synchronization before solve and after update is part of the
operator contract, not redundant communication.  The implementation must
remain observationally pure for a PETSc `MatMult`: temporary state is restored
after applying the operator and cannot leak into the time-integrated fields.

## Safe change checklist

Before changing the solver, answer all of these questions:

1. Is the changed array a scalar, a Cartesian vector, a primal edge 1-form,
   or a primal face 2-form?
2. Does the operation see a regular panel edge, a pole, or an extraordinary
   physical edge?  The latter two require their dedicated constraints.
3. Does it run before or after canonical owner reduction and boundary/coupling
   work?  The order is usually part of the numerical method.
4. Does a change affect induced `B`, background `Badd`, or their reconstructed
   total?  They have intentionally different roles.
5. If communication or field registration changes, is the change represented
   in `FieldSpecs()` so every MPI rank obtains the same frozen registry?

The available repository validation builds `MercuryZ4` and exercises the
generic multi-rank Z0 infrastructure.  A real Mercury time-advance comparison
still requires a decomposed `Z4_Mercury/CASE` input and reference result; do
not claim physical regression acceptance from a successful compilation alone.
