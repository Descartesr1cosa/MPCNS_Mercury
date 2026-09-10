# MarsZ5 parallel multi-ion solver

`MarsZ5` ports the ten-ion Mars model in `fortran-mars` to the MPCNS distributed
grid, field, topology, halo and I/O infrastructure. It advances Ar+, CO2+, CO+,
C+, H+, N2+, N+, NO+, O+ and O2+ as separate five-equation Euler fluids. The
face flux is the same metric-aware local Lax-Friedrichs (Rusanov) flux used by
`Z4_Mercury`. The source operator retains the Fortran photo-chemistry network,
ion-neutral drag, electron pressure, Lorentz and Hall coupling.

Magnetic flux is stored as an oriented face 2-form and advanced by the curl of
an edge electric 1-form (CT). The electric field contains ideal, Hall and
electron-pressure terms. MarsZ5 intentionally has no PETSc dependency and no
Mercury internal-resistance implicit solve.

An optional crustal initial field is read from `CASE/setup/Dipoles`. Thousands
of dipoles are summed as edge vector-potential integrals, synchronized, and
curled with the existing CT operator into fixed `Badd_*` face fluxes. Fresh
runs add these to the total magnetic field; restarts retain their saved total
flux. The supplied 4840-source fixture is synthetic and disabled by default.
Input units, coordinates, literature and data provenance are documented in
[`docs/mars/MarsCrustalDipoles.md`](../docs/mars/MarsCrustalDipoles.md).

The ten-species photochemistry is advanced by a cell-local, coupled backward-
Euler solve with an active density lower bound. Chemical momentum and energy
losses use positive diagonal point-implicit updates. Consequently chemistry no
longer limits the global time step. The explicit step is the minimum of the
per-species fluid CFL and the total-charged-fluid fast-magnetosonic CFL; no Hall
`h^2` timestep constraint is applied. Cells that fail the full local nonlinear solve
retry with up to sixteen chemistry-only substeps without reducing the global
fluid step.

Build Mars independently from the MPCNS repository root. Disabling the Mercury
target also means that configuring Mars does not require PETSc:

```sh
cmake -S . -B build-mars -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPCNS_BUILD_MERCURY=OFF -DMPCNS_BUILD_MARS=ON
cmake --build build-mars --target MarsZ5 -j
```

Place the decomposed MPCNS grid/connectivity files in `Z5_Mars/CASE/geometry`.
The setup directory is prepared automatically from `9999setup`.

The neutral atmosphere, fourteen-channel EUVAC-equivalent photoionization and
the lower ion boundary are evaluated analytically from each cell's physical
position. No fixed-grid atmosphere file is read. The solid inner boundary uses
the complete Fortran reaction network in local photochemical equilibrium; H+
retains the original 0.3-upstream-density condition. Model equations,
coefficients, references, parameters and limitations are documented in
[`docs/mars/MarsAtmospherePhotoionization.md`](../docs/mars/MarsAtmospherePhotoionization.md).
The Mercury-compatible ideal UCT, Hall Rusanov, electron-density limiter,
timestep policy and remaining discretization differences are documented in
[`docs/mars/MarsMercuryEMDiscretization.md`](../docs/mars/MarsMercuryEMDiscretization.md).
