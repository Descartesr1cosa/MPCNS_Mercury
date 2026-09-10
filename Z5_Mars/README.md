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
curled with the existing CT operator. `Badd_*` contains IMF plus these dipole
fluxes; `B_*` stores only the induced field, initialized to zero. Both fresh
runs and restarts rebuild Badd from the current setup (IMF remains when dipoles
are disabled). Restarts retain the evolved induced field. Current is computed
from the induced field; physical EMFs and forces use the total field.
Version-2 checkpoints save the ten ion states, three induced fluxes and three
complete Badd fluxes, so postprocessing can recover total flux directly.
On restart, saved Badd is overwritten using the current setup and dipole cache;
it never alters the induced field. Old Mars checkpoint formats are rejected;
start a fresh run after this update. The magnetic-convention header is explicit.
`Read_Badd_Dipoles=1` reuses validated dipole-only face-flux caches in
`DATA/Badd_Dipoles/dipoleXXXX.bin`; a stale/missing cache triggers recomputation. The supplied 4840-source `Dipoles` is converted from the user-provided Langlais
(2004) M23/-20/14 solution using the accompanying `forward.c` conventions.
It retains body-fixed coordinates; crustal fields remain disabled by default.
`Crustal_Dipoles_Cutoff_km` controls the source distance cutoff on vector potential
(0 disables it; the current example uses 1500 km).
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

The default build generates only `MarsZ5`; it has no test targets and does not
require PETSc. From the MPCNS repository root:

```sh
cmake -S . -B build-mars -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPCNS_BUILD_MERCURY=OFF -DMPCNS_BUILD_MARS=ON
cmake --build build-mars -j
```

The explicit solver options above also override any Mercury setting cached by
an earlier configuration.

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
