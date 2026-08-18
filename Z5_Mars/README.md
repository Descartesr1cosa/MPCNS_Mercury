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

Build from the MPCNS repository root:

```sh
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build --target MarsZ5 -j
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
