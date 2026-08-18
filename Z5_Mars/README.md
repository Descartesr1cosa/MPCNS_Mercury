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

The neutral atmosphere and EUVAC input are optional per-rank ASCII files named
`CASE/input/atmosphere_rankXXXXXX.dat`. Records follow local block order and
then `(i,j,k)` inner-cell order. Each record contains the ten reactive neutral
densities `Ar CO2 CO C H N2 N NO O O2` in cm^-3, followed by the two auxiliary
neutral columns retained by the Fortran input,
and then the fourteen Fortran EUVAC production rates in cm^-3 s^-1. The two
auxiliary columns are transported for output but do not enter `source_species`,
matching the serial code. Missing files select a source-free atmosphere,
which is useful for parallel transport/CT validation.
