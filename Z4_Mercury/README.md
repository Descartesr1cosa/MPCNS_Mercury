# MercuryZ4 cubic-sphere application

`Z4_Mercury` ports the Mercury two-ion MHD solver onto the current MPCNS
topology, field catalog, halo, owner-sync, metric, and DEC infrastructure.
Panel interfaces are represented by topology equivalence classes; oriented
edge/face forms are synchronized through their canonical owners.

The grid preprocessing order is important for a cubic sphere: extrapolate
coordinates, exchange panel/corner ghost coordinates, then calculate metrics.
At startup MercuryZ4 validates every physical and ghost metric point and stops
with the rank, block, and logical index if a virtual cell is non-finite or
degenerate.

Configure from the repository root with MPI compiler wrappers, then build:

```sh
cmake -S . -B build -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build --target MercuryZ4 -j
```

Put the decomposed grid and connectivity under `Z4_Mercury/CASE/geometry`.
The reference controls are in `Z4_Mercury/9999setup`; the executable prepares
`CASE/setup` from them when needed. Set `continue_calc` to `0` for a fresh run.
