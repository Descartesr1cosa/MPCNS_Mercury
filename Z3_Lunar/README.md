# LunarZ3 cubic-sphere application

`Z3_Lunar` reuses the MPCNS topology, field, halo, metric, DEC, reconstruction,
Hall, ambipolar, CT, and I/O infrastructure. Its plasma model contains one H+
fluid only. It has no second-ion conservation system or related production and
collision sources, no interior physical-resistivity solve, and no material
coupling channels. Every loaded grid block is treated as `Fluid`.

Build from the repository root:

```sh
cmake -S . -B build
cmake --build build --target LunarZ3 -j
```

The independent input example is in `Z3_Lunar/9999setup`. The executable
prepares `Z3_Lunar/CASE/setup` from it when needed. Set `continue_calc` to `0`
for a fresh run.
