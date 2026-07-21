# MPCNS offline post-data interface

The exporter is opt-in and is connected to the Z4 solver lifecycle. It does
not change the existing restart/Tecplot/VTK cadence.

The normal CASE interface is:

```text
bool    post_static_output      1
string  post_output_path        ./DATA_bin
```

`post_static_output` writes the static files once at the end of solver setup,
after initialization, field synchronization, and reconstruction-weight setup.
It defaults to disabled when absent, preserving older CASE files. There is no
post-data call in the timestep/checkpoint output path.

## Solver-facing API

The following explicit API remains available for custom profiles or call sites:

```cpp
POST::WriteOptions options;
options.constant_fields = {"Badd_xi", "Badd_eta", "Badd_zeta",
                           "Na", "Photo_rate"};
options.normalization = {
    {"length_ref", L_ref},
    {"time_ref", time_ref},
    {"density_ref", density_ref},
    {"velocity_ref", U_ref},
    {"magnetic_field_ref", B_ref}
};
options.physical_constants = {{"gamma", gamma}};
options.species = {"H", "Na"};

solver.WritePostStaticData("./DATA_bin", options);
```

Lower-level users can call `IOModule::SetPostDataContext`,
`IOModule::WritePostStaticData` directly, or construct
`POST::PostDataWriter` over existing read-only objects.

The primary restart data is not duplicated. `manifest.json` describes the existing
`DATA/flow_fieldNNNN.bin` checkpoint format (`MPCNSRST`), including its header,
field/block records, ghost-layer behavior, value loop order, field descriptors,
and the relation to the block-local topology maps.

`OUTPUT_DEC_JEDGE` defaults to `OFF` and is intended only for debug builds. If
it is enabled at CMake time and setup parameter `output_dec_jedge` is also
true, `J_xi/J_eta/J_zeta` are appended to the normal
`DATA/flow_fieldNNNN.bin` restart whitelist. They follow Bface's existing
output cadence; there is no separate file or frequency control. Values are
written after implicit convergence, owner/alias reconciliation, physical
boundary handling, and `Jedge` sync, never from an intermediate RK/SNES stage.
In this debug mode the restart reader auto-detects either the legacy format
without Jedge or a complete `J_xi/J_eta/J_zeta` triplet. A partial triplet is
rejected, and all MPI rank files must agree on whether the triplet is present.
If the CMake option is disabled or `output_dec_jedge` is false, any
restart containing a Jedge component is rejected and only the traditional
restart layout is accepted.

Jedge is the oriented covariant Edge 1-form `J·dr` computed from induced B
only; prescribed `Badd` is excluded. Its physical scale is
`current_density_ref * length_ref`, not `current_density_ref`.

## MPI and files

Every rank writes `geometry_NNNN.bin`, `topology_NNNN.bin`,
`reconstruction_NNNN.bin`, and `constant_field_NNNN.bin`. Rank 0 writes the
single `manifest.json` after a barrier. A serial reader merges chunks by the
global int64 IDs. No field or geometry gather is performed.

Every binary file starts with the 80-byte `MPCNSBIN` v1 header. It is followed
by 56-byte section headers and tightly packed little-endian payloads. Section
headers contain name, scalar type, component count, item count, and byte count.

Block-local map sections use the prefix `bRRRR_L`, where `RRRR` is the
rank-local block index and `L` is the entity location number:

| L | entity location |
|---:|---|
| 0 | Node |
| 1 | Cell |
| 2–4 | Edge Xi/Eta/Zeta |
| 5–7 | Face Xi/Eta/Zeta |

Each prefix has `_shape`, `_gid`, `_sign`, and `_owner` sections. Values are
linearized with `i + ni * (j + nj * k)`.

Each topology chunk also has `block_metadata[rank, local_block, physics_code]`.
Physics codes and Cell flag bits are declared in `manifest.json`. Fluid fields
such as `U_H` and `U_Na` are intentionally inactive on Solid blocks; readers
must use the physics domain/Cell flags rather than treating those records as
missing Fluid state.

Edge and Face output IDs always come from the complete quotient ID maps
(`qkey_to_qid`). The compact owner-sync ID tables are a separate internal ID
space and are never serialized as entity IDs. With validation enabled, output
performs an MPI multiplicity check requiring every quotient Node/Edge/Face/Cell
ID to be emitted exactly once.

`Bcell_weights` already includes both the cell-outward sign and the local to
global face-orientation sign. Consequently a reader only performs the CSR
multiply-add against global owner-oriented face values.

PostData v3 extends the `BfaceJedge_*` scalar CSR input space so it exactly
matches the complete Bface arrays stored in restart files. IDs below the
quotient Face count remain owner-oriented physical Face IDs. Every allocated
Face slot not represented by quotient topology, including physical-boundary,
interface, parallel, and coupling ghosts, has an auxiliary globally unique
storage ID. These auxiliary IDs are decomposition-scoped and are valid only
with dynamic files having the same `mesh_uuid` and rank/block decomposition.

The topology sections `Bstore_global_id`, `Bstore_quotient_id`,
`Bstore_address`, `Bstore_sign`, `Bstore_owner`, and `Bstore_flags` map every
solver Bface storage slot to this mixed column space. `Bstore_address` contains
`[rank, rank_local_block, restart_location_code, i, j, k]`. For a quotient
slot, use the owner-oriented value (`local * Bstore_sign`). For an auxiliary
ghost slot, use the raw solver-native restart value with no orientation
transform. The geometry chunks contain matching `Bghost_*` centers, area
vectors, measures, flags, and storage addresses. `physical_patch_meta`,
`physical_patch_names`, and `interface_patch_meta` preserve the physical and
inner/parallel/coupling patch descriptions used to classify those slots.

`BfaceJedge_*` is routed to global Edge row owners and merges all shared-edge
alias contributions; singular edges use the registry's corrected inverse
Hodge and both the source-alias and final owner orientation signs used by the
registry's owner writeback. Its coefficients now include every real or ghost
Face actually read by the solver curl stencil. `JedgeJcell_*` remains a three-output CSR from global
Edge IDs to global Cell IDs, containing the final 36-component cell weights
and the Pole-ring least-squares override. Python must not reconstruct either
operator from geometry or reapply boundary conditions.

Long field names are represented by stable per-file section prefixes
`field_NNNN`; `manifest.json` maps constant-field names to these prefixes.
