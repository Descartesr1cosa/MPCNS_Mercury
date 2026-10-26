# Preprocess cube foundation gate

## Purpose

This gate creates a small periodic structured cube with the separate Preprocess checkout, then reads it through the real Mercury Grid, topology, metric, DEC, and halo paths. It is the fast acceptance case for future MHD operators before they are run on a production Mercury mesh.

Preprocess is a read-only dependency. The gate never edits its source tree, input file, build tree, or Git state. Its build directory and generated files live under /tmp.

## Run

From the Mercury repository root:

    NP=3 BUILD_JOBS=4 bash Z0_NULL/run_preprocess_cube_validation.sh

The script enforces exactly three MPI ranks and at most four compiler jobs. Set PREPROCESS_DIR only when the sibling Preprocess checkout is elsewhere. Set KEEP_PREPROCESS_CUBE_FIXTURE=1 to retain a generated CASE for inspection.

## Fixture contract

| Setting | Value | Purpose |
| --- | --- | --- |
| periodic cube | enabled | Known periodic topology in all directions |
| cube nodes | 9 x 9 x 9 | Small local execution cost |
| split limit | 2 | Four blocks assigned over three ranks |
| generated ranks | 3 | Matches the MPI gate |
| grid output | ASCII | Uses Mercury Read_Method=0 |
| cubed sphere | disabled | Isolates structured-cube behavior |

The Preprocess executable is invoked once without mpirun, although it is built with MPI support. Its serial output path emits all three rank geometry files deterministically. Mercury then validates those files with three MPI ranks. This is intentional: during initial integration, the Preprocess parallel writer raised std::bad_alloc for this fixture after partial output; that separate code base is not modified here.

## Assertions

| Case | Coverage |
| --- | --- |
| all_vertex | Grid loading, periodic interfaces, topology quotient, DEC chain, metric equivalence, face/edge/corner halo stages, physical IO |
| all_owner_vertex | Adds reduce-to-owner and broadcast semantics |
| owner_only | Validates owner-only storage after synchronization |

The Z0 DEC chain retains exact comparisons against the legacy metric and CT stencils. Therefore this is simultaneously a generated-grid integration test and a numerical-equivalence gate.

## Next layers for MHD

1. Manufactured periodic 0/1/2-form derivatives and incidence identities.
2. Hodge and d1-transpose action tests including MPI owner reductions.
3. A short manufactured ideal-MHD advance checking conservation and divergence preservation.
4. The first acceptance case for the extracted implicit resistive operator and a future non-diagonal Hodge.
