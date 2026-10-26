#!/usr/bin/env bash
# Generate a small periodic cube with the read-only Preprocess checkout, then
# run the fundamental Grid/topology/DEC/MPI gate against that generated CASE.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
BUILD_DIR="${ROOT_DIR}/build/Z0_NULL"
PREPROCESS_DIR="${PREPROCESS_DIR:-${ROOT_DIR}/../Preprocess}"
PREPROCESS_BUILD_DIR="${PREPROCESS_BUILD_DIR:-/tmp/mpcns-preprocess-cube-build}"
NP="${NP:-3}"
BUILD_JOBS="${BUILD_JOBS:-4}"
TIMEOUT_SECONDS="${PREPROCESS_CUBE_TIMEOUT:-120}"
LOG_DIR="${BUILD_DIR}/preprocess_cube_validation_logs"

if [[ "${NP}" != "3" ]]; then
    echo "[preprocess cube] this fixture is generated for exactly 3 ranks (got: ${NP})" >&2
    exit 2
fi
if ! [[ "${BUILD_JOBS}" =~ ^[1-4]$ ]]; then
    echo "[preprocess cube] BUILD_JOBS must be in [1,4] (got: ${BUILD_JOBS})" >&2
    exit 2
fi
if ! [[ "${TIMEOUT_SECONDS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "[preprocess cube] PREPROCESS_CUBE_TIMEOUT must be positive (got: ${TIMEOUT_SECONDS})" >&2
    exit 2
fi
if [[ ! -f "${PREPROCESS_DIR}/CMakeLists.txt" ||
      ! -f "${PREPROCESS_DIR}/INPUT/setup.txt" ]]; then
    echo "[preprocess cube] invalid PREPROCESS_DIR: ${PREPROCESS_DIR}" >&2
    exit 2
fi

# The build and generated files intentionally stay outside the Preprocess
# checkout. Its MPI writer is not used: serial output emits every rank file
# deterministically, while downstream Mercury validation remains MPI-based.
cmake -S "${PREPROCESS_DIR}" -B "${PREPROCESS_BUILD_DIR}" -DMPCNS_ENABLE_MPI=ON
cmake --build "${PREPROCESS_BUILD_DIR}" --target MAIN --parallel "${BUILD_JOBS}"
cmake -S "${ROOT_DIR}/Z0_NULL" -B "${BUILD_DIR}"
cmake --build "${BUILD_DIR}" --target z0_null --parallel "${BUILD_JOBS}"

WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/mpcns-preprocess-cube.XXXXXX")
if [[ "${KEEP_PREPROCESS_CUBE_FIXTURE:-0}" != "1" ]]; then
    trap 'test -n "${WORK_DIR:-}" && test "${WORK_DIR#${TMPDIR:-/tmp}/mpcns-preprocess-cube.}" != "${WORK_DIR}" && rm -rf -- "${WORK_DIR}"' EXIT
fi
mkdir -p "${WORK_DIR}/INPUT" "${WORK_DIR}/CASE/setup"

cp "${PREPROCESS_DIR}/INPUT/setup.txt" "${WORK_DIR}/INPUT/setup.txt"
perl -0pi -e '
    s/int\s+proc_num\s+\d+/int\tproc_num\t3/;
    s/string\s+gfilename\s+\S+/string\tgfilename\tcube.grd/;
    s/string\s+bfilename\s+\S+/string\tbfilename\tcube.inp/;
    s/string\s+ffilename\s+\S+/string\tffilename\tcube.fvbnd/;
    s/int\s+split_limit\s+\d+/int\tsplit_limit\t2/;
    s/int\s+refine_level\s+\d+/int refine_level 0/;
    s/int\s+grd_writetype\s+\d+/int\tgrd_writetype\t0/;
    s/bool\s+generate_cubic_periodic_grid\s+=\s+\d+/bool generate_cubic_periodic_grid = 1/;
    s/bool\s+generate_cubic_sphere_grid\s+=\s+\d+/bool generate_cubic_sphere_grid = 0/;
' "${WORK_DIR}/INPUT/setup.txt"

(
    cd "${WORK_DIR}"
    "${PREPROCESS_BUILD_DIR}/MAIN" > preprocess.log 2>&1
)

for rank in 0 1 2; do
    for item in "grid   ${rank}.grd" \
                "boundary_condition/boundary   ${rank}.txt" \
                "boundary_condition/inner   ${rank}.txt" \
                "boundary_condition/parallel   ${rank}.txt"; do
        if [[ ! -f "${WORK_DIR}/OUTPUT/geometry/${item}" ]]; then
            echo "[preprocess cube] missing generated geometry file: ${item}" >&2
            exit 1
        fi
    done
done

cp "${ROOT_DIR}/Z0_NULL/CASE/setup/"* "${WORK_DIR}/CASE/setup/"
cp -a "${WORK_DIR}/OUTPUT/geometry" "${WORK_DIR}/CASE/geometry"
perl -0pi -e 's/int\s+Read_Method\s+\d+/int\tRead_Method\t0/' \
    "${WORK_DIR}/CASE/setup/0_General_Setting.txt"

mkdir -p "${LOG_DIR}"
run_case() {
    local name="$1"
    shift
    local log="${LOG_DIR}/${name}.log"

    echo "[preprocess cube] running ${name}"
    (
        cd "${WORK_DIR}"
        timeout --preserve-status --signal=TERM "${TIMEOUT_SECONDS}s" \
            env Z0_TEST_DRIVEN_SYNC=1 "$@" mpirun -np "${NP}" "${BUILD_DIR}/z0_null"
    ) > "${log}" 2>&1

    if grep -q "\\[FAIL\\]" "${log}" ||
       ! grep -q "Z0_NULL finished normally" "${log}"; then
        echo "[preprocess cube] ${name} failed. See ${log}" >&2
        tail -n 100 "${log}" >&2
        exit 1
    fi
}

# all/vertex reaches face, edge, and corner halo stages. The latter two cases
# additionally exercise reduce-to-owner and owner-only storage semantics.
run_case all_vertex       Z0_HALO_MODE=all Z0_HALO_LEVEL=vertex
run_case all_owner_vertex Z0_HALO_MODE=all Z0_HALO_LEVEL=vertex Z0_OWNER_SYNC=1
run_case owner_only       Z0_HALO_MODE=all Z0_HALO_LEVEL=vertex Z0_OWNER_SYNC=1 Z0_OWNER_ONLY=1

echo "[preprocess cube] all Grid/topology/DEC/MPI cases passed"
echo "[preprocess cube] logs: ${LOG_DIR}"
if [[ "${KEEP_PREPROCESS_CUBE_FIXTURE:-0}" == "1" ]]; then
    echo "[preprocess cube] retained fixture: ${WORK_DIR}"
fi
