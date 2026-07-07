#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
BUILD_DIR="${ROOT_DIR}/build/Z0_NULL"
LOG_DIR="${BUILD_DIR}/validation_logs"
NP="${NP:-3}"

if [[ ! -d "${BUILD_DIR}" ]]; then
    echo "[Z0 validation] missing build directory: ${BUILD_DIR}" >&2
    echo "Run cmake first, for example: cmake -S . -B build" >&2
    exit 1
fi

cmake --build "${BUILD_DIR}" --target z0_null

if [[ ! -d "${BUILD_DIR}/CASE" ]]; then
    echo "[Z0 validation] missing CASE directory: ${BUILD_DIR}/CASE" >&2
    exit 1
fi

mkdir -p "${LOG_DIR}"

run_case() {
    local name="$1"
    shift
    local log="${LOG_DIR}/${name}.log"

    echo "[Z0 validation] running ${name}"
    (
        cd "${BUILD_DIR}"
        env Z0_TEST_DRIVEN_SYNC=1 "$@" mpirun -np "${NP}" ./z0_null
    ) >"${log}" 2>&1

    if grep -q "\\[FAIL\\]" "${log}"; then
        echo "[Z0 validation] ${name} failed. See ${log}" >&2
        tail -n 80 "${log}" >&2
        exit 1
    fi

    if ! grep -q "Z0_NULL finished normally" "${log}"; then
        echo "[Z0 validation] ${name} did not finish normally. See ${log}" >&2
        tail -n 80 "${log}" >&2
        exit 1
    fi
}

# Mode isolation verifies registration filters; all_vertex verifies every
# active halo kind over face, edge, and vertex stages in one run.
run_case cell_vertex      Z0_HALO_MODE=cell      Z0_HALO_LEVEL=vertex
run_case node_vertex      Z0_HALO_MODE=node      Z0_HALO_LEVEL=vertex
run_case edgeforms_vertex Z0_HALO_MODE=edgeforms Z0_HALO_LEVEL=vertex
run_case faceforms_vertex Z0_HALO_MODE=faceforms Z0_HALO_LEVEL=vertex
run_case forms_vertex     Z0_HALO_MODE=forms     Z0_HALO_LEVEL=vertex
run_case all_face         Z0_HALO_MODE=all       Z0_HALO_LEVEL=face
run_case all_edge         Z0_HALO_MODE=all       Z0_HALO_LEVEL=edge
run_case all_vertex       Z0_HALO_MODE=all       Z0_HALO_LEVEL=vertex
run_case all_owner_vertex Z0_HALO_MODE=all       Z0_HALO_LEVEL=vertex Z0_OWNER_SYNC=1
run_case owner_only       Z0_HALO_MODE=all       Z0_HALO_LEVEL=vertex Z0_OWNER_SYNC=1 Z0_OWNER_ONLY=1

echo "[Z0 validation] all halo validation cases passed"
echo "[Z0 validation] logs: ${LOG_DIR}"
