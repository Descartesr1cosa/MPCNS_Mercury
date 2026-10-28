#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
BUILD_DIR="${ROOT_DIR}/build/Z0_NULL"
LOG_DIR="${BUILD_DIR}/validation_logs"
NP="${NP:-3}"
TIMEOUT_SECONDS="${Z0_VALIDATION_TIMEOUT:-180}"

if [[ "${NP}" != "3" ]]; then
    echo "[Z0 Coupling] bundled geometry requires NP=3 (got: ${NP})" >&2
    exit 2
fi
if ! [[ "${TIMEOUT_SECONDS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "[Z0 Coupling] timeout must be a positive integer" >&2
    exit 2
fi

cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}"
cmake --build "${BUILD_DIR}" --target z0_null --parallel 4

CASE_DIR="${BUILD_DIR}/CASE"
GEOMETRY_SOURCE="${CASE_DIR}/geometry_3"
GEOMETRY_ARCHIVE="${CASE_DIR}/geometry.tar.gz"
GEOMETRY_LINK="${CASE_DIR}/geometry"
if [[ ! -d "${GEOMETRY_SOURCE}" ]]; then
    if [[ ! -f "${GEOMETRY_ARCHIVE}" ]]; then
        echo "[Z0 Coupling] missing geometry archive" >&2
        exit 1
    fi
    tar -xzf "${GEOMETRY_ARCHIVE}" -C "${CASE_DIR}"
fi
if [[ ! -e "${GEOMETRY_LINK}" && ! -L "${GEOMETRY_LINK}" ]]; then
    ln -s geometry_3 "${GEOMETRY_LINK}"
elif [[ ! -L "${GEOMETRY_LINK}" ||
        "$(readlink "${GEOMETRY_LINK}")" != "geometry_3" ]]; then
    echo "[Z0 Coupling] CASE/geometry is not the managed geometry_3 link" >&2
    exit 1
fi

mkdir -p "${LOG_DIR}"
LOG="${LOG_DIR}/coupling_integration.log"
echo "[Z0 Coupling] running prepared/Legacy integration validation"
(
    cd "${BUILD_DIR}"
    timeout --preserve-status --signal=TERM "${TIMEOUT_SECONDS}s" \
        env Z0_TEST_DRIVEN_SYNC=1 \
            Z0_COUPLING_TEST=1 \
            Z0_HALO_MODE=all \
            Z0_HALO_LEVEL=vertex \
            mpirun -np "${NP}" ./z0_null
) >"${LOG}" 2>&1

if grep -q "\[FAIL\]" "${LOG}"; then
    echo "[Z0 Coupling] validation failed: ${LOG}" >&2
    tail -n 100 "${LOG}" >&2
    exit 1
fi
if ! grep -q "Coupling Legacy/Prepared exact full exchange" "${LOG}" ||
   ! grep -q "Coupling repeated execution reuses prepared workspace" "${LOG}" ||
   ! grep -q "Z0_NULL finished normally" "${LOG}"; then
    echo "[Z0 Coupling] validation did not complete: ${LOG}" >&2
    tail -n 100 "${LOG}" >&2
    exit 1
fi

echo "[Z0 Coupling] validation passed"
echo "[Z0 Coupling] log: ${LOG}"
