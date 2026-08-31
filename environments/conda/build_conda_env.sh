#!/usr/bin/env bash
# Build (or refresh) one funannotate conda environment from a manifest in
# environments/conda/ into the shared, node-visible conda root.
#
# Pre-built (frozen) conda envs are the deployment unit — the 1.8.17 toolchain
# and the pip-installed master/1.9.0-beta.10 funannotate are all standard
# bioconda/conda-forge packages, so there is nothing to conda-build except the
# funannotate python module itself (which is pip, not conda-build). Building the
# whole env under a stable path makes it readable by every SLURM compute node and
# idempotently re-usable, instead of re-solving multi-GB envs per login node.
#
# Usage:
#   ./build_conda_env.sh funannotate-1.8.17            # build one env
#   ./build_conda_env.sh --all                         # build every manifest
#   ./build_conda_env.sh funannotate-1.8.17 --refresh  # rebuild, overwrite lock
#
# Envs land in $CONDA_ENVS_ROOT (default /bigdata/stajichlab/shared/condaenv).

set -euo pipefail

ENVS_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"
MANIFEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# mamba (or micromamba) is far faster for large solves than conda.
MAMBA="${MAMBA:-mamba}"
if ! command -v "$MAMBA" >/dev/null 2>&1; then
    echo "ERROR: $MAMBA not found on PATH. Set MAMBA=/path/to/mamba or conda." >&2
    exit 1
fi

build_one() {
    local manifest="$1"; shift
    local refresh=0
    [[ "${1:-}" == "--refresh" ]] && refresh=1

    local name
    name="$(basename "${manifest}" .yml)"
    local prefix="${ENVS_ROOT}/${name}"

    echo "============================================================"
    echo "[build] ${name} -> ${prefix}"
    echo "============================================================"

    if [[ -d "${prefix}" && ${refresh} -eq 0 ]]; then
        echo "[build] ${prefix} already exists; nothing to do (--refresh to rebuild)."
        return 0
    fi

    if [[ -d "${prefix}" ]]; then
        echo "[build] Refreshing ${prefix} (partial rebuild)..."
    fi

    mkdir -p "${ENVS_ROOT}"

    # -p installs to the exact shared prefix (no janky 'envs/' nesting in a
    # group dir), and -f pulls package specs + any pip block from the yaml.
    # --yes avoids interactive prompts on the compute-side build node.
    "${MAMBA}" env create --yes --prefix "${prefix}" \
        --file "${manifest}"

    echo "[build] ${name} installed at ${prefix}"
}

if [[ "${1:-}" == "--all" ]]; then
    for m in "${MANIFEST_DIR}"/funannotate-*.yml; do
        build_one "${m}"
    done
elif [[ -n "${1:-}" ]]; then
    build_one "${MANIFEST_DIR}/${1}.yml"
else
    echo "Usage: $0 <env-name> [--refresh] | --all" >&2
    exit 2
fi