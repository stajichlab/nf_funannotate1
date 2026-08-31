#!/usr/bin/env bash
# Build the funannotate 1.9.0-beta.10 conda env -- PERL EVM backend, no Rust.
#
# This is the "no rust" 1.9 variant: every external is a conda package
# (trinity, pasa, evidencemodeler, augustus, snap, ...) installed by the
# manifest funannotate-1.9.0-beta.10.yml, with the 1.9.0-beta.10 funannotate
# pip-installed on top (no conda package exists for the beta -- bioconda tops
# at 1.8.17).
#
# Unlike the rust variant (build_1.9_rust.sh) there are NO source builds and
# NO hand-written activate.d snippet: the conda trinity/pasa/evidencemodeler
# packages set EVM_HOME/PASAHOME/TRINITY_HOME (and AUGUSTUS_*, QUARRY_PATH,
# ZOE) through their own etc/conda/activate.d scripts on activation.
#
# Why the explicit perl-engine checks below: as of 1.9, funannotate auto-
# detects the EVM engine when FUNANNOTATE_EVM_ENGINE is unset -- if a bare
# `evidence_modeler` binary is on PATH it silently uses the Rust engine
# (funannotate/predict.py:448, aux_scripts/funannotate-runEVM.py:19). The
# conda perl EVM only installs `evidence_modeler.pl`, so the auto-detect
# lands on perl, but this script verifies that and pins the engine to perl so
# a stray rust binary in the env can never flip the backend silently.
#
# Usage:
#   ./build_1.9.sh                       # build env, verify perl stack
#   ./build_1.9.sh --refresh             # rebuild conda manifest
#   ./build_1.9.sh --sbatch              # self-submit as one SLURM job
#
# Env var knobs:
#   FUNANNOTATE_ENV     env name to build (default funannotate-1.9.0-beta.10)
#   CONDA_ENVS_ROOT     shared conda root (default /bigdata/stajichlab/shared/condaenv)
#   SLURM_PARTITION     --sbatch queue   (default stajichlab)
#   SLURM_CPUS, SLURM_MEM, SLURM_TIME    --sbatch resources (defaults below)

set -euo pipefail

# Resolve the repo root. On the login node BASH_SOURCE[0] is valid, so a plain
# interactive run works without any env vars. But `sbatch` runs the script from
# a spool copy where BASH_SOURCE[0] does NOT point at the real checkout, so a
# job must receive PROJ_ROOT via --export (the --sbatch branch below does this).
# Prefer PROJ_ROOT when set; fall back to BASH_SOURCE[0] only when the result
# actually looks like this repo (has environments/conda). Otherwise fail loudly
# instead of silently creating LOG_DIR under some spool path.
REPO_ROOT="${PROJ_ROOT:-}"
if [[ -z "${REPO_ROOT}" ]]; then
    manifest_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_ROOT="$(cd "${manifest_dir}/../.." && pwd)"
fi
if [[ ! -d "${REPO_ROOT}/environments/conda" ]]; then
    echo "ERROR: project root '${REPO_ROOT}' has no environments/conda." >&2
    echo "       On SLURM, BASH_SOURCE[0] is invalid; submit via:  ${0##*/} --sbatch" >&2
    echo "       or export PROJ_ROOT=<abs path to nf_funannotate1>" >&2
    exit 1
fi
MANIFEST_DIR="${REPO_ROOT}/environments/conda"
LOG_DIR="${REPO_ROOT}/logs/conda_builds"
mkdir -p "${LOG_DIR}"

ENV_NAME="${FUNANNOTATE_ENV:-funannotate-1.9.0-beta.10}"
ENVS_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"
PREFIX="${ENVS_ROOT}/${ENV_NAME}"

# ── Optional SLURM self-submission ───────────────────────────────────────────
if [[ "${1:-}" == "--sbatch" && -z "${SLURM_JOB_ID:-}" ]]; then
    shift
    exec sbatch --export=ALL,PROJ_ROOT="${REPO_ROOT}" \
        --job-name="build-${ENV_NAME}" \
        --partition="${SLURM_PARTITION:-stajichlab}" \
        --cpus-per-task="${SLURM_CPUS:-8}" \
        --mem="${SLURM_MEM:-32G}" \
        --time="${SLURM_TIME:-12:00:00}" \
        --output="${LOG_DIR}/${ENV_NAME}.sbatch.out" \
        --error="${LOG_DIR}/${ENV_NAME}.sbatch.err" \
        "$(realpath "$0")" "$@"
fi

REFRESH=0
[[ "${1:-}" == "--refresh" ]] && REFRESH=1

echo "============================================================"
echo "[1.9] env     : ${ENV_NAME} -> ${PREFIX}"
echo "[1.9] refresh : ${REFRESH}"
echo "============================================================"

# ── Phase 1: conda env (all externals conda-packaged, no source builds) ─────
"${MANIFEST_DIR}/build_conda_env.sh" "${ENV_NAME}" $([[ ${REFRESH} -eq 1 ]] && echo --refresh)

export CONDA_PREFIX="${PREFIX}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ── Phase 2: verify the conda perl-EVM stack, not rust ───────────────────────
echo "[1.9] verifying conda-built externals (perl EVM)..."
test -x "${PREFIX}/bin/evidence_modeler.pl" || {
    echo "[1.9] ERROR: evidence_modeler.pl missing -- conda evidencemodeler package not installed?" >&2
    exit 1
}
if [[ -x "${PREFIX}/bin/evidence_modeler" ]]; then
    echo "[1.9] ERROR: bare 'evidence_modeler' found in bin/ -- that is the RUST"
    echo "        binary; with it on PATH 1.9 auto-selects the rust engine even when"
    echo "        FUNANNOTATE_EVM_ENGINE is unset. Remove it or use the -rust env." >&2
    exit 1
fi
test -x "${PREFIX}/bin/Trinity" || {
    echo "[1.9] ERROR: Trinity not on PATH -- conda trinity package missing?" >&2
    exit 1
}
if ! compgen -G "${PREFIX}/opt/pasa*" >/dev/null; then
    echo "[1.9] ERROR: no opt/pasa* -- conda pasa package missing?" >&2
    exit 1
fi
echo "[1.9] externals OK: evidence_modeler.pl / Trinity / opt/pasa* present; no rust evidence_modeler"

# ── Phase 3: GLIBC note (low risk for conda-only env) ────────────────────────
# Solving on this host is self-constraining: the solver records the host's
# __glibc (2.28 on Rocky 8) as a virtual package and only picks packages that
# satisfy it, so unlike source-built rust binaries there is no bypass that can
# produce a GLIBC_2.38 dependency. This differs from build_1.9_rust.sh, whose
# cargo/make artifacts are not solver-checked. Nothing further to test here.
echo "[1.9] GLIBC: conda-only solve on this host is __glibc-constrained (${PREFIX} built for host glibc)"

# ── Phase 4: verify funannotate reports the perl stack ───────────────────────
echo "[1.9] running funannotate check --show-versions (FUNANNOTATE_EVM_ENGINE=perl)..."
if command -v conda >/dev/null 2>&1; then
    conda run -p "${PREFIX}" env FUNANNOTATE_EVM_ENGINE=perl funannotate check --show-versions
else
    echo "WARN: conda not on PATH -- skipping funannotate check (run manually after activating)." >&2
fi

echo ""
echo "[1.9] DONE. Activate with the pipeline:"
echo "  nextflow run ... -profile annotate,slurm,conda --conda_env ${ENV_NAME}"
echo "Note: set FUNANNOTATE_EVM_ENGINE=perl at runtime for the perl backend (1.9 auto-"
echo "      detects rust if a bare 'evidence_modeler' binary ever appears on PATH)."
