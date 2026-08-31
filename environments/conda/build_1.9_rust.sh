#!/usr/bin/env bash
# Build the funannotate 1.9.0-beta.10 "Rust EVM/PASA/Trinity" conda env end-to-end.
#
# Two phases:
#   1) Create the pure-conda env (funannotate-1.9.0-beta.10-rust.yml, no
#      trinity/pasa/evidencemodeler packages) via build_conda_env.sh, same as
#      every other frozen env.
#   2) Source-build the hyphaltip Rust `rust_optimize` forks INTO that env,
#      using the exact same install_scripts that funannotate-live's other
#      packaging paths use (pixi [task] install-externals, Dockerfile.base,
#      conda-recipe/build.sh) so the three paths can't drift apart.
#
# This scripts the steps the 1.9.0-beta.10-rust.yml header documents by hand:
#   bowtie2 -> trinity -> evm -> pasa (Dockerfile build order; bowtie2 first
#   so a regression there fails in seconds instead of after Trinity's
#   multi-minute make), then installs an etc/conda/activate.d snippet that
#   exports the rust runtime env vars, then a GLIBC floor smoke test, then
#   `funannotate check --show-versions`.
#
# Runtime env vars: AUGUSTUS_*, QUARRY_PATH, ZOE are set by the conda
# augustus/codingquarry/snap packages' OWN activate.d scripts (verified in the
# 1.8.17 env), so this only exports what the source-built forks need and the
# conda packages don't provide: FUNANNOTATE_EVM_ENGINE, EVM_HOME, PASAHOME,
# TRINITY_HOME/TRINITYHOME, PERL5LIB (PASA hooks + PerlLib), and the opt/pasa
# bin on PATH -- mirrors Dockerfile.base's ENV block for the same variables.
#
# Usage:
#   ./build_1.9_rust.sh                       # build env + rust forks
#   ./build_1.9_rust.sh --refresh             # rebuild conda manifest, redo forks
#   ./build_1.9_rust.sh --sbatch              # self-submit as one SLURM job
#
# Env var knobs:
#   FUNANNOTATE_LIVE    dir with install_scripts/ (default ~/projects/funannotate/funannotate-live)
#   FUNANNOTATE_ENV     env name to build    (default funannotate-1.9.0-beta.10-rust)
#   CONDA_ENVS_ROOT     shared conda root    (default /bigdata/stajichlab/shared/condaenv)
#   SLURM_PARTITION     --sbatch queue       (default stajichlab)
#   SLURM_CPUS, SLURM_MEM, SLURM_TIME        --sbatch resources (defaults below)

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

ENV_NAME="${FUNANNOTATE_ENV:-funannotate-1.9.0-beta.10-rust}"
ENVS_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"
PREFIX="${ENVS_ROOT}/${ENV_NAME}"
FA_LIVE="${FUNANNOTATE_LIVE:-${HOME}/projects/funannotate/funannotate-live}"
SCRIPTS="${FA_LIVE}/install_scripts"

# Four source-build scripts + the symlink helper, in Dockerfile build order.
INSTALL_SCRIPTS=(pixi_install_bowtie2.sh pixi_install_trinity.sh pixi_install_evm.sh pixi_install_pasa.sh)
SYMLINK_SCRIPT="pixi_setup_symlinks.sh"

# ── Optional SLURM self-submission ───────────────────────────────────────────
if [[ "${1:-}" == "--sbatch" && -z "${SLURM_JOB_ID:-}" ]]; then
    shift
    mkdir -p "${LOG_DIR}"
    exec sbatch --export=ALL,PROJ_ROOT="${REPO_ROOT}" \
        --job-name="build-${ENV_NAME}" \
        --partition="${SLURM_PARTITION:-stajichlab}" \
        --cpus-per-task="${SLURM_CPUS:-32}" \
        --mem="${SLURM_MEM:-64G}" \
        --time="${SLURM_TIME:-17:00:00}" \
        --output="${LOG_DIR}/${ENV_NAME}.sbatch.out" \
        --error="${LOG_DIR}/${ENV_NAME}.sbatch.err" \
        "$(realpath "$0")" "$@"
fi

REFRESH=0
[[ "${1:-}" == "--refresh" ]] && REFRESH=1

echo "============================================================"
echo "[1.9-rust] env        : ${ENV_NAME} -> ${PREFIX}"
echo "[1.9-rust] repo base  : ${FA_LIVE} (scripts from ${SCRIPTS})"
echo "[1.9-rust] refresh    : ${REFRESH}"
echo "============================================================"

for s in "${INSTALL_SCRIPTS[@]}" "${SYMLINK_SCRIPT}"; do
    if [[ ! -f "${SCRIPTS}/${s}" ]]; then
        echo "ERROR: ${SCRIPTS}/${s} not found. Set FUNANNOTATE_LIVE to the funannotate-live checkout." >&2
        exit 1
    fi
done

# ── Phase 1: pure conda env ──────────────────────────────────────────────────
"${MANIFEST_DIR}/build_conda_env.sh" "${ENV_NAME}" $([[ ${REFRESH} -eq 1 ]] && echo --refresh)

# ── Activate bound to this prefix (scripts need CONDA_PREFIX + toolchain) ────
export CONDA_PREFIX="${PREFIX}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# ── Phase 2: source-build the Rust forks, Dockerfile order ───────────────────
for s in "${INSTALL_SCRIPTS[@]}"; do
    log="${LOG_DIR}/${ENV_NAME}.$(basename "${s}" .sh).log"
    echo "[1.9-rust] running ${s} (log: ${log})..."
    if ! bash "${SCRIPTS}/${s}" >"${log}" 2>&1; then
        echo "[1.9-rust] ERROR: ${s} failed. Last 50 lines of ${log}:" >&2
        tail -50 "${log}" >&2 || true
        exit 1
    fi
    echo "[1.9-rust] ${s} OK"
done

# ── Phase 3: symlink helper (fasta) + verification of installed artifacts ────
bash "${SCRIPTS}/${SYMLINK_SCRIPT}"

echo "[1.9-rust] verifying installed artifacts..."
test -x "${PREFIX}/bin/bowtie2" && test -x "${PREFIX}/bin/bowtie2-align-s-v256" && test -x "${PREFIX}/bin/bowtie2-align-l-v256"
test -x "${PREFIX}/bin/Trinity" && test -x "${PREFIX}/bin/sam_to_read_coords"
test -x "${PREFIX}/bin/evidence_modeler" && test -f "${PREFIX}/opt/evm/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl"
test -x "${PREFIX}/opt/pasa/src/Launch_PASA_pipeline.pl" && test -d "${PREFIX}/opt/pasa/src"

# ── Phase 4: activate.d snippet (rust runtime env vars) ──────────────────────
AD_DIR="${PREFIX}/etc/conda/activate.d"
DD_DIR="${PREFIX}/etc/conda/deactivate.d"
mkdir -p "${AD_DIR}" "${DD_DIR}"

cat > "${AD_DIR}/funannotate-rust.sh" <<'EOF'
# Rust EVM / PASA / Trinity source-built into this env (install_scripts from
# funannotate-live). Sets the runtime env vars that the missing conda
# trinity/pasa/evidencemodeler packages would otherwise supply via their own
# activate.d. AUGUSTUS_* / QUARRY_PATH / ZOE come from the conda augustus /
# codingquarry / snap packages and are NOT redefined here.
export FUNANNOTATE_EVM_ENGINE="rust"
export EVM_HOME="${CONDA_PREFIX}/opt/evm"
export PASAHOME="${CONDA_PREFIX}/opt/pasa/src"
export TRINITYHOME="${CONDA_PREFIX}/opt/trinityrnaseq"
export TRINITY_HOME="${CONDA_PREFIX}/opt/trinityrnaseq"
export PERL5LIB="${CONDA_PREFIX}/opt/pasa/src/SAMPLE_HOOKS:${CONDA_PREFIX}/opt/pasa/src/PerlLib${PERL5LIB:+:${PERL5LIB}}"
export PATH="${CONDA_PREFIX}/bin:${CONDA_PREFIX}/opt/pasa/bin${PATH:+:${PATH}}"
EOF

cat > "${DD_DIR}/funannotate-rust.sh" <<'EOF'
unset FUNANNOTATE_EVM_ENGINE EVM_HOME PASAHOME TRINITYHOME TRINITY_HOME
EOF

echo "[1.9-rust] activate.d snippet written to ${AD_DIR}/funannotate-rust.sh"

# ── Phase 5: GLIBC floor smoke test ──────────────────────────────────────────
# The source-built forks are compiled with the conda-forge c/cxx toolchain and
# can require GLIBC_2.38 (see Dockerfile.base: on glibc 2.28 RHEL8 hosts every
# compiled external dies at exec -- bamsifter, ParaFly, Inchworm, Chrysalis).
# Check the max GLIBC_* symbol each opt/ binary needs vs the HOST glibc.
glibc_smoke() {
    local host_ver
    host_ver="$(ldd --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+$' || echo "0.0")"
    echo "[1.9-rust] host glibc: ${host_ver}"
    if ! command -v objdump >/dev/null 2>&1; then
        echo "WARN: objdump not found -- GLIBC smoke test skipped (binutils needed)." >&2
        return 0
    fi
    local fail=0 bin req
    while IFS= read -r bin; do
        req="$(objdump -T "${bin}" 2>/dev/null | grep -oE 'GLIBC_[0-9]+\.[0-9]+' | sort -Vu | tail -1 | cut -d_ -f2)"
        [[ -z "${req}" ]] && continue
        if awk -v a="${host_ver}" -v b="${req}" 'BEGIN{exit !(a < b)}'; then
            echo "  FAIL ${bin##*/} needs glibc ${req} > host ${host_ver}" >&2
            fail=1
        fi
    done < <(find "${PREFIX}/opt" -type f -executable 2>/dev/null; \
             find "${PREFIX}/bin" -maxdepth 1 \
                  \( -name 'evidence_modeler' -o -name 'bowtie2*' -o -name 'sam_to_read_coords' \) 2>/dev/null)
    if [[ ${fail} -eq 1 ]]; then
        echo "ERROR: build produced binaries that need a newer glibc than this host provides." >&2
        echo "  The frozen env will not run on ${host_ver} hosts -- either rebuild on a" >&2
        echo "  newer OS or run jobs on a node with glibc >= the required version." >&2
        return 1
    fi
    echo "[1.9-rust] GLIBC smoke test passed (all binaries <= host glibc ${host_ver})"
}
glibc_smoke

# ── Phase 6: verify funannotate reports the rust stack ───────────────────────
echo "[1.9-rust] running funannotate check --show-versions..."
if command -v conda >/dev/null 2>&1; then
    conda run -p "${PREFIX}" funannotate check --show-versions
else
    echo "WARN: conda not on PATH -- skipping funannotate check (run manually after activating)." >&2
fi

echo ""
echo "[1.9-rust] DONE. Activate with the pipeline:"
echo "  nextflow run ... -profile annotate,slurm,conda --conda_env ${ENV_NAME}"
