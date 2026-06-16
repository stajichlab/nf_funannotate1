#!/usr/bin/env bash
# build_tools.sh — build the two Rust helper binaries used by the SRA/RNA-seq
# steps of funannotate.nf, installing them into nextflow/tools/bin/.
#
# These are NOT committed to the repo (they are dynamically-linked, platform-
# specific ELFs). Run this once on a build node before running the pipeline with
# the SRA path enabled (--run_sra_fetch true). params.fastq_hdr_script and
# params.readlen_script (conf/profile_annotate.config) point at the results.
#
#   bash nextflow/scripts/build_tools.sh
#
# Requires a Rust toolchain (cargo). On the UCR HPCC: `module load rust`.
# Pins are overridable, e.g.:  FIXHDR_REV=<sha> ENFORCE_REV=<sha> bash ... build_tools.sh

set -euo pipefail

# ── pinned source revisions ───────────────────────────────────────────────────
FIXHDR_URL="${FIXHDR_URL:-https://github.com/hyphaltip/fix_fastq_header_trinity}"
FIXHDR_REV="${FIXHDR_REV:-ddfe9bdfa5d69c151247764aec702871b10ae291}"
ENFORCE_URL="${ENFORCE_URL:-https://github.com/hyphaltip/enforce_seqpair_readlen}"
ENFORCE_REV="${ENFORCE_REV:-a2bc290fe127dacc5b910dedc25bc40d1165b0cd}"

# ── locate nextflow/tools (script lives in nextflow/scripts/) ─────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)/tools"
BIN_DIR="${TOOLS_DIR}/bin"
mkdir -p "${BIN_DIR}"

# ── ensure cargo is available ─────────────────────────────────────────────────
if ! command -v cargo >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load rust 2>/dev/null || true
fi
command -v cargo >/dev/null 2>&1 || {
    echo "[build_tools] ERROR: cargo not found. Install Rust or 'module load rust'." >&2
    exit 1
}
echo "[build_tools] using $(cargo --version)"

# Build each tool independently so one failure doesn't block the other.
rc=0

# ── fix_fastq_header_trinity → bin/fix_fastq_header_trinity ───────────────────
# Needs a Rust toolchain supporting edition 2024 (rust >= 1.85 / module load rust).
echo "[build_tools] building fix_fastq_header_trinity @ ${FIXHDR_REV}"
if cargo install --git "${FIXHDR_URL}" --rev "${FIXHDR_REV}" --root "${TOOLS_DIR}" --force; then
    echo "[build_tools] OK fix_fastq_header_trinity"
else
    echo "[build_tools] WARN: fix_fastq_header_trinity failed to build (see edition note above)." >&2
    rc=1
fi

# ── enforce_seqpair_readlen → bin/enforce_seqpair_readlen ─────────────────────
# cargo installs the binary already named `enforce_seqpair_readlen`.
echo "[build_tools] building enforce_seqpair_readlen @ ${ENFORCE_REV}"
if cargo install --git "${ENFORCE_URL}" --rev "${ENFORCE_REV}" --root "${TOOLS_DIR}" --force; then
    echo "[build_tools] OK enforce_seqpair_readlen"
else
    echo "[build_tools] WARN: enforce_seqpair_readlen failed to build." >&2
    rc=1
fi

echo "[build_tools] installed in ${BIN_DIR}:"
ls -l "${BIN_DIR}" 2>/dev/null || true
echo "[build_tools] done (rc=${rc}). funannotate.nf uses these via params.fastq_hdr_script / params.readlen_script."
exit ${rc}
