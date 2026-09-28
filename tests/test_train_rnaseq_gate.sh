#!/usr/bin/bash -l
# Regression test: FUNANNOTATE_TRAIN must degrade to ab-initio, not hard-fail,
# when funannotate train's RNA-seq concordance gate rejects the reads.
#
# Background: funannotate 1.9.0-rc.3 samples the RNA-seq reads, maps them to
# the genome, and stops train with exit 3 when fewer than --min_rnaseq_map_rate
# (10%) map ("RNA-seq concordance gate FAILED"). That outcome is deterministic:
# a retry with more memory cannot change it. Before this fix the module treated
# exit 3 as an infra failure, retried it, and after the last retry the global
# errorStrategy 'finish' stopped the whole run from submitting new tasks.
# Found 2026-09-26 in BFD/Funannotate_benchmarking's v1.9.0-rc3_container_rust
# cell (Pyrenophora_teres_0-1: 2.73% of sampled reads mapped).
#
# Like test_train_retry_cleanup.sh, this EXTRACTS the real failure branch from
# modules/local/funannotate_train.nf and runs it against fixtures, so it tracks
# the module rather than a copy of its logic.
#
# Usage: bash tests/test_train_rnaseq_gate.sh
# No SLURM/funannotate/PASA needed -- pure bash.

set -uo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ ! -f "${PROJECT_DIR}/modules/local/funannotate_train.nf" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_DIR="$SLURM_SUBMIT_DIR"
fi
MODULE="${PROJECT_DIR}/modules/local/funannotate_train.nf"
[ -s "$MODULE" ] || { echo "FAIL: module not found: $MODULE" >&2; exit 1; }

BLOCK=$(awk '
    /if \[ "\\\$TRAIN_STATUS" -ne 0 \]; then/ { grabbing=1; depth=1; print; next }
    grabbing {
        print
        if ($0 ~ /^[[:space:]]*if[[:space:]].*then[[:space:]]*$/) depth++
        if ($0 ~ /^[[:space:]]*fi[[:space:]]*$/) {
            depth--
            if (depth == 0) exit
        }
    }
' "$MODULE")
[ -n "$BLOCK" ] || { echo "FAIL: could not extract the TRAIN_STATUS failure branch from $MODULE" >&2; exit 1; }

# Run the extracted branch once. Args: <TRAIN_STATUS> <capture-log text>.
# Sets HARNESS_EXIT and WORKDIR/TRAINDIR for the caller's checks.
run_branch() {
    WORKDIR=$(mktemp -d)
    TRAINDIR="${WORKDIR}/genome_annotation_training/TESTSTRAIN/training"
    mkdir -p "$TRAINDIR"
    local rendered
    rendered=$(echo "$BLOCK" \
        | sed 's/\\\$/\$/g' \
        | sed 's/\\\\/\\/g' \
        | sed "s#\${params.training_target}#${WORKDIR//#/\\#}/genome_annotation_training#g" \
        | sed 's/\${out}/TESTSTRAIN/g' \
        | sed 's/\${species}/Test species/g' \
        | sed 's/\${pasa_tier}/stringent/g' \
        | sed 's/\${params.pasa_mysql}/false/g')
    {
        echo '#!/usr/bin/env bash'
        echo 'set -u'
        echo 'stop_mysqldb() { :; }'
        echo 'cd "'"$WORKDIR"'"'
        echo "$rendered"
    } > "${WORKDIR}/harness.sh"
    printf '%s\n' "$2" > "${WORKDIR}/funannotate_train_capture.log"
    set +e
    TRAIN_STATUS="$1" bash "${WORKDIR}/harness.sh" > "${WORKDIR}/harness.out" 2>&1
    HARNESS_EXIT=$?
    set -e
}

FAILS=0
GATE_MSG='[Sep 26 07:24 PM]: RNA-seq concordance gate FAILED: 5,461 of 200,000 sampled reads (2.7%) map to the genome, below --min_rnaseq_map_rate (10%).'

# 1. gate failure -> exit 0, not-trainable marker, failure TSV
run_branch 3 "$GATE_MSG"
if [ "$HARNESS_EXIT" -eq 0 ]; then echo "PASS: gate failure exits 0 (no retry)"
else echo "FAIL: gate failure exited $HARNESS_EXIT, expected 0" >&2; FAILS=$((FAILS+1)); fi
if [ -f "${TRAINDIR}/.pasa_train_failed" ]; then echo "PASS: gate failure writes .pasa_train_failed"
else echo "FAIL: gate failure did not write .pasa_train_failed" >&2; FAILS=$((FAILS+1)); fi
if grep -q 'rnaseq_gate' "${WORKDIR}/TESTSTRAIN.pasa_train_failed.tsv" 2>/dev/null; then echo "PASS: failure TSV records rnaseq_gate"
else echo "FAIL: TESTSTRAIN.pasa_train_failed.tsv missing or lacks rnaseq_gate" >&2; FAILS=$((FAILS+1)); fi
[ "${KEEP:-0}" = 1 ] && cat "${WORKDIR}/harness.out" >&2; rm -rf "$WORKDIR"

# 2. exit 3 WITHOUT the gate message -> still a hard failure (exit 3 propagates)
run_branch 3 "some other failure that happens to exit 3"
if [ "$HARNESS_EXIT" -eq 3 ]; then echo "PASS: exit 3 without gate message still hard-fails"
else echo "FAIL: exit 3 without gate message exited $HARNESS_EXIT, expected 3" >&2; FAILS=$((FAILS+1)); fi
[ "${KEEP:-0}" = 1 ] && cat "${WORKDIR}/harness.out" >&2; rm -rf "$WORKDIR"

# 3. gate message but a different exit code -> not treated as the gate
run_branch 1 "$GATE_MSG"
if [ "$HARNESS_EXIT" -eq 1 ]; then echo "PASS: gate message with exit 1 still hard-fails"
else echo "FAIL: gate message with exit 1 exited $HARNESS_EXIT, expected 1" >&2; FAILS=$((FAILS+1)); fi
[ "${KEEP:-0}" = 1 ] && cat "${WORKDIR}/harness.out" >&2; rm -rf "$WORKDIR"

if [ "$FAILS" -gt 0 ]; then echo "=== FAIL: $FAILS check(s) failed ===" >&2; exit 1; fi
echo "=== PASS: RNA-seq concordance gate degrades to ab-initio ==="
