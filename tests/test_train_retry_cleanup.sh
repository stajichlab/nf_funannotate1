#!/usr/bin/bash -l
# Regression test for the FUNANNOTATE_TRAIN retry-cache-poisoning bug (fixed
# 2026-09-11, commit "fix(train): clear trinity_gg/hisat2 checkpoints on
# retry, not just success").
#
# Background: Trinity's genome-guided assembly checkpoints under
# $TRAINDIR/trinity_gg/ hardcode the node-local $SCRATCH path of the SLURM
# job that generated them. A retried FUNANNOTATE_TRAIN task runs as a
# brand-new job with a brand-new $SCRATCH, so resuming from those checkpoints
# died with "mkdir: cannot create directory '/scratch/.../<dead job id>':
# Permission denied" and silently reported "0 transcripts derived from
# Trinity" -- an infinite failure loop across retries, discovered in
# BFD/Funannotate_benchmarking's conda-provisioned cells. $TRAINDIR/hisat2/
# has the same staleness risk if alignment was interrupted mid-write. The fix
# clears both on the hard-failure exit path, not only on success.
#
# Rather than duplicate that shell logic here (which would drift from the
# real module), this test EXTRACTS the actual failure-branch block from
# modules/local/funannotate_train.nf, substitutes its Groovy ${...}
# interpolations with test values, and executes the extracted block for
# real against fixture directories -- so it fails if the real fix is ever
# reverted or edited away, not just if this test's own copy of the logic
# would.
#
# Usage: bash tests/test_train_retry_cleanup.sh
# No SLURM/funannotate/PASA needed -- pure bash, runs in well under a second.

set -uo pipefail

# See test_pasa_transdecoder_script.sh for why BASH_SOURCE is tried first --
# an ambient $SLURM_SUBMIT_DIR from an unrelated enclosing job allocation
# must not silently win over the actual script location.
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ ! -f "${PROJECT_DIR}/modules/local/funannotate_train.nf" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_DIR="$SLURM_SUBMIT_DIR"
fi
MODULE="${PROJECT_DIR}/modules/local/funannotate_train.nf"

if [ ! -s "$MODULE" ]; then
    echo "FAIL: module not found: $MODULE" >&2
    exit 1
fi

# Extract the failure branch: from the TRAIN_STATUS check through its
# closing "fi" at the same (4-space) indentation. Depth-tracks nested
# if/fi so the inner "PASA completed alignment" branch doesn't confuse the
# extraction boundary.
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

if [ -z "$BLOCK" ]; then
    echo "FAIL: could not extract the TRAIN_STATUS failure branch from $MODULE -- has the module structure changed?" >&2
    exit 1
fi

if ! grep -q 'rm -rf "\\\$TRAINDIR/trinity_gg"' <<<"$BLOCK"; then
    echo "FAIL: extracted failure branch no longer clears trinity_gg on failure -- the 2026-09-11 fix appears to have been reverted" >&2
    echo "--- extracted block ---" >&2
    echo "$BLOCK" >&2
    exit 1
fi
if ! grep -q 'rm -rf "\\\$TRAINDIR/hisat2"' <<<"$BLOCK"; then
    echo "FAIL: extracted failure branch no longer clears hisat2 on failure -- the 2026-09-11 fix appears to have been reverted" >&2
    exit 1
fi
echo "PASS: failure branch in $MODULE still clears trinity_gg/hisat2 (static check)"

# Turn the extracted block into a standalone, runnable shell script: unescape
# Groovy's \$ -> $ and substitute its ${...} interpolations with test values.
WORKDIR=$(mktemp -d)
TRAINDIR="${WORKDIR}/genome_annotation_training/TESTSTRAIN/training"
mkdir -p "${TRAINDIR}/hisat2" "${TRAINDIR}/trinity_gg"
: > "${TRAINDIR}/hisat2/marker.bam"
: > "${TRAINDIR}/trinity_gg/marker.cmds"

RENDERED=$(echo "$BLOCK" \
    | sed 's/\\\$/\$/g' \
    | sed "s#\${params.training_target}#${WORKDIR//#/\\#}/genome_annotation_training#g" \
    | sed 's/\${out}/TESTSTRAIN/g' \
    | sed 's/\${species}/Test species/g' \
    | sed 's/\${pasa_tier}/stringent/g' \
    | sed 's/\${params.pasa_mysql}/false/g')

HARNESS="${WORKDIR}/harness.sh"
{
    echo '#!/usr/bin/env bash'
    echo 'set -u'
    echo 'stop_mysqldb() { :; }'  # not called (pasa_mysql=false) but stubbed defensively
    echo 'cd "'"$WORKDIR"'"'
    echo "$RENDERED"
} > "$HARNESS"
chmod +x "$HARNESS"

# funannotate_train_capture.log deliberately lacks "PASA assigned ... loci"
# so this exercises the real hard-fail path, not the graceful-degrade path.
echo "some unrelated training output, no PASA summary line here" > "${WORKDIR}/funannotate_train_capture.log"

set +e
TRAIN_STATUS=1 bash "$HARNESS"
HARNESS_EXIT=$?
set -e

FAILS=0

if [ "$HARNESS_EXIT" -ne 1 ]; then
    echo "FAIL: harness exited $HARNESS_EXIT, expected 1 (the original TRAIN_STATUS must still propagate so Nextflow retries)" >&2
    FAILS=$((FAILS+1))
else
    echo "PASS: harness exit code (1) still propagates TRAIN_STATUS -- retry still triggers"
fi

if [ -d "${TRAINDIR}/trinity_gg" ]; then
    echo "FAIL: trinity_gg/ still present after a failed train attempt -- a retry would reuse its stale \$SCRATCH-tied checkpoints" >&2
    FAILS=$((FAILS+1))
else
    echo "PASS: trinity_gg/ removed on failure"
fi

if [ -d "${TRAINDIR}/hisat2" ]; then
    echo "FAIL: hisat2/ still present after a failed train attempt" >&2
    FAILS=$((FAILS+1))
else
    echo "PASS: hisat2/ removed on failure"
fi

rm -rf "$WORKDIR"

if [ "$FAILS" -gt 0 ]; then
    echo "=== FAIL: $FAILS check(s) failed ===" >&2
    exit 1
fi
echo "=== PASS: retry-cache cleanup fires on the hard-failure exit path ==="
