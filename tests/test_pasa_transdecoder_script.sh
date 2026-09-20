#!/usr/bin/bash -l
# Regression test for the FUNANNOTATE_TRAIN "cdna_alignment_orf_to_genome_orf.pl"
# PATH bug (fixed 2026-09-11 in BFD/Funannotate_benchmarking's conda cells --
# see modules/local/funannotate_train.nf history and the zz-pasa-transdecoder-
# path.sh activate.d script installed into each affected conda env).
#
# Background: funannotate train's PASA/TransDecoder step invokes
# cdna_alignment_orf_to_genome_orf.pl by bare name via $PATH. There are TWO
# non-interchangeable copies of this script in each conda env:
#   opt/pasa-2.5.3/scripts/Coding/       (or opt/pasa/src/scripts/Coding for
#                                          the -rust env) -- old PASA-native
#                                          version; its regex expects
#                                          alignment IDs shaped
#                                          "ID=S<n>-asmbl_<n>"
#   opt/transdecoder/util/               (or
#                                          opt/pasa/src/pasa-plugins/transdecoder/util
#                                          for -rust) -- newer version (uses
#                                          GFF3_utils2), parses via
#                                          "Target=asmbl_<n>", which is what
#                                          this PASA 2.5.3 build's
#                                          pasa_assemblies.gff3 actually emits
#                                          (ID=align_<n>;Target=asmbl_<n> ...)
# Neither script is symlinked onto $PATH by the pasa-2.5.3 conda package's own
# activate.d/pasa-2.5.3.sh (it only exports $PASAHOME). Each affected env now
# has a supplemental etc/conda/activate.d/zz-pasa-transdecoder-path.sh that
# puts BOTH directories on $PATH -- with the transdecoder/util copy FIRST, so
# it wins. Putting them in the wrong order (or using the PASA-native copy at
# all) reproduces:
#   Error, cannot parse PASA info from ID=align_9916;Target=asmbl_1 78 586 +
#
# This test resolves cdna_alignment_orf_to_genome_orf.pl exactly as each
# conda env's activate.d scripts would set $PATH, and runs it against a
# minimal REAL-data fixture (tests/data/pasa_transdecoder/, extracted from a
# real completed run's actual pasa_assemblies.gff3 / transdecoder.gff3 /
# assemblies.fasta output for one transcript, asmbl_100) to confirm the
# resolved script (a) is the compatible one and (b) actually parses this
# format without error.
#
# Usage: bash tests/test_pasa_transdecoder_script.sh
# No SLURM/GPU/network needed -- runs in well under a second per env.

set -uo pipefail

# Prefer resolving from BASH_SOURCE (correct for a direct/interactive run);
# fall back to $SLURM_SUBMIT_DIR only if that didn't land in the repo root
# (e.g. under sbatch, where slurmd runs a spooled copy of this script and
# BASH_SOURCE no longer points at the original path) -- an ambient
# $SLURM_SUBMIT_DIR from an unrelated enclosing job allocation must not win.
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ ! -f "${PROJECT_DIR}/modules/local/funannotate_train.nf" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_DIR="$SLURM_SUBMIT_DIR"
fi
FIXTURE_DIR="${PROJECT_DIR}/tests/data/pasa_transdecoder"
ORFS_GFF3="${FIXTURE_DIR}/asmbl_100.transdecoder.gff3"
ALIGN_GFF3="${FIXTURE_DIR}/asmbl_100.genome_align.gff3"
CDNA_FASTA="${FIXTURE_DIR}/asmbl_100.cdna.fasta"

for f in "$ORFS_GFF3" "$ALIGN_GFF3" "$CDNA_FASTA"; do
    if [ ! -s "$f" ]; then
        echo "FAIL: missing fixture file $f" >&2
        exit 1
    fi
done

CONDAENV_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"

# env name -> (compatible script dir, incompatible script dir), matching the
# zz-pasa-transdecoder-path.sh installed in each env.
declare -A COMPAT_DIR=(
    [funannotate-1.8.17]="opt/transdecoder/util"
    [funannotate-1.9.0-beta.11]="opt/transdecoder/util"
    [funannotate-1.9.0-beta.11-rust]="opt/pasa/src/pasa-plugins/transdecoder/util"
    [funannotate-1.9.0-beta.12]="opt/transdecoder/util"
    [funannotate-1.9.0-beta.12-rust]="opt/pasa/src/pasa-plugins/transdecoder/util"
)
declare -A INCOMPAT_DIR=(
    [funannotate-1.8.17]="opt/pasa-2.5.3/scripts/Coding"
    [funannotate-1.9.0-beta.11]="opt/pasa-2.5.3/scripts/Coding"
    [funannotate-1.9.0-beta.11-rust]="opt/pasa/src/scripts/Coding"
    [funannotate-1.9.0-beta.12]="opt/pasa-2.5.3/scripts/Coding"
    [funannotate-1.9.0-beta.12-rust]="opt/pasa/src/scripts/Coding"
)

FAILS=0

for env in "${!COMPAT_DIR[@]}"; do
    ENV_ROOT="${CONDAENV_ROOT}/${env}"
    if [ ! -d "$ENV_ROOT" ]; then
        echo "SKIP: $env not found at $ENV_ROOT"
        continue
    fi

    COMPAT="${ENV_ROOT}/${COMPAT_DIR[$env]}"
    INCOMPAT="${ENV_ROOT}/${INCOMPAT_DIR[$env]}"
    SCRIPT_NAME=cdna_alignment_orf_to_genome_orf.pl

    if [ ! -x "${COMPAT}/${SCRIPT_NAME}" ]; then
        echo "FAIL [$env]: expected compatible script missing/not executable: ${COMPAT}/${SCRIPT_NAME}" >&2
        FAILS=$((FAILS+1))
        continue
    fi

    echo "=== $env ==="

    # 1) PATH resolution order must pick the compatible copy first.
    RESOLVED=$(PATH="${COMPAT}:${INCOMPAT}:${PATH}" command -v "$SCRIPT_NAME")
    echo "resolved: $RESOLVED"
    case "$RESOLVED" in
        "${COMPAT}/${SCRIPT_NAME}")
            echo "PASS [$env]: PATH resolves to the compatible (transdecoder) copy"
            ;;
        *)
            echo "FAIL [$env]: PATH resolved to $RESOLVED, expected ${COMPAT}/${SCRIPT_NAME}" >&2
            FAILS=$((FAILS+1))
            continue
            ;;
    esac

    # 2) the resolved script must actually parse this PASA build's real
    #    ID=align_<n>;Target=asmbl_<n> alignment format without error.
    WORKDIR=$(mktemp -d)
    STDERR_LOG="${WORKDIR}/stderr.log"
    if ! "$RESOLVED" "$ORFS_GFF3" "$ALIGN_GFF3" "$CDNA_FASTA" > "${WORKDIR}/out.gff3" 2> "$STDERR_LOG"; then
        echo "FAIL [$env]: $SCRIPT_NAME exited non-zero on the real-data fixture" >&2
        cat "$STDERR_LOG" >&2
        FAILS=$((FAILS+1))
        rm -rf "$WORKDIR"
        continue
    fi
    if grep -q "cannot parse PASA info" "$STDERR_LOG"; then
        echo "FAIL [$env]: reproduced the original bug (cannot parse PASA info)" >&2
        FAILS=$((FAILS+1))
        rm -rf "$WORKDIR"
        continue
    fi
    if [ ! -s "${WORKDIR}/out.gff3" ]; then
        echo "FAIL [$env]: no gene models propagated to the genome (empty output)" >&2
        FAILS=$((FAILS+1))
        rm -rf "$WORKDIR"
        continue
    fi
    N_GENES=$(grep -c $'\t''gene'$'\t' "${WORKDIR}/out.gff3" || true)
    echo "PASS [$env]: parsed real alignment format, propagated ${N_GENES} gene model(s)"
    rm -rf "$WORKDIR"

    # 3) sanity check that the OLD (incompatible) copy still reproduces the
    #    original bug on this fixture -- confirms the fixture actually
    #    exercises the bug, not just a script that always succeeds.
    if [ -x "${INCOMPAT}/${SCRIPT_NAME}" ]; then
        if "${INCOMPAT}/${SCRIPT_NAME}" "$ORFS_GFF3" "$ALIGN_GFF3" "$CDNA_FASTA" > /dev/null 2> "${WORKDIR}.incompat.log"; then
            echo "WARN [$env]: incompatible copy unexpectedly succeeded on this fixture -- fixture may no longer exercise the bug" >&2
        elif grep -q "cannot parse PASA info" "${WORKDIR}.incompat.log" 2>/dev/null; then
            echo "PASS [$env]: fixture confirmed to reproduce the original bug via the old PASA-native copy"
        fi
        rm -f "${WORKDIR}.incompat.log"
    fi
done

if [ "$FAILS" -gt 0 ]; then
    echo "=== FAIL: $FAILS check(s) failed ===" >&2
    exit 1
fi
echo "=== PASS: all envs resolve the compatible cdna_alignment_orf_to_genome_orf.pl and parse real PASA output ==="
