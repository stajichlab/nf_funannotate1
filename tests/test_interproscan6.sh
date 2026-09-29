#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 2 --mem 8gb --time 1-00:00:00
#SBATCH --job-name=test_interproscan6
#SBATCH --output=logs/test_interproscan6.%j.log

# Real-data test for INTERPROSCAN_RUN (InterProScan 6, nested run) plus the
# funannotate side of the hand-off:
#   1. stage an existing predicted proteome into a private test target dir
#   2. run the real module via tests/nf/interproscan_run_test.nf under
#      -profile annotate,slurm,ucr_hpcc (so the task goes through the real
#      queue/clusterOptions/resource config and SLURM submission)
#   3. run funannotate's own iprscan2annotations.py (from the pipeline's
#      funannotate image) on the resulting iprscan.xml.gz and count the
#      InterPro / GO lines it extracts
#
# Usage (from the repo root):
#   sbatch tests/test_interproscan6.sh [genome_id] [source_target_dir]
#     genome_id          default Ordospora_colligata_OC4 (1,864 proteins)
#     source_target_dir  default genome_annotation (must hold
#                        <genome_id>/predict_results/<genome_id>.proteins.fa)
# Needs scripts/setup_interproscan6.sh to have been run once.

set -euo pipefail

# $BASH_SOURCE is slurmd's spooled copy under sbatch; use the submit dir.
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${PROJECT_DIR}"
GENOME="${1:-Ordospora_colligata_OC4}"
SRC_TARGET="$(readlink -f "${2:-genome_annotation}")"
RUN_ID="${SLURM_JOB_ID:-manual_$(date +%Y%m%d%H%M%S)}"
OUT="${PROJECT_DIR}/tests/output/interproscan6_test/${RUN_ID}"
TARGET="${OUT}/genome_annotation"
FUN_SIF="${FUN_SIF:-/bigdata/stajichlab/shared/lib/singularity_cache/funannotate-1.9.0-rc.1.sif}"
FUN_DB="${FUN_DB:-/bigdata/stajichlab/shared/lib/funannotate_db}"

PROT="${SRC_TARGET}/${GENOME}/predict_results/${GENOME}.proteins.fa"
[ -s "${PROT}" ] || { echo "ERROR: no proteome at ${PROT}" >&2; exit 1; }
mkdir -p "${TARGET}/${GENOME}/predict_results"
ln -sf "${PROT}" "${TARGET}/${GENOME}/predict_results/${GENOME}.proteins.fa"
echo "[test] genome=${GENOME} proteins=$(grep -c '>' "${PROT}") out=${OUT}"

module load apptainer 2>/dev/null || true

# ── 2. module under the real pipeline configs ───────────────────────────────
START=$(date +%s)
nextflow -log "${OUT}/nextflow.log" run tests/nf/interproscan_run_test.nf \
    -profile annotate,slurm,ucr_hpcc \
    -work-dir "${OUT}/work" \
    -with-trace "${OUT}/trace.tsv" \
    --target "${TARGET}" \
    --test_ids "${GENOME}" \
    -ansi-log false
END=$(date +%s)
echo "[test] pipeline wall time: $(( END - START )) s"

MISC="${TARGET}/${GENOME}/annotate_misc"
ls -la "${MISC}"
[ -s "${MISC}/iprscan.xml.gz" ] || { echo "FAIL: no iprscan.xml.gz" >&2; exit 1; }
[ -s "${MISC}/iprscan.tsv.gz" ] || { echo "FAIL: no iprscan.tsv.gz" >&2; exit 1; }

# ── 3. funannotate's parser on the IPS6 XML ──────────────────────────────────
ANN="${OUT}/annotations.iprscan.txt"
apptainer exec -B /bigdata:/bigdata --env FUNANNOTATE_DB="${FUN_DB}" "${FUN_SIF}" bash -c '
    S=$(python -c "import funannotate, os; print(os.path.join(os.path.dirname(funannotate.__file__), \"aux_scripts\", \"iprscan2annotations.py\"))")
    python "$S" "$0" "$1"' "${MISC}/iprscan.xml.gz" "${ANN}"

N_PROT=$(grep -c '>' "${PROT}")
N_TSV_PROT=$(zcat "${MISC}/iprscan.tsv.gz" | cut -f1 | sort -u | wc -l)
N_IPR=$(grep -c $'\tdb_xref\tInterPro:' "${ANN}" || true)
N_IPR_PROT=$(grep $'\tdb_xref\tInterPro:' "${ANN}" | cut -f1 | sort -u | wc -l)
N_GO=$(grep -c $'\tgo_' "${ANN}" || true)
N_GO_PROT=$(grep $'\tgo_' "${ANN}" | cut -f1 | sort -u | wc -l)
echo "[test] input proteins:                     ${N_PROT}"
echo "[test] proteins with any IPS6 match (TSV): ${N_TSV_PROT}"
echo "[test] funannotate InterPro lines:         ${N_IPR} (${N_IPR_PROT} proteins)"
echo "[test] funannotate GO lines:               ${N_GO} (${N_GO_PROT} proteins)"
if [ "${N_IPR}" -gt 0 ]; then
    echo "PASS: funannotate parsed InterPro terms from the IPS6 XML"
else
    echo "FAIL: funannotate extracted no InterPro terms" >&2
    exit 1
fi
