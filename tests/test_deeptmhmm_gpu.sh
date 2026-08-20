#!/usr/bin/bash -l
#SBATCH -p short_gpu -N 1 -n 1 -c 4 --mem 16gb --gres=gpu:1 --time 0-01:00:00
#SBATCH --job-name=test_deeptmhmm_gpu
#SBATCH --output=logs/test_deeptmhmm_gpu.%j.log

# Standalone smoke test for the DeepTMHMM-1.0.sif container (built from
# ~/projects/containers/deeptmhmm_container/deeptmhmm.def, interpro/deeptmhmm:1.0
# base with licensed code+weights baked in at /opt/deeptmhmm) on a real GPU
# node -- independent of the funannotate.nf pipeline wiring (DEEPTMHMM_ANNOTATION
# isn't exercised here; this only proves the container itself runs on our GPU
# partition and produces real TMRs.gff3 output).
#
# Usage: sbatch tests/test_deeptmhmm_gpu.sh [protein_fasta]
#   protein_fasta defaults to tests/data/deeptmhmm_test_proteins.faa
#   (200 real Swiss-Prot fungal proteins, a real subset of lib/swissprot_fungi.faa)

set -euo pipefail

SIF="${DEEPTMHMM_SIF:-/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0.sif}"
FASTA="${1:-tests/data/deeptmhmm_test_proteins.faa}"
OUTDIR="tests/output/deeptmhmm_gpu_test/${SLURM_JOB_ID:-manual}"
# /bigdata is auto-mounted system-wide by apptainer on this cluster's compute
# nodes (verified directly with no -B at all on a short-partition job) -- this
# explicit bind is just portability insurance for sites where that isn't
# true. NOTE: $BASH_SOURCE points at slurmd's spooled copy of this script
# under an sbatch job, not the original repo path, so PROJECT_DIR must come
# from $SLURM_SUBMIT_DIR (falling back to BASH_SOURCE for a manual/non-batch run).
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

FASTA_ABS="$(readlink -f "${FASTA}")"
# predict.py creates --output-dir itself and errors if it already exists, so
# OUTDIR itself must NOT be pre-created (readlink -f below resolves the
# not-yet-existing leaf fine) -- host-side log redirects below go to a
# sibling LOGDIR instead, since those DO need to exist before the shell
# opens them.
LOGDIR="${OUTDIR}.logs"
mkdir -p "$(dirname "${OUTDIR}")" "${LOGDIR}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer

echo "=== host / GPU visibility ==="
hostname
nvidia-smi --query-gpu=name,memory.total,memory.used,utilization.gpu --format=csv || echo "nvidia-smi not available"

echo "=== container: $SIF ==="
echo "=== input: $FASTA_ABS ($(grep -c '^>' "$FASTA_ABS") sequences) ==="

# Code + weights are baked into the image at /opt/deeptmhmm; predict.py loads
# its model files by path relative to cwd, so it must be run with cwd there
# regardless of the caller's cwd. --output-dir writes results directly (no
# copy-back dance needed, unlike the old biswasaneel image's predict.py).
OUTDIR_ABS="$(readlink -f "${OUTDIR}")"
/usr/bin/time -v apptainer exec --nv -B "${PROJECT_DIR}" "$SIF" \
    bash -c "cd /opt/deeptmhmm && python3 predict.py --fasta ${FASTA_ABS} --output-dir ${OUTDIR_ABS}" \
    2> "${LOGDIR}/time.log" | tee "${LOGDIR}/run.log"

echo "=== outputs ==="
ls -la "${OUTDIR}"

echo "=== TMRs.gff3 sanity check ==="
GFF3="${OUTDIR}/TMRs.gff3"
if [ ! -s "$GFF3" ]; then
    echo "FAIL: TMRs.gff3 missing or empty" >&2
    exit 1
fi
N_PROTEINS_IN=$(grep -c '^>' "$FASTA_ABS")
N_PROTEINS_OUT=$(grep -oP '^\S+' "$GFF3" | sort -u | wc -l)
echo "proteins in:  ${N_PROTEINS_IN}"
echo "proteins out: ${N_PROTEINS_OUT}"
if [ "$N_PROTEINS_OUT" -lt 1 ]; then
    echo "FAIL: no proteins found in TMRs.gff3" >&2
    exit 1
fi
head -20 "$GFF3"

echo "=== timing / peak RSS (from /usr/bin/time -v) ==="
grep -E "Elapsed|Maximum resident|Percent of CPU" "${LOGDIR}/time.log" || cat "${LOGDIR}/time.log"

echo "=== PASS: DeepTMHMM GPU smoke test completed: ${OUTDIR} ==="
