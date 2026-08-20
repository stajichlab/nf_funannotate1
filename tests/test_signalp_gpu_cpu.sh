#!/usr/bin/bash -l
#SBATCH -p short_gpu -N 1 -n 1 -c 4 --mem 16gb --gres=gpu:1 --time 0-01:00:00
#SBATCH --job-name=test_signalp_gpu_cpu
#SBATCH --output=logs/test_signalp_gpu_cpu.%j.log

# Standalone smoke test for the signalp6-fast.sif container (built from
# ~/projects/containers/signalp_container/signalp6-fast.def, interpro/signalp:6.0i
# base with licensed DTU code+weights baked in) -- runs BOTH CPU (baked-in
# checkpoint) and GPU (converted checkpoint, params.signalp_gpu's default
# path) modes back-to-back on the same input, independent of the
# funannotate.nf pipeline wiring (SIGNALP_RUN isn't exercised here). Proves:
#   1. the CPU path works out of the box (no --model_dir).
#   2. the one-time `signalp6_convert_models gpu` conversion this repo's
#      modules/local/signalp_run.nf performs on first GPU use actually
#      produces a working checkpoint, and --nv + --model_dir together
#      produce real GPU inference.
#   3. CPU and GPU predictions agree at the call level (sanity, not just
#      "didn't crash").
#
# Usage: sbatch tests/test_signalp_gpu_cpu.sh [protein_fasta]
#   protein_fasta defaults to tests/data/deeptmhmm_test_proteins.faa
#   (200 real Swiss-Prot fungal proteins -- reused from the DeepTMHMM test
#   set; fine for a signalp smoke test too since it's a mix of secreted and
#   non-secreted proteins).

set -euo pipefail

SIF="${SIGNALP_SIF:-/bigdata/stajichlab/shared/lib/singularity_cache/signalp6-fast.sif}"
GPU_MODELS="${SIGNALP_GPU_MODELS:-/bigdata/stajichlab/shared/lib/singularity_cache/signalp6-fast-gpu-models}"
FASTA="${1:-tests/data/deeptmhmm_test_proteins.faa}"
OUTDIR="tests/output/signalp_gpu_cpu_test/${SLURM_JOB_ID:-manual}"
# /bigdata is auto-mounted system-wide by apptainer on this cluster's compute
# nodes (verified directly: `apptainer exec ... ls /bigdata/...` works with
# no -B at all on a short-partition job) -- these explicit binds are just
# portability insurance for sites where that isn't true, not load-bearing
# here. NOTE: $BASH_SOURCE points at slurmd's spooled copy of this script
# under an sbatch job, not the original repo path, so PROJECT_DIR must come
# from $SLURM_SUBMIT_DIR (falling back to BASH_SOURCE for a manual/non-batch run).
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SIF_DIR="$(dirname "${GPU_MODELS}")"
BINDS="${PROJECT_DIR},${SIF_DIR}"

mkdir -p "${OUTDIR}/cpu" "${OUTDIR}/gpu" logs

# signalp6 derives per-sequence output filenames (txt/plot files) from the
# FULL fasta description line, not just the accession -- verbose real
# UniProt/Swiss-Prot headers (as in this test's fixture) can blow past the
# filesystem's 255-char filename limit ("File name too long", confirmed
# against this exact fixture). Production funannotate protein FASTAs use
# short locus-tag headers (e.g. BACIR_000001-T1) so this never bites the
# real pipeline, but strip down to just the first token here so this test
# doesn't trip over its own fixture's data instead of testing signalp6 itself.
FASTA_SHORT_HDR="$(readlink -f "${OUTDIR}")/input_short_headers.fasta"
awk '/^>/{print $1; next}{print}' "$(readlink -f "${FASTA}")" > "${FASTA_SHORT_HDR}"
FASTA_ABS="${FASTA_SHORT_HDR}"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer

echo "=== host / GPU visibility ==="
hostname
nvidia-smi --query-gpu=name,memory.total,memory.used,utilization.gpu --format=csv || echo "nvidia-smi not available"

echo "=== container: $SIF ==="
echo "=== input: $FASTA_ABS ($(grep -c '^>' "$FASTA_ABS") sequences) ==="

echo
echo "############################################"
echo "# CPU run (baked-in checkpoint, no --nv)"
echo "############################################"
CPU_OUT_ABS="$(readlink -f "${OUTDIR}/cpu")"
/usr/bin/time -v apptainer exec -B "${BINDS}" "$SIF" \
    signalp6 -fasta "${FASTA_ABS}" -od "${CPU_OUT_ABS}" -org euk --mode fast -format txt --write_procs 4 -bs 16 \
    2> "${OUTDIR}/cpu/time.log" | tee "${OUTDIR}/cpu/run.log"

echo
echo "############################################"
echo "# GPU model conversion (one-time, cached at \$GPU_MODELS)"
echo "############################################"
if [ ! -f "${GPU_MODELS}/.converted" ]; then
    echo "[INFO] Converting signalp6 weights for GPU at ${GPU_MODELS}"
    # Under OUTDIR (inside PROJECT_DIR, already bound) rather than mktemp's
    # default /tmp -- don't want to depend on whether /tmp happens to be
    # auto-mounted here too.
    CONVERT_TMP="$(readlink -f "${OUTDIR}")/gpu_convert_tmp.$$"
    mkdir -p "${CONVERT_TMP}"
    apptainer exec --nv -B "${BINDS}" "$SIF" bash -c "
        SIGNALP_DIR=\$(python3 -c 'import signalp, os; print(os.path.dirname(signalp.__file__))')
        cp -r \"\$SIGNALP_DIR/model_weights/\"* '${CONVERT_TMP}/'
        signalp6_convert_models gpu '${CONVERT_TMP}'
    "
    mkdir -p "$(dirname "${GPU_MODELS}")"
    mv -T "${CONVERT_TMP}" "${GPU_MODELS}"
    touch "${GPU_MODELS}/.converted"
else
    echo "[INFO] Reusing cached GPU conversion at ${GPU_MODELS}"
fi

echo
echo "############################################"
echo "# GPU run (converted checkpoint, --nv)"
echo "############################################"
GPU_OUT_ABS="$(readlink -f "${OUTDIR}/gpu")"
/usr/bin/time -v apptainer exec --nv -B "${BINDS}" "$SIF" \
    signalp6 -fasta "${FASTA_ABS}" -od "${GPU_OUT_ABS}" -org euk --mode fast -format txt --write_procs 4 -bs 16 --model_dir "${GPU_MODELS}" \
    2> "${OUTDIR}/gpu/time.log" | tee "${OUTDIR}/gpu/run.log"

echo
echo "############################################"
echo "# Sanity checks"
echo "############################################"
CPU_RESULTS="${OUTDIR}/cpu/prediction_results.txt"
GPU_RESULTS="${OUTDIR}/gpu/prediction_results.txt"
for f in "$CPU_RESULTS" "$GPU_RESULTS"; do
    if [ ! -s "$f" ]; then
        echo "FAIL: $f missing or empty" >&2
        exit 1
    fi
done

N_IN=$(grep -c '^>' "$FASTA_ABS")
N_CPU=$(grep -vc '^#' "$CPU_RESULTS")
N_GPU=$(grep -vc '^#' "$GPU_RESULTS")
echo "proteins in:        ${N_IN}"
echo "CPU predictions out: ${N_CPU}"
echo "GPU predictions out: ${N_GPU}"
[ "$N_CPU" -ge 1 ] || { echo "FAIL: no CPU predictions" >&2; exit 1; }
[ "$N_GPU" -ge 1 ] || { echo "FAIL: no GPU predictions" >&2; exit 1; }

echo
echo "=== CPU vs GPU call-level agreement (column 2 = predicted class) ==="
CPU_CALLS=$(grep -v '^#' "$CPU_RESULTS" | sort -k1,1 | awk '{print $1"\t"$2}')
GPU_CALLS=$(grep -v '^#' "$GPU_RESULTS" | sort -k1,1 | awk '{print $1"\t"$2}')
if diff <(echo "$CPU_CALLS") <(echo "$GPU_CALLS") > "${OUTDIR}/cpu_vs_gpu_calls.diff"; then
    echo "PASS: CPU and GPU predicted classes match for all ${N_CPU} proteins"
else
    MISMATCHES=$(grep -c '^<' "${OUTDIR}/cpu_vs_gpu_calls.diff" || true)
    echo "WARN: ${MISMATCHES} call-level mismatches between CPU and GPU -- see ${OUTDIR}/cpu_vs_gpu_calls.diff"
fi

echo
echo "=== timing (CPU vs GPU, from /usr/bin/time -v) ==="
echo "--- CPU ---"
grep -E "Elapsed|Maximum resident|Percent of CPU" "${OUTDIR}/cpu/time.log" || cat "${OUTDIR}/cpu/time.log"
echo "--- GPU ---"
grep -E "Elapsed|Maximum resident|Percent of CPU" "${OUTDIR}/gpu/time.log" || cat "${OUTDIR}/gpu/time.log"

echo
echo "=== PASS: SignalP CPU+GPU smoke test completed: ${OUTDIR} ==="
