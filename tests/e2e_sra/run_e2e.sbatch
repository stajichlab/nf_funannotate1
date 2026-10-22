#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 32gb --time 3-00:00:00
#SBATCH --job-name=sra_e2e
#SBATCH --output=logs/sra_e2e.%j.log
# End-to-end SRA test driver (Yarrowia lipolytica + Cryptococcus deneoformans). One job runs the
# whole sequence from the launch directory; each Nextflow launch submits its own SLURM tasks.
#
#   1. fetch the two genomes                          5. rerun with the key removed
#   2. run 1: SRA query + fetch + train + predict        -> the species must be searched again
#   3. rerun after rm -rf work .nextflow              Each step is followed by a check.sh mode;
#   4. own-reads: header-only CSV + key kept             a failed check stops the job.
#
# Usage (from an EMPTY launch directory; sbatch passes the environment on):
#   mkdir -p /path/to/sra_e2e/logs && cd /path/to/sra_e2e
#   REPO=/path/to/nf_funannotate1 sbatch $REPO/tests/e2e_sra/run_e2e.sbatch
# Optional env:
#   SITE_CONFIG=/path/mysite.config   partitions / sif_dir for a site other than UCR HPCC
#   PROTEINS=/path/swissprot_fungi.faa   protein evidence (default: ./lib/swissprot_fungi.faa)
#   SIF_DIR=/path/singularity_cache   shared container cache, exported as NXF_APPTAINER_CACHE. Set it: a
#                                     sif_dir given only in a site config does not reach container_funannotate
#                                     (nextflow.config evaluates it first), so the funannotate image would be
#                                     looked up in <launch dir>/singularity_cache. UCR HPCC:
#                                     /bigdata/stajichlab/shared/lib/singularity_cache
#   PROVISION / EXECUTOR / SITE_PROFILE   as for run_annotate.sh (default: slurm + singularity;
#                                         the ucr_hpcc profile is added automatically on UCR HPCC)
#   REVISION=<tag|branch|sha>         pin a revision when PIPELINE is a project name
#   STOP_AFTER=1..5                   run only the first N steps
#   START_AT=1..5                     skip the steps before N (resume after a failed check; the earlier
#                                     steps' outputs must still be in the launch directory)
#
# Time: the first run is the long one (SRA download, Trinity, PASA, predict for two ~20 Mb
# genomes; hours). The other steps reuse its caches and take minutes plus queue time.
set -euo pipefail

: "${REPO:?set REPO=/path/to/nf_funannotate1}"
[ -f "$REPO/run_annotate.sh" ] || { echo "REPO=$REPO has no run_annotate.sh" >&2; exit 2; }
HERE="$REPO/tests/e2e_sra"
STOP_AFTER="${STOP_AFTER:-5}"
START_AT="${START_AT:-1}"
mkdir -p logs

# Inputs for the launch directory.
for f in samples.csv params.yaml fetch_genomes.sh check.sh; do
    [ -e "$f" ] || cp "$HERE/$f" .
done
PROTEINS="${PROTEINS:-$PWD/lib/swissprot_fungi.faa}"
[ -s "$PROTEINS" ] || { echo "protein evidence not found: $PROTEINS (set PROTEINS=...)" >&2; exit 2; }

[ -z "${SIF_DIR:-}" ] || export NXF_APPTAINER_CACHE="$SIF_DIR"
export PIPELINE="$REPO"                         # run this checkout, not the published project
RUN=(bash "$REPO/run_annotate.sh" -params-file params.yaml --proteins "$PROTEINS")
step() { echo; echo "=== [$(date '+%F %T')] step $1: $2"; }
skip() { [ "$START_AT" -gt "$1" ]; }       # true when step $1 is before START_AT
done_after() { if [ "$STOP_AFTER" -le "$1" ]; then echo "STOP_AFTER=$STOP_AFTER reached"; exit 0; fi; }

f=rnaseq_reads/sra_query/Yarrowia_lipolytica.sra_query

if ! skip 1; then
    step 1 "fetch genomes"
    if [ ! -s genomes/Yarrowia_lipolytica_CLIB122.fna.gz ] || [ ! -s genomes/Cryptococcus_deneoformans_JEC21.fna.gz ]; then
        bash fetch_genomes.sh
    fi
fi
done_after 1

if ! skip 2; then
    step 2 "run 1: query, fetch, train, predict"
    "${RUN[@]}"
    bash check.sh run1
fi
done_after 2

if ! skip 3; then
    step 3 "rerun with no Nextflow cache: SRA steps must come from the caches"
    rm -rf work .nextflow
    "${RUN[@]}"
    bash check.sh rerun
fi
done_after 3

if ! skip 4; then
    step 4 "own-reads path: Yarrowia header-only CSV + key kept, no NCBI query"
    [ -e "$f.csv.orig" ] || cp "$f.csv" "$f.csv.orig"
    head -1 "$f.csv.orig" > "$f.csv"
    rm -rf work .nextflow
    "${RUN[@]}"
    bash check.sh own-reads
fi
done_after 4

step 5 "key removed: Yarrowia must be searched again"
rm -f "$f.key"
rm -rf work .nextflow
"${RUN[@]}"
bash check.sh key-removed

echo; echo "=== [$(date '+%F %T')] all steps passed"
