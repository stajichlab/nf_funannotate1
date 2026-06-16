#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_annotate
#SBATCH --output=logs/annotate_launch.%j.log

# Launch the eukaryotic genome annotation pipeline (funannotate.nf).
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_annotate.sh
#
# Default provisioning is Lmod modules on SLURM. Swap axes with --   nope, set the
# profile env vars below or edit the -profile line:
#   PROVISION=singularity sbatch nextflow/run_annotate.sh
#   EXECUTOR=local        sbatch nextflow/run_annotate.sh   # head + tasks local
#
# Common overrides (passed straight through to nextflow):
#   sbatch nextflow/run_annotate.sh --run_annotate true --run_antismash true
#   sbatch nextflow/run_annotate.sh --n_test 2 --run_sra_fetch false
#   sbatch nextflow/run_annotate.sh --asmid GCA_000001405.15
#   sbatch nextflow/run_annotate.sh --taxon PHYLUM:Ascomycota

set -euo pipefail

module load nextflow

EXECUTOR="${EXECUTOR:-slurm}"
PROVISION="${PROVISION:-module}"

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/funannotate.nf \
    -c nextflow/nextflow.config \
    -profile annotate,${EXECUTOR},${PROVISION} \
    -resume \
    "$@"
