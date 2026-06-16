#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_earlgrey
#SBATCH --output=logs/earlgrey_launch.%j.log

# Launch the EarlGrey curated repeat-masking pipeline (DEFERRED / not yet tuned
# for this project — see plan). Builds a curated TE library once per species on
# the best representative (> cutoff_mb) and applies it with RepeatMasker to every
# conspecific strain, writing input_clean_genomes/<asmid>.masked.fasta — the file
# funannotate.nf's tantan step consumes (its storeDir skips wherever this exists).
#
# The earlgrey profile self-provisions via module loads in its beforeScript, so no
# separate provisioning profile is needed; just pick the executor.
#
# Submit from the PROJECT ROOT directory:
#   sbatch nextflow/run_earlgrey.sh
#   sbatch nextflow/run_earlgrey.sh --n_test 1 --cutoff_mb 150 --repeat_taxon eukarya

set -euo pipefail

module load nextflow

EXECUTOR="${EXECUTOR:-slurm}"

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/earlgrey_mask.nf \
    -c nextflow/nextflow.config \
    -profile earlgrey,${EXECUTOR} \
    -resume \
    "$@"
