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
# Submit from a launch directory containing samples.csv:
#   sbatch run_earlgrey.sh
#   sbatch run_earlgrey.sh --n_test 1 --cutoff_mb 150 --repeat_taxon eukarya
#
# earlgrey_mask.nf is a SECONDARY entry script (the project default is funannotate.nf),
# so it is selected explicitly with `-main-script`. The pipeline is resolved by PROJECT
# NAME from Nextflow's asset cache, NOT this script's path — sbatch copies the script to
# a spool dir, so $BASH_SOURCE is useless here.
#
#   PIPELINE   override the source (default: the published GitHub project)
#              - a local checkout for development:  PIPELINE=/path/to/checkout
#   REVISION   git branch / tag / commit to run (default: pipeline default branch)

set -euo pipefail

module load nextflow

PIPELINE="${PIPELINE:-stajichlab/nf_funannotate1}"
REVISION="${REVISION:-}"
EXECUTOR="${EXECUTOR:-slurm}"

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run "${PIPELINE}" ${REVISION:+-r "${REVISION}"} \
    -main-script earlgrey_mask.nf \
    -profile earlgrey,${EXECUTOR} \
    -resume \
    "$@"
