#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_annotate
#SBATCH --output=logs/annotate_launch.%j.log

# Launch the eukaryotic genome annotation pipeline (funannotate.nf).
# Submit from a launch directory containing samples.csv (and lib/ assets):
#   sbatch run_annotate.sh
#
# The pipeline is resolved by PROJECT NAME from Nextflow's asset cache, NOT by this
# script's path — sbatch copies the script to a spool dir, so $BASH_SOURCE is useless
# here. `nextflow run <project>` clones/updates ~/.nextflow/assets/<project> itself.
#
#   PIPELINE   override the source (default: the published GitHub project)
#              - a local checkout for development:  PIPELINE=/path/to/checkout
#              - or the current dir:                PIPELINE=$PWD
#   REVISION   git branch / tag / commit to run (default: pipeline default branch)
#
# Default provisioning is Lmod modules on SLURM. Swap axes via env vars:
#   PROVISION=singularity sbatch run_annotate.sh
#   EXECUTOR=local        sbatch run_annotate.sh   # head + tasks local
#   REVISION=v0.1.0       sbatch run_annotate.sh   # pin a release
#
# Common overrides (passed straight through to nextflow):
#   sbatch run_annotate.sh --run_annotate true --run_antismash true
#   sbatch run_annotate.sh --n_test 2 --run_sra_fetch false
#   sbatch run_annotate.sh --asmid GCA_000001405.15
#   sbatch run_annotate.sh --taxon PHYLUM:Ascomycota

set -euo pipefail

module load nextflow

PIPELINE="${PIPELINE:-stajichlab/nf_funannotate1}"
REVISION="${REVISION:-}"
EXECUTOR="${EXECUTOR:-slurm}"
PROVISION="${PROVISION:-module}"

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run "${PIPELINE}" ${REVISION:+-r "${REVISION}"} \
    -profile annotate,${EXECUTOR},${PROVISION} \
    -resume \
    "$@"
