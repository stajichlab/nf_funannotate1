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
# Default provisioning is the UCR HPCC institutional profile (Lmod modules) on
# SLURM. Swap axes via env vars:
#   PROVISION=singularity sbatch run_annotate.sh   # portable containers
#   PROVISION=conda       sbatch run_annotate.sh   # shared conda envs
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

# Conda provisioning axis (PROVISION=conda): conf/provision_conda.config resolves
# `conda_envs_root` from $CONDA_ENVS_ROOT (falling back to ~/.conda). The shared
# envs are built under this root by environments/conda/build_conda_env.sh; the
# directory is on /bigdata so every SLURM node can read it. Harmless for the
# ucr_hpcc/singularity axes, which ignore it.
export CONDA_ENVS_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"

PIPELINE="${PIPELINE:-stajichlab/nf_funannotate1}"
REVISION="${REVISION:-}"
EXECUTOR="${EXECUTOR:-slurm}"
PROVISION="${PROVISION:-ucr_hpcc}"

# conf/provision_ucr_hpcc.config (loaded only by the `ucr_hpcc` profile) carries
# two unrelated things bundled together: Lmod module-loading beforeScripts, AND
# the per-process SLURM clusterOptions/queue safety net (queue selection, and
# critically `--export=ALL,SCRATCH=,TMPDIR=` so each task gets its OWN
# SLURM-prolog-created node-local $SCRATCH instead of a leaked/absent one from
# the submitting shell -- see GENOME_CLEAN's comment in provision_ucr_hpcc.config
# for the "GENEMARK_RUN epidemic" this guards against). Swapping PROVISION to
# `conda`/`singularity` used to REPLACE `ucr_hpcc` in the profile list, silently
# dropping that safety net on every UCR HPCC SLURM run and reintroducing the
# node-local /tmp exhaustion failures. Always keep `ucr_hpcc` in the profile
# list (it must come before the provisioning axis so conda's/singularity's own
# process.shell / beforeScript overrides still win), and layer the requested
# PROVISION on top only when it differs.
PROFILES="annotate,${EXECUTOR},ucr_hpcc"
if [ "${PROVISION}" != "ucr_hpcc" ]; then
    PROFILES="${PROFILES},${PROVISION}"
fi

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run "${PIPELINE}" ${REVISION:+-r "${REVISION}"} \
    -profile "${PROFILES}" \
    -resume \
    "$@"
