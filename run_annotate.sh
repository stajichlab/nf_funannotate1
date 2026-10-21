#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 32gb --time 7-00:00:00
#SBATCH --job-name=nxf_annotate
#SBATCH --output=logs/annotate_launch.%j.log
# --mem 32gb (was 8gb): under PROVISION=singularity, Nextflow's own image
# pull/build (e.g. SETUP_AUGUSTUS_CONFIG's `apptainer pull` of the funannotate
# image) runs in THIS head-job process, not a dispatched SLURM task -- 8GB
# OOM-killed mksquashfs mid-build converting the funannotate OCI image to SIF
# (confirmed 2026-09-06).

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
# Default provisioning is singularity (params.container_funannotate, pinned to
# funannotate v1.9.0) on SLURM. On UCR HPCC the `ucr_hpcc` profile is also loaded
# (auto-detected; see SITE_PROFILE below); elsewhere only the portable profiles are. Swap axes via env vars:
#   PROVISION=ucr_hpcc    sbatch run_annotate.sh   # UCR Lmod modules (funannotate/1.8 default module)
#   PROVISION=conda       sbatch run_annotate.sh   # shared conda envs
#   EXECUTOR=local        sbatch run_annotate.sh   # head + tasks local
#   REVISION=v0.3.0       sbatch run_annotate.sh   # pin a release
#   SITE_CONFIG=$PWD/mysite.config sbatch run_annotate.sh   # any other SLURM site (partitions, sif_dir)
#   SITE_PROFILE=         sbatch run_annotate.sh   # force no institutional profile
#
# Common overrides (passed straight through to nextflow):
#   sbatch run_annotate.sh --run_annotate true --run_antismash true
#   sbatch run_annotate.sh --n_test 2 --run_sra_fetch false
#   sbatch run_annotate.sh --asmid GCA_000001405.15
#   sbatch run_annotate.sh --taxon PHYLUM:Ascomycota

set -euo pipefail

# Site setup: `module` only exists on Lmod/Environment-Modules clusters. Elsewhere
# have nextflow (and apptainer) on PATH.
if command -v module >/dev/null 2>&1; then
    module load nextflow 2>/dev/null || true
fi
# Nextflow's own image-pull/caching step (e.g. SETUP_AUGUSTUS_CONFIG's
# `apptainer pull`) runs in THIS head-job process, not a dispatched SLURM
# task, so it never sees the per-label beforeScript module loads that cover
# task execution under the singularity provisioning axis. `module load
# nextflow` only drags in the old `singularity` module as a dependency (no
# `apptainer` binary) -- load apptainer explicitly so head-job pulls work.
if command -v module >/dev/null 2>&1; then module load apptainer 2>/dev/null || true; fi

# Conda provisioning axis (PROVISION=conda): conf/provision_conda.config resolves
# `conda_envs_root` from $CONDA_ENVS_ROOT (falling back to ~/.conda). The shared
# envs are built under this root by environments/conda/build_conda_env.sh; the
# directory is on /bigdata so every SLURM node can read it. Harmless for the
# ucr_hpcc/singularity axes, which ignore it.
if [ -d /bigdata/stajichlab ]; then
    export CONDA_ENVS_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"
fi

PIPELINE="${PIPELINE:-stajichlab/nf_funannotate1}"
REVISION="${REVISION:-}"
EXECUTOR="${EXECUTOR:-slurm}"
PROVISION="${PROVISION:-singularity}"

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
# SITE_PROFILE: an institutional profile layered after the executor (UCR: ucr_hpcc).
# Auto-selected as ucr_hpcc on UCR HPCC (detected by /bigdata/stajichlab); empty
# everywhere else. SITE_CONFIG: your own site config (partitions, sif_dir, DB paths;
# start from conf/site_template.config) passed as `-c`.
if [ -z "${SITE_PROFILE+x}" ]; then
    if [ -d /bigdata/stajichlab ]; then SITE_PROFILE=ucr_hpcc; else SITE_PROFILE=""; fi
fi
SITE_CONFIG="${SITE_CONFIG:-}"

PROFILES="annotate,${EXECUTOR}"
if [ -n "${SITE_PROFILE}" ]; then
    PROFILES="${PROFILES},${SITE_PROFILE}"
fi
if [ "${PROVISION}" != "${SITE_PROFILE}" ]; then
    PROFILES="${PROFILES},${PROVISION}"
fi

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run "${PIPELINE}" ${REVISION:+-r "${REVISION}"} \
    -profile "${PROFILES}" \
    ${SITE_CONFIG:+-c "${SITE_CONFIG}"} \
    -resume \
    "$@"
