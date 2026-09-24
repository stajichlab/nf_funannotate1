#!/usr/bin/env bash
# Real (non-stub) FUNANNOTATE_TRAIN test with PASA on the MySQL (MariaDB)
# backend, one Penicillium citrinum genome (PECIT1), on UCR HPCC.
#
# It checks the MariaDB path that changed on 2026-09-24 (container_mariadb now
# defaults to container_funannotate). Two arms cover the two code paths in
# modules/local/funannotate_train.nf:
#
#   singularity  train runs inside the funannotate image -> "in-image" branch
#                (the image has mariadbd + mariadb-install-db). Also runs
#                SETUP_MARIADB_DATADIR with the funannotate image.
#   conda        train runs in the conda env, which has no mariadb-install-db
#                -> "sidecar" branch: MariaDB runs from params.container_mariadb
#                (= the funannotate image) via `singularity instance start`.
#
# The funannotate version is one knob. For the official release, run e.g.:
#   FUNANNOTATE_VERSION=1.9.0 bash tests/submit_train_pasa_mysql.sh
# singularity uses docker://ghcr.io/nextgenusfs/funannotate:${FUNANNOTATE_VERSION}
# conda       uses the env   ${CONDA_ENVS_ROOT}/funannotate-${FUNANNOTATE_VERSION}
# The script stops before submitting if that image tag or env does not exist.
#
# Run on a login node, from the repository checkout:
#   bash tests/submit_train_pasa_mysql.sh
# Options (environment variables):
#   FUNANNOTATE_VERSION  default 1.9.0-rc.1
#   ARMS                 default "singularity conda"
#   TEST_ROOT            default /bigdata/stajichlab/jstajich/projects/nf/test_runs/train_pasa_mysql
#   PREVIEW=1            run `nextflow run -preview` per arm instead of sbatch
#                        (checks config + workflow graph, runs no task)
# Each arm gets its own launch directory:
#   ${TEST_ROOT}/${FUNANNOTATE_VERSION}/${arm}/
# Check results afterwards with tests/check_train_pasa_mysql.sh <launch dir>.
#
# Scope: clean -> mask -> train -> predict. There is no parameter that stops
# the pipeline after train, so predict also runs; the train result is
# available as soon as FUNANNOTATE_TRAIN finishes.

set -euo pipefail

# Run from the checkout (no BASH_SOURCE: see the global SLURM note).
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
if [ ! -f "${PROJECT_DIR}/modules/local/funannotate_train.nf" ]; then
    echo "ERROR: run this from the nf_funannotate1 checkout (or set PROJECT_DIR)" >&2
    exit 1
fi

FUNANNOTATE_VERSION="${FUNANNOTATE_VERSION:-1.9.0-rc.1}"
ARMS="${ARMS:-singularity conda}"
TEST_ROOT="${TEST_ROOT:-/bigdata/stajichlab/jstajich/projects/nf/test_runs/train_pasa_mysql}"
# Same default as run_annotate.sh; exported so PREVIEW=1 resolves the same env root.
export CONDA_ENVS_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"
PREVIEW="${PREVIEW:-0}"

IMAGE="docker://ghcr.io/nextgenusfs/funannotate:${FUNANNOTATE_VERSION}"
CONDA_ENV="funannotate-${FUNANNOTATE_VERSION}"

# Test inputs: one genome, plus RNA-seq reads and the Trinity-GG assembly that
# an earlier run of this pipeline already fetched/built for this species
# (reused so the test does not download SRA data or run Trinity again).
SPECIES_TAG="Penicillium_citrinum"
SAMPLE_ROW="Penicillium citrinum,IBT 23319,PECIT1,PECIT1,eurotiomycetes_odb10,1,5077,${PROJECT_DIR}/realtest/genomes/Penicillium_citrinum_IBT_23319.GCA_028827155.1.fasta.gz"
READS_SRC="${PROJECT_DIR}/rnaseq_reads"
TRINITY_SRC="${PROJECT_DIR}/rnaseq_data/${SPECIES_TAG}.trinity-GG.fasta"

# ── Preflight ────────────────────────────────────────────────────────────────
for f in "${READS_SRC}/${SPECIES_TAG}_norm_R1.fastq.gz" "${READS_SRC}/${SPECIES_TAG}_norm_R2.fastq.gz" \
         "${TRINITY_SRC}" "${SAMPLE_ROW##*,}"; do
    [ -s "$f" ] || { echo "ERROR: missing test input: $f" >&2; exit 1; }
done

# Both arms need the image: singularity runs train in it, conda runs the
# MariaDB sidecar from it (container_mariadb follows container_funannotate).
# Anonymous GHCR tag lookup: stops here if the release tag is not published yet.
token=$(curl -fsS "https://ghcr.io/token?scope=repository:nextgenusfs/funannotate:pull" \
        | python3 -c 'import sys,json;print(json.load(sys.stdin)["token"])')
code=$(curl -s -o /dev/null -w '%{http_code}' -H "Authorization: Bearer ${token}" \
       -H 'Accept: application/vnd.oci.image.index.v1+json,application/vnd.docker.distribution.manifest.list.v2+json,application/vnd.docker.distribution.manifest.v2+json' \
       "https://ghcr.io/v2/nextgenusfs/funannotate/manifests/${FUNANNOTATE_VERSION}")
[ "$code" = "200" ] || { echo "ERROR: image tag not found on GHCR: ${IMAGE} (HTTP ${code})" >&2; exit 1; }

for arm in $ARMS; do
    case "$arm" in
        singularity) ;;
        conda)
            [ -d "${CONDA_ENVS_ROOT}/${CONDA_ENV}" ] || { echo "ERROR: conda env not found: ${CONDA_ENVS_ROOT}/${CONDA_ENV}" >&2; exit 1; }
            ;;
        *) echo "ERROR: unknown arm '$arm' (use singularity and/or conda)" >&2; exit 1 ;;
    esac
done

# ── Per-arm launch directories ───────────────────────────────────────────────
for arm in $ARMS; do
    LAUNCH="${TEST_ROOT}/${FUNANNOTATE_VERSION}/${arm}"
    if [ -e "${LAUNCH}/.nextflow.log" ] && [ "$PREVIEW" != "1" ]; then
        echo "ERROR: ${LAUNCH} already has a run. Remove it or set TEST_ROOT to start fresh." >&2
        exit 1
    fi
    mkdir -p "${LAUNCH}/rnaseq_reads" "${LAUNCH}/rnaseq_data" "${LAUNCH}/logs"

    printf '%s\n' "SPECIES,STRAIN,ASMID,LOCUSTAG,BUSCO_LINEAGE,TRANSL_TABLE,NCBI_TAXONID,GENOME" \
        "${SAMPLE_ROW}" > "${LAUNCH}/samples.csv"
    # Real files (hard links; copy if on another filesystem), not symlinks:
    # the containers only see the launch directory.
    for f in "${READS_SRC}/${SPECIES_TAG}_norm_R1.fastq.gz" "${READS_SRC}/${SPECIES_TAG}_norm_R2.fastq.gz" \
             "${READS_SRC}/${SPECIES_TAG}_norm_SE.fastq.gz"; do
        [ -e "${LAUNCH}/rnaseq_reads/$(basename "$f")" ] || cp -l "$f" "${LAUNCH}/rnaseq_reads/" 2>/dev/null || cp "$f" "${LAUNCH}/rnaseq_reads/"
    done
    [ -e "${LAUNCH}/rnaseq_data/$(basename "${TRINITY_SRC}")" ] || cp -l "${TRINITY_SRC}" "${LAUNCH}/rnaseq_data/" 2>/dev/null || cp "${TRINITY_SRC}" "${LAUNCH}/rnaseq_data/"

    cat > "${LAUNCH}/test_train_pasa_mysql.config" <<EOF
// Written by tests/submit_train_pasa_mysql.sh (arm: ${arm}, funannotate ${FUNANNOTATE_VERSION})
params {
    samples          = "${LAUNCH}/samples.csv"
    pasa_mysql       = true      // the backend under test
    run_repeatmasker = true
    run_earlgrey     = false
    run_sra_fetch    = true      // keeps RNA-seq training on; local reads are used, nothing is downloaded
    run_update       = false
    run_annotate     = false
    run_antismash    = false
    run_interpro     = false
    run_signalp      = false
    run_deeptmhmm    = false
    skip_fcs         = true      // FCS-GX is not part of this test
}
EOF

    # Record exactly what ran.
    {
        echo "date:        $(date -Is)"
        echo "arm:         ${arm}"
        echo "funannotate: ${FUNANNOTATE_VERSION}"
        echo "image:       ${IMAGE}"
        if [ "$arm" = conda ]; then echo "conda_env:   ${CONDA_ENVS_ROOT}/${CONDA_ENV}"; fi
        echo "pipeline:    ${PROJECT_DIR} @ $(git -C "${PROJECT_DIR}" rev-parse HEAD)"
        echo "uncommitted: $(git -C "${PROJECT_DIR}" status --porcelain --untracked-files=no -- modules lib conf nextflow.config subworkflows workflows funannotate.nf | wc -l) file(s) in pipeline code"
    } > "${LAUNCH}/run_info.txt"

    CONFIGS=(-c "${LAUNCH}/test_train_pasa_mysql.config" -c "${PROJECT_DIR}/conf/site_ucr_hpcc.config")
    case "$arm" in
        singularity)
            CONFIGS+=(-c "${PROJECT_DIR}/conf/site_ucr_hpcc_singularity.config")
            VERSION_ARGS=(--container_funannotate "${IMAGE}") ;;
        conda)
            VERSION_ARGS=(--conda_env "${CONDA_ENV}" --container_funannotate "${IMAGE}") ;;
    esac
    # conda arm: --container_funannotate sets the sidecar MariaDB image
    # (container_mariadb follows it).

    echo "== ${arm}: ${LAUNCH}"
    if [ "$PREVIEW" = "1" ]; then
        (
            cd "${LAUNCH}"
            source /etc/profile.d/modules.sh 2>/dev/null || true
            module load nextflow 2>/dev/null || true
            PROFILES="annotate,slurm,ucr_hpcc,${arm}"
            nextflow run "${PROJECT_DIR}" -profile "${PROFILES}" "${CONFIGS[@]}" "${VERSION_ARGS[@]}" -preview
        )
    else
        (
            cd "${LAUNCH}"
            PIPELINE="${PROJECT_DIR}" PROVISION="${arm}" \
                sbatch --job-name="nxf_trainmysql_${arm}" "${PROJECT_DIR}/run_annotate.sh" \
                "${CONFIGS[@]}" "${VERSION_ARGS[@]}"
        )
    fi
done
