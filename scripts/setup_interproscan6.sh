#!/usr/bin/env bash
# setup_interproscan6.sh — one-time, per-site setup for the InterProScan 6
# step (INTERPROSCAN_RUN, params.interproscan_engine = 'ips6').
#
# InterProScan 6 is itself a Nextflow pipeline (ebi-pf-team/interproscan6)
# that runs every member-database search in its own container. INTERPROSCAN_RUN
# launches it once per genome as a nested `nextflow run`. Three shared
# resources must exist before that, or each genome's task would fetch its own
# copy (and parallel tasks would race each other doing it):
#
#   1. a pinned checkout of the IPS6 workflow   -> params.iprscan6_pipeline
#   2. its container images, pre-pulled         -> params.sif_dir
#   3. the InterPro + member-DB data            -> params.iprscan6_datadir
#
# This script does 1 and 2, and optionally 3 (--data). Re-running is safe: it
# skips whatever already exists.
#
#   bash scripts/setup_interproscan6.sh [--data] [--licensed]
#
# --licensed also pulls the images of IPS6's licensed apps (DeepTMHMM,
# Phobius, SignalP, InterPro-N). They are skipped by default because the
# pipeline runs IPS6 without them. interpro/interpro-n is not a public image.
#
# Overridable via env:
#   IPS6_VERSION   workflow release tag            (default 6.0.2.2)
#   IPS6_DIR       where to clone the workflow     (default /bigdata/stajichlab/shared/lib/interproscan6/$IPS6_VERSION)
#   SIF_DIR        apptainer image cache           (default $NXF_APPTAINER_CACHE / NXF_SINGULARITY_CACHE / APPTAINER_CACHE, else the UCR shared cache)
#   IPR_VERSION    InterPro data release (--data)  (default 110.0)
#   IPS6_DATA_ROOT data root (--data); data lands in $IPS6_DATA_ROOT/$IPR_VERSION (default /srv/projects/db/interproscan/6.0.0)
#
# Needs git, apptainer (or singularity) and outbound HTTPS. On UCR HPCC:
# `module load apptainer`.

set -euo pipefail

IPS6_VERSION="${IPS6_VERSION:-6.0.2.2}"
IPS6_DIR="${IPS6_DIR:-/bigdata/stajichlab/shared/lib/interproscan6/${IPS6_VERSION}}"
SIF_DIR="${SIF_DIR:-${NXF_APPTAINER_CACHE:-${NXF_SINGULARITY_CACHE:-${APPTAINER_CACHE:-/bigdata/stajichlab/shared/lib/singularity_cache}}}}"
IPR_VERSION="${IPR_VERSION:-110.0}"
IPS6_DATA_ROOT="${IPS6_DATA_ROOT:-/srv/projects/db/interproscan/6.0.0}"
WANT_DATA=false
WANT_LICENSED=false
for a in "$@"; do
    case "$a" in
        --data)     WANT_DATA=true ;;
        --licensed) WANT_LICENSED=true ;;
        *) echo "ERROR: unknown option $a" >&2; exit 1 ;;
    esac
done
LICENSED_RE='^interpro/(deeptmhmm|interpro-n|phobius|signalp):'

source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer 2>/dev/null || true
CT=$(command -v apptainer || command -v singularity || true)
[ -n "$CT" ] || { echo "ERROR: apptainer/singularity not on PATH" >&2; exit 1; }

# ── 1. workflow checkout ──────────────────────────────────────────────────────
if [ -f "${IPS6_DIR}/main.nf" ]; then
    echo "[ips6] workflow already at ${IPS6_DIR} ($(git -C "${IPS6_DIR}" describe --tags 2>/dev/null || echo '?'))"
else
    mkdir -p "$(dirname "${IPS6_DIR}")"
    git clone --quiet --depth 1 --branch "${IPS6_VERSION}" \
        https://github.com/ebi-pf-team/interproscan6.git "${IPS6_DIR}"
    echo "[ips6] cloned ${IPS6_VERSION} -> ${IPS6_DIR}"
fi

# ── 2. container images ───────────────────────────────────────────────────────
# Read the image list from the checkout itself, so it always matches the
# pinned release. Names follow Nextflow's own cache naming
# (interpro/hmmer:3.3 -> interpro-hmmer-3.3.img), so the nested run finds them.
mkdir -p "${SIF_DIR}"
PULL_FAILED=0
while read -r img; do
    name="$(echo "${img}" | tr '/:' '--').img"
    if ! $WANT_LICENSED && [[ "${img}" =~ ${LICENSED_RE} ]]; then
        echo "[ips6] skip   ${img} (licensed app; use --licensed)"
    elif [ -s "${SIF_DIR}/${name}" ]; then
        echo "[ips6] have   ${name}"
    else
        echo "[ips6] pull   docker://${img} -> ${SIF_DIR}/${name}"
        rm -f "${SIF_DIR}/${name}.part"
        if "$CT" pull --name "${SIF_DIR}/${name}.part" "docker://${img}" </dev/null; then
            mv "${SIF_DIR}/${name}.part" "${SIF_DIR}/${name}"
        else
            rm -f "${SIF_DIR}/${name}.part"
            echo "[ips6] WARNING: pull failed for ${img}" >&2
            PULL_FAILED=1
        fi
    fi
done < <(grep -rhoE "container +'[^']+'" "${IPS6_DIR}/modules" | sed -E "s/container +'([^']+)'/\1/" | sort -u)
[ "$PULL_FAILED" = 0 ] || { echo "ERROR: some image pulls failed (see above)" >&2; exit 1; }

# ── 3. data (optional) ────────────────────────────────────────────────────────
# IPS6 expects <datadir>/<db>/<version>/..., e.g. <datadir>/interpro/110.0/.
if $WANT_DATA; then
    DEST="${IPS6_DATA_ROOT}/${IPR_VERSION}"
    if [ -d "${DEST}/interpro/${IPR_VERSION}" ]; then
        echo "[ips6] data already at ${DEST}"
    else
        mkdir -p "${DEST}"
        cd "${DEST}"
        # Each archive is <name>-<version>.tar.gz under
        # https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/<ips major.minor>/.
        # Let IPS6 itself download: its DOWNLOAD step (maxForks 1) fetches and
        # md5-checks every missing archive. -stub-run is not used because the
        # download is a real step; a 1-sequence input keeps the scan trivial.
        printf '>probe\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n' > probe.faa
        nextflow run "${IPS6_DIR}" -profile "$(basename "$CT")" \
            --input probe.faa --datadir "${DEST}" --interpro "${IPR_VERSION}" \
            --formats tsv --outdir "${DEST}/.probe_out" -w "${DEST}/.probe_work"
        rm -rf probe.faa "${DEST}/.probe_out" "${DEST}/.probe_work" .nextflow .nextflow.log*
        echo "[ips6] data ready at ${DEST}"
    fi
fi

echo "[ips6] done. Pipeline params:"
echo "    --iprscan6_pipeline ${IPS6_DIR}"
echo "    --iprscan6_datadir  ${IPS6_DATA_ROOT}/${IPR_VERSION}  --iprscan6_interpro ${IPR_VERSION}"
