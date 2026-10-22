#!/usr/bin/env bash
# LIVE test of SRA_QUERY_BATCH against real NCBI (not run in CI: needs network, and NCBI
# answers change). Runs the real process through tests/nf/sra_query_cache_test.nf with real
# edirect, on species chosen for the cases the Fungi_BFD ports were written for:
#
#   Lentinula_edodes            strain-level first taxon ID used to find nothing        -> >= 1 run
#   Cryptococcus_deneoformans   ~215 matches: efetch past the end timed out             -> completes, >= 1 run
#   Trichophyton_rubrum         top picks were THP-1 macrophage co-culture              -> SRR12096667/8 not chosen
#   Yarrowia_lipolytica         ordinary species                                         -> >= 1 run
#   Nonexistent_organismus      no matches, no taxon ID                                  -> header-only CSV + key
#
# then a second launch that must reuse every cached query (no NCBI call), and format checks
# on real NCBI output the module's awk/regexes assume (runinfo col 26 = BioSample, col 30 =
# SampleName; biosample docsum text; xtract present).
#
# edirect is taken from PATH, or from a container when EDIRECT_CONTAINER is set, e.g.
#   EDIRECT_CONTAINER=https://depot.galaxyproject.org/singularity/entrez-direct:22.4--he881be0_0
# (needs apptainer/singularity). NCBI_API_KEY in the environment raises the rate limit.
#
# Usage: bash tests/live_sra_query.sh        (several minutes; keep outputs: KEEP_T=/some/dir)
set -uo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HARNESS="${PROJECT_DIR}/tests/nf/sra_query_cache_test.nf"
command -v nextflow >/dev/null || { echo "SKIP: nextflow not found"; exit 0; }

T="${KEEP_T:-$(mktemp -d)}"; mkdir -p "$T"; [ -n "${KEEP_T:-}" ] || trap 'rm -rf "$T"' EXIT
LAUNCH="$T/launch"; SQ="$LAUNCH/rnaseq_reads/sra_query"; mkdir -p "$SQ"

# ── edirect: PATH or container ───────────────────────────────────────────────────────────
NFCFG=()
if command -v esearch >/dev/null; then
    ED() { "$@"; }
elif [ -n "${EDIRECT_CONTAINER:-}" ]; then
    RT=$(command -v apptainer || command -v singularity) || { echo "FAIL: no apptainer/singularity"; exit 1; }
    ED() { "$RT" exec "$EDIRECT_CONTAINER" "$@"; }
    cat > "$T/ed.config" <<EOC
process.container = '${EDIRECT_CONTAINER}'
apptainer.enabled = ${RT##*/} == 'apptainer'
singularity.enabled = ${RT##*/} == 'singularity'
apptainer.autoMounts = true
singularity.autoMounts = true
EOC
    NFCFG=(-c "$T/ed.config")
else
    echo "SKIP: esearch not on PATH and EDIRECT_CONTAINER not set"; exit 0
fi

# macOS has no `timeout`; use a passthrough so the module runs (no time limit then).
if ! command -v timeout >/dev/null && ! command -v gtimeout >/dev/null; then
    mkdir -p "$T/shim"; printf '#!/bin/sh\nshift\nexec "$@"\n' > "$T/shim/timeout"; chmod +x "$T/shim/timeout"
    export PATH="$T/shim:$PATH"
elif ! command -v timeout >/dev/null; then
    mkdir -p "$T/shim"; ln -sf "$(command -v gtimeout)" "$T/shim/timeout"; export PATH="$T/shim:$PATH"
fi

pass=0; fail=0; warn=0
ok()   { pass=$((pass+1)); echo "ok   $1"; }
bad()  { fail=$((fail+1)); echo "FAIL $1"; }
note() { warn=$((warn+1)); echo "WARN $1"; }
check(){ if eval "$2"; then ok "$1"; else bad "$1"; fi; }

# ── resolve taxon IDs live (no hardcoded IDs to rot) ─────────────────────────────────────
taxid() { ED esearch -db taxonomy -query "\"$1\"[Scientific Name]" </dev/null | ED efetch -format uid | head -1; }
declare -a TAGS; TAGS=(Lentinula_edodes Cryptococcus_deneoformans Trichophyton_rubrum Yarrowia_lipolytica Nonexistent_organismus)
LIST=""
for t in "${TAGS[@]}"; do
    if [ "$t" = Nonexistent_organismus ]; then id=""; else
        id=$(taxid "${t//_/ }"); [ -n "$id" ] || { echo "FAIL: no taxid for $t"; exit 1; }
    fi
    LIST="${LIST:+$LIST,}$t:$id"
done
echo "species list: $LIST"

# ── run 1: real queries ──────────────────────────────────────────────────────────────────
runnf() { ( cd "$LAUNCH" && nextflow run "$HARNESS" ${NFCFG[@]+"${NFCFG[@]}"} -ansi-log false -w "$T/work" --species_list "$LIST" > "$T/$1.log" 2>&1 ); }
start=$(date +%s); runnf run1; rc=$?; echo "run 1: exit $rc in $(( $(date +%s) - start )) s"
check "run 1 completed" '[ "$rc" -eq 0 ]'
WL="$T/work_log.txt"; cat "$T"/work/*/*/.command.log > "$WL" 2>/dev/null

runs()  { sed 1d "$SQ/$1.sra_query.csv" 2>/dev/null | cut -d, -f3; }
nruns() { runs "$1" | grep -c . ; }

for t in Lentinula_edodes Yarrowia_lipolytica; do
    check "$t: at least one run" '[ "$(nruns '"$t"')" -ge 1 ]'
    check "$t: key written"      '[ -s "$SQ/'"$t"'.sra_query.key" ]'
done
check "Cryptococcus_deneoformans: query did not fail all attempts" '! grep -q "All 3 attempts failed for Cryptococcus_deneoformans" "$WL"'
check "Cryptococcus_deneoformans: at least one run" '[ "$(nruns Cryptococcus_deneoformans)" -ge 1 ]'
check "Cryptococcus_deneoformans: key written"      '[ -s "$SQ/Cryptococcus_deneoformans.sra_query.key" ]'
if runs Trichophyton_rubrum | grep -qE 'SRR12096667|SRR12096668'; then
    bad "Trichophyton_rubrum: THP-1 co-culture runs SRR12096667/8 were chosen"
else
    ok "Trichophyton_rubrum: THP-1 co-culture runs not chosen"
fi
grep -h "excluding" "$WL" | sed 's/^/     /' | head -20

# A name and no taxon ID that match nothing: a completed empty result (Count 0), so a
# header-only CSV AND a key, not a failure.
check "Nonexistent_organismus: header-only CSV" '[ "$(nruns Nonexistent_organismus)" -eq 0 ] && [ -s "$SQ/Nonexistent_organismus.sra_query.csv" ]'
check "Nonexistent_organismus: key written (completed, not failed)" '[ "$(cat "$SQ/Nonexistent_organismus.sra_query.key" 2>/dev/null)" = "v2|" ]'

echo "--- chosen runs ---"
for t in "${TAGS[@]}"; do printf '%-28s %s\n' "$t" "$(runs "$t" | paste -sd' ' -)"; done

# ── run 2: everything cached, no NCBI call ───────────────────────────────────────────────
rm -rf "$T/work"
runnf run2; WL2="$T/work_log2.txt"; cat "$T"/work/*/*/.command.log > "$WL2" 2>/dev/null
for t in Lentinula_edodes Yarrowia_lipolytica Cryptococcus_deneoformans; do
    check "run 2: $t reused from cache" 'grep -q "Reusing cached result for '"$t"'" "$WL2"'
done
check "run 2: no esearch for a cached species" '! grep -q "Found .* paired-end accessions for \(Lentinula_edodes\|Yarrowia_lipolytica\|Cryptococcus_deneoformans\)" "$WL2"'

# ── format checks on real NCBI output the module depends on ──────────────────────────────
if ED sh -c 'command -v xtract' >/dev/null 2>&1; then ok "xtract available (needed for the Count lookup)"; else bad "xtract missing from the edirect environment"; fi
acc=$(runs Yarrowia_lipolytica | head -1)
if [ -n "$acc" ]; then
    RI=$(ED esearch -db sra -query "$acc" </dev/null | ED efetch -format runinfo | sed -n 2p)
    f26=$(printf '%s' "$RI" | awk -F, '{print $26}'); f30=$(printf '%s' "$RI" | awk -F, '{print $30}')
    check "runinfo column 26 is a BioSample accession ($f26)" '[[ "$f26" =~ ^SAM ]]'
    [ -n "$f30" ] && ok "runinfo column 30 (SampleName) present" || note "runinfo column 30 empty for $acc"
    DOC=$(ED esearch -db biosample -query "${f26}[Accession]" </dev/null | ED efetch -format docsum)
    check "BioSample docsum for $f26 is non-empty" '[ -n "$DOC" ]'
else
    note "no Yarrowia run to check runinfo columns against"
fi

echo "passed=$pass failed=$fail warnings=$warn   (outputs: $T${KEEP_T:+, kept})"
[ "$fail" -eq 0 ] || { echo "--- run 1 log tail ---"; tail -20 "$T/run1.log"; exit 1; }
