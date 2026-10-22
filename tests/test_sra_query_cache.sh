#!/usr/bin/env bash
# Regression test for the SRA query cache rule in SRA_QUERY_BATCH (the .key
# sidecar, ported from Fungi_BFD). Runs the REAL process with fake esearch/efetch
# on PATH (no NCBI access, no containers), in a scratch launch directory:
#
#   1. completed query with a hit          -> CSV has a run row, key written
#   2. every attempt fails                 -> header-only CSV, NO key
#   3. header-only CSV, no key             -> queried again (the stale-empty bug)
#   4. header-only CSV, matching key       -> reused, no query
#   5. header-only CSV, key for other taxon-> queried again
#   6. CSV with a run row                  -> reused, no query
#   7. 3 matches                           -> efetch asks for -stop 3 (not 250: NCBI 400s past the end)
#   8. 900 matches                         -> efetch capped at -stop 250
#  10. taxon list 162425+227321           -> query ORs the species name and both txids; key lists both; CSV taxonid = first
#  11. same list, header-only CSV + key    -> reused
#  12. tag Fusarium_sp. (not a binomial)   -> query has only the txid terms
#   9. 0 matches                           -> efetch not called; completed query: header-only CSV + key
#
# Usage: bash tests/test_sra_query_cache.sh     (needs nextflow; ~1 min)
set -uo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HARNESS="${PROJECT_DIR}/tests/nf/sra_query_cache_test.nf"
command -v nextflow >/dev/null || { echo "SKIP: nextflow not found"; exit 0; }

T="${KEEP_T:-$(mktemp -d)}"; [ -n "${KEEP_T:-}" ] || trap 'rm -rf "$T"' EXIT
BIN="$T/bin"; mkdir -p "$BIN"
TAG=Aspergillus_nidulans; TID=162425
SQ="$T/launch/rnaseq_reads/sra_query"
CALLS="$T/calls.log"

# esearch passes through; efetch logs the call, then fails or returns one runinfo row.
printf '#!/bin/sh\ncat >/dev/null 2>&1; echo "esearch $*" >> "%s"; echo ESEARCH\n' "$CALLS" > "$BIN/esearch"
cat > "$BIN/efetch" <<EOS
#!/bin/sh
echo "efetch \$*" >> "$CALLS"
[ -n "\${FAKE_FAIL:-}" ] && exit 1
echo "Run,ReleaseDate,LoadDate,spots,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySource,LibrarySelection,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Sample,BioSample,SampleName"
EOS
# runinfo columns used by the module: 1 Run, 2 date, 4 spots, 13 strategy, 16 layout, 19 platform; pad to 30.
echo 'echo "SRR1234567,2024-01-01,x,500000,x,x,x,x,x,x,x,mycelium,RNA-Seq,x,x,PAIRED,x,x,ILLUMINA,x,SRP1,PRJ1,x,x,x,x,x,x,x,sampleA"' >> "$BIN/efetch"
printf '#!/bin/sh\necho "${FAKE_COUNT:-1}"\n' > "$BIN/xtract"   # esearch Count
printf '#!/bin/sh\nexit 0\n' > "$BIN/sleep"       # no 30/60 s retry backoff
command -v timeout >/dev/null || printf '#!/bin/sh\nshift\nexec "$@"\n' > "$BIN/timeout"   # macOS has none
chmod +x "$BIN"/*

pass=0; fail=0
check() { if eval "$2"; then pass=$((pass+1)); echo "ok   $1"; else fail=$((fail+1)); echo "FAIL $1"; fi; }
rows() { grep -c . "$SQ/$TAG.sra_query.csv" 2>/dev/null || echo 0; }
run() {  # run [FAKE_FAIL=1] ; reports number of efetch calls in $NCALL
    : > "$CALLS"
    ( cd "$T/launch" && env PATH="$BIN:$PATH" "$@" nextflow run "$HARNESS" -ansi-log false -w "$T/work" \
        --species_tag "${TAGV:-$TAG}" --taxonid "${TIDV:-$TID}" > "$T/nf.log" 2>&1 )
    NCALL=$(grep -c efetch "$CALLS" 2>/dev/null || true); NCALL=${NCALL:-0}
}
reset() { rm -rf "$T/launch" "$T/work"; mkdir -p "$SQ"; }
HDR='species_tag,taxonid,sra_accession,spots,platform,layout'

reset; run FAKE_OK=1
check "1 hit: CSV has a run row"      '[ "$(rows)" -eq 2 ]'
check "1 hit: key written"            '[ "$(cat "$SQ/$TAG.sra_query.key" 2>/dev/null)" = "v2|$TID" ]'

reset; run FAKE_FAIL=1
check "2 failed query: header-only CSV" '[ "$(rows)" -eq 1 ]'
check "2 failed query: no key"          '[ ! -e "$SQ/$TAG.sra_query.key" ]'

echo "$HDR" > "$SQ/$TAG.sra_query.csv"; rm -rf "$T/work"; run FAKE_OK=1
check "3 header-only, no key: re-queried" '[ "$NCALL" -ge 1 ]'
check "3 ...and the CSV was replaced"     '[ "$(rows)" -eq 2 ]'

reset; echo "$HDR" > "$SQ/$TAG.sra_query.csv"; echo "v2|$TID" > "$SQ/$TAG.sra_query.key"; run FAKE_OK=1
check "4 header-only + matching key: reused" '[ "$NCALL" -eq 0 ]'
check "4 ...CSV unchanged"                   '[ "$(rows)" -eq 1 ]'

reset; echo "$HDR" > "$SQ/$TAG.sra_query.csv"; echo "v2|999" > "$SQ/$TAG.sra_query.key"; run FAKE_OK=1
check "5 header-only + other taxon key: re-queried" '[ "$NCALL" -ge 1 ]'

reset; printf '%s\n%s,%s,SRR9,1000000,ILLUMINA,PAIRED\n' "$HDR" "$TAG" "$TID" > "$SQ/$TAG.sra_query.csv"; run FAKE_OK=1
check "6 CSV with a run: reused" '[ "$NCALL" -eq 0 ]'

reset; run FAKE_COUNT=3
check "7 3 matches: -stop 3"  'grep -q -- "-stop 3\$" "$CALLS"'
reset; run FAKE_COUNT=900
check "8 900 matches: -stop 250" 'grep -q -- "-stop 250\$" "$CALLS"'
reset; run FAKE_COUNT=0
check "9 0 matches: efetch not called" '[ "$NCALL" -eq 0 ]'
check "9 0 matches: header-only CSV + key" '[ "$(rows)" -eq 1 ] && [ "$(cat "$SQ/$TAG.sra_query.key" 2>/dev/null)" = "v2|$TID" ]'

reset; TIDV='162425+227321' run FAKE_COUNT=1
check "10 list: species name in query" 'grep -q "\"Aspergillus nidulans\"\[Organism\] OR txid162425\[Organism\] OR txid227321\[Organism\]" "$CALLS"'
check "10 list: key lists both"        '[ "$(cat "$SQ/$TAG.sra_query.key" 2>/dev/null)" = "v2|162425+227321" ]'
check "10 list: CSV taxonid is the first" 'sed -n 2p "$SQ/$TAG.sra_query.csv" | grep -q "^$TAG,162425,"'
reset; echo "$HDR" > "$SQ/$TAG.sra_query.csv"; echo "v2|162425+227321" > "$SQ/$TAG.sra_query.key"; TIDV='162425+227321' run FAKE_COUNT=1
check "11 list + matching key: reused" '[ "$NCALL" -eq 0 ]'
reset; TAGV=Fusarium_sp. TIDV=5506 run FAKE_COUNT=1
check "12 non-binomial tag: txid only" 'grep -q "(txid5506\[Organism\])" "$CALLS" && ! grep -q "Fusarium sp" "$CALLS"'

echo "passed=$pass failed=$fail"
[ "$fail" -eq 0 ] || { echo "--- last nextflow log ---"; tail -30 "$T/nf.log"; exit 1; }
