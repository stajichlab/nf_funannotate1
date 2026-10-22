#!/usr/bin/env bash
# Check the outputs of an end-to-end SRA run (run it in the launch directory).
#   bash check.sh             after the first run: query, reads, Trinity, annotation exist
#   bash check.sh rerun       after a second launch: the SRA steps were reused, not repeated
#   bash check.sh own-reads   after the own-reads launch (README.md step 4): no query was made
#   bash check.sh key-removed after deleting the key too: the species was searched again
set -uo pipefail
MODE="${1:-run1}"; LOGDIR="logs/nextflow"
SPECIES="Yarrowia_lipolytica Cryptococcus_deneoformans"
pass=0; fail=0
check() { if eval "$2"; then pass=$((pass+1)); echo "ok   $1"; else fail=$((fail+1)); echo "FAIL $1"; fi; }
latest_log() { ls -t nf_run*.log "$LOGDIR"/*.log .nextflow.log 2>/dev/null | head -1; }

case "$MODE" in
run1)
    for sp in $SPECIES; do
        q="rnaseq_reads/sra_query/$sp.sra_query"
        check "$sp: query CSV lists a run"       '[ "$(grep -c . "'"$q"'.csv" 2>/dev/null)" -ge 2 ]'
        check "$sp: query key written"           '[ -s "'"$q"'.key" ]'
        check "$sp: reads R1 non-empty"          '[ -s "rnaseq_reads/'"$sp"'_norm_R1.fastq.gz" ]'
        check "$sp: reads R2 non-empty"          '[ -s "rnaseq_reads/'"$sp"'_norm_R2.fastq.gz" ]'
        check "$sp: Trinity transcripts exist"   '[ -s "rnaseq_data/'"$sp"'.trinity-GG.fasta" ]'
    done
    n=$(find genome_annotation -name '*.gbk' -size +0c 2>/dev/null | wc -l)
    check "GBK per genome (found $n, want 2)"    '[ "$n" -ge 2 ]'
    find rnaseq_reads -maxdepth 1 -name '*_norm_*' -exec stat -c '%n %Y %s' {} + | sort > .e2e_reads.snapshot
    echo "--- chosen runs ---"; for sp in $SPECIES; do printf '%-28s %s\n' "$sp" "$(sed 1d rnaseq_reads/sra_query/$sp.sra_query.csv | cut -d, -f3 | paste -sd' ' -)"; done ;;
rerun)
    # Run again after `rm -rf work .nextflow`, so nothing is resumed: the SRA outputs must come
    # from the caches in rnaseq_reads/ (storeDir) and rnaseq_reads/sra_query/ (query CSV + key).
    L=$(latest_log); echo "log: $L"
    cat work/*/*/.command.log 2>/dev/null > /tmp/e2e_cmdlog.$$
    for sp in $SPECIES; do
        check "$sp: query reused, not repeated" 'grep -q "Reusing cached result for '"$sp"'" /tmp/e2e_cmdlog.$$ && ! grep -q "Found .* paired-end accessions for '"$sp"'" /tmp/e2e_cmdlog.$$'
    done
    check "read files untouched (same size and mtime as after run 1)" 'find rnaseq_reads -maxdepth 1 -name "*_norm_*" -exec stat -c "%n %Y %s" {} + | sort | cmp -s - .e2e_reads.snapshot'
    check "no download ran (no fasterq/prefetch)" '! grep -qiE "prefetch|fasterq-dump|fastq-dump" /tmp/e2e_cmdlog.$$'
    rm -f /tmp/e2e_cmdlog.$$ ;;
own-reads)
    cat work/*/*/.command.log 2>/dev/null > /tmp/e2e_cmdlog.$$
    check "Yarrowia query reused (header-only CSV + key)" 'grep -q "Reusing cached result for Yarrowia_lipolytica" /tmp/e2e_cmdlog.$$'
    check "no NCBI query for Yarrowia"                    '! grep -q "Found .* paired-end accessions for Yarrowia_lipolytica" /tmp/e2e_cmdlog.$$'
    rm -f /tmp/e2e_cmdlog.$$ ;;
key-removed)
    f=rnaseq_reads/sra_query/Yarrowia_lipolytica.sra_query
    cat work/*/*/.command.log 2>/dev/null > /tmp/e2e_cmdlog.$$
    check "Yarrowia searched again (no key)"   'grep -q "Found .* paired-end accessions for Yarrowia_lipolytica" /tmp/e2e_cmdlog.$$'
    check "CSV lists runs again"               '[ "$(grep -c . "'"$f"'.csv")" -ge 2 ]'
    check "key written again"                  '[ -s "'"$f"'.key" ]'
    rm -f /tmp/e2e_cmdlog.$$ ;;
esac
echo "passed=$pass failed=$fail"; [ "$fail" -eq 0 ]
