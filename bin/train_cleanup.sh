#!/usr/bin/env bash
# train_cleanup.sh TRAINDIR KEEP_PASA
#
# Remove funannotate-train intermediates that funannotate predict never reads,
# once a genome's training has resolved (funannotate_train.pasa.gff3 exists).
# Called by FUNANNOTATE_TRAIN after a successful train and by
# FUNANNOTATE_PREDICT after a successful predict (which covers genomes trained
# before params.train_cleanup was turned on). Idempotent.
#
# predict reads only training/funannotate_train.* (pasa.gff3, transcripts.gff3,
# coordSorted.bam, stringtie.gtf, trinity-GG.fasta); the pipeline also reads
# transcript.alignments.bam. Several of those are symlinks into this directory
# (funannotate_train.trinity-GG.fasta -> trinity.fasta, ...coordSorted.bam ->
# trinity.alignments.bam), so anything a symlink here points at is kept.
#
# pasa/ (PASA assemblies, TransDecoder dirs, the PASA sqlite DB) is only
# removed with KEEP_PASA=0: funannotate update reuses it. Measured on Bd
# (2026-09-28): ~1.5 of 1.7 GB per genome with pasa/, ~0.6 GB without.
set -uo pipefail

d="$1"
keep_pasa="${2:-1}"
[ -d "$d" ] || exit 0
if [ ! -s "$d/funannotate_train.pasa.gff3" ]; then
    echo "[INFO] train_cleanup: no funannotate_train.pasa.gff3 in $d; left as is"
    exit 0
fi

keep=()
for l in "$d"/* "$d"/.[!.]*; do
    [ -L "$l" ] && keep+=("$(readlink -f "$l")")
done

targets=(getBestModel genome.fasta.gmap genome.fasta.cidx trinity.fasta.clean
         trinity.fasta.clean.cidx trinity.fasta.clean.fai trinity.fasta.cln
         trinity.fasta.cidx outparts_cln.sort pasa.step1.gff3)
[ "$keep_pasa" = 0 ] && targets+=(pasa)

freed=0
for t in "${targets[@]}"; do
    p="$d/$t"
    { [ -e "$p" ] && [ ! -L "$p" ]; } || continue
    rp=$(readlink -f "$p")
    held=""
    for k in "${keep[@]}"; do
        case "$k" in "$rp" | "$rp"/*) held="$k"; break ;; esac
    done
    if [ -n "$held" ]; then
        echo "[INFO] train_cleanup: keeping $t (a symlink points at $held)"
        continue
    fi
    kb=$(du -sk "$p" 2>/dev/null | cut -f1)
    rm -rf "$p" && freed=$((freed + ${kb:-0}))
done
echo "[INFO] train_cleanup: freed $((freed / 1024)) MB in $d"
