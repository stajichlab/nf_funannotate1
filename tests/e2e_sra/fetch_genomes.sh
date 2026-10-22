#!/usr/bin/env bash
# Download the two test genomes (NCBI RefSeq) into ./genomes/, named as samples.csv expects.
# Usage: bash fetch_genomes.sh        (from the launch directory holding samples.csv)
set -euo pipefail
mkdir -p genomes
get() {  # <GCF accession, e.g. GCF_000002525.2> <output name>
    local acc="$1" out="$2" p dir
    p="${acc:0:3}/${acc:4:3}/${acc:7:3}/${acc:10:3}"
    dir=$(curl -fsSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/${p}/" \
            | grep -o "href=\"${acc}_[^\"/]*/\"" | head -1 | sed 's/href="//; s/\/"//')
    [ -n "$dir" ] || { echo "no FTP directory for $acc" >&2; exit 1; }
    curl -fsSL -o "genomes/${out}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/${p}/${dir}/${dir}_genomic.fna.gz"
    gzip -t "genomes/${out}"
    echo "genomes/${out}"
}
get GCF_000002525.2 Yarrowia_lipolytica_CLIB122.fna.gz
get GCF_000091045.1 Cryptococcus_deneoformans_JEC21.fna.gz
