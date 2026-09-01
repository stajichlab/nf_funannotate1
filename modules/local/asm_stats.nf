process ASM_STATS {
    label 'setup'
    label 'process_low'

    storeDir { params.tables_dir }

    input:
    path samples
    path genome_dir

    output:
    path 'asm_stats.tsv.gz', emit: stats
    path 'versions.yml',     emit: versions

    script:
    """
    set -euo pipefail

    TMPFILE=\$(mktemp)
    trap 'rm -f \$TMPFILE' EXIT

    printf 'ASMID\\ttotal_length_bp\\tN50_bp\\tcontig_count\\n' > \$TMPFILE

    awk -F',' 'NR>1 {print \$3}' ${samples} | sort -u | while read asmid; do
        [ -z "\$asmid" ] && continue
        asmid="\$(echo "\$asmid" | xargs)"

        if [ -f "${genome_dir}/\${asmid}.fa.gz" ]; then
            genome="${genome_dir}/\${asmid}.fa.gz"
        elif [ -f "${genome_dir}/\${asmid}.fa" ]; then
            genome="${genome_dir}/\${asmid}.fa"
        elif [ -f "${genome_dir}/\${asmid}.masked.fasta.gz" ]; then
            genome="${genome_dir}/\${asmid}.masked.fasta.gz"
        elif [ -f "${genome_dir}/\${asmid}.masked.fasta" ]; then
            genome="${genome_dir}/\${asmid}.masked.fasta"
        else
            echo "[WARN] No genome file found for \${asmid} in ${genome_dir}" >&2
            continue
        fi

        lens=\$(gzip -cd -f "\$genome" | awk '
            /^>/     { if (seqlen) print seqlen; seqlen=0; next }
            { seqlen += length(\$0) }
            END      { if (seqlen) print seqlen }
        ' | sort -rn)

        total_bp=\$(echo "\$lens" | awk '{s+=\$1} END{print s+0}')
        contigs=\$(echo "\$lens" | wc -l)
        n50=\$(echo "\$lens" | awk -v total="\$total_bp" 'BEGIN{sum=0} {sum+=\$1; if (sum >= total/2) {print \$1; exit}}')

        [ -z "\$total_bp" ] && total_bp="0"
        [ -z "\$n50" ] && n50="0"
        [ -z "\$contigs" ] && contigs="0"

        printf '%s\\t%s\\t%s\\t%s\\n' "\$asmid" "\$total_bp" "\$n50" "\$contigs" >> \$TMPFILE
    done

    gzip -c \$TMPFILE > asm_stats.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
        gzip: \$(gzip --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf 'ASMID\\ttotal_length_bp\\tN50_bp\\tcontig_count\\n' | gzip -c > asm_stats.tsv.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: 1.3.4
    END_VERSIONS
    """
}
