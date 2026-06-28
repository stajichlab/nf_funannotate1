/*
 * asm_stats — Generate assembly statistics for clean genomes
 *
 * This module generates asm_stats.tsv with columns: ASMID, total_length_bp, N50_bp, contig_count.
 * Stats are used by earlgrey_mask.nf to select representative genomes per species (SELECT_REPS).
 *
 * Include in your workflow:
 *   include { ASM_STATS } from './modules/asm_stats'
 *   ASM_STATS(samples_csv, genome_dir)
 */

process ASM_STATS {
    label 'setup'

    storeDir { params.tables_dir }

    cpus   4
    memory '8 GB'
    time   '2h'

    input:
    path samples
    path genome_dir

    output:
    path 'asm_stats.tsv.gz', emit: stats

    script:
    """
    set -euo pipefail

    TMPFILE=\$(mktemp)
    trap 'rm -f \$TMPFILE' EXIT

    printf 'ASMID\\ttotal_length_bp\\tN50_bp\\tcontig_count\\n' > \$TMPFILE

    # Extract ASMIDs from samples.csv
    awk -F',' 'NR>1 {print \$2}' ${samples} | sort -u | while read asmid; do
        [ -z "\$asmid" ] && continue
        asmid="\$(echo "\$asmid" | xargs)"  # trim whitespace

        # Look for genome file: prefer .fa.gz, fall back to .fa, then .masked.fasta.gz
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

        # Use seqkit to compute stats
        total_bp=\$(seqkit stats -T "\$genome" 2>/dev/null | tail -n 1 | awk '{print \$4}')
        n50=\$(seqkit fx2tab -l "\$genome" 2>/dev/null | sort -rn -k2 | \\
            awk -v total="\$total_bp" 'BEGIN{sum=0} {sum+=\$2; if(sum >= total/2) {print \$2; exit}}')
        contigs=\$(seqkit stats -T "\$genome" 2>/dev/null | tail -n 1 | awk '{print \$3}')

        [ -z "\$total_bp" ] && total_bp="0"
        [ -z "\$n50" ] && n50="0"
        [ -z "\$contigs" ] && contigs="0"

        printf '%s\\t%s\\t%s\\t%s\\n' "\$asmid" "\$total_bp" "\$n50" "\$contigs" >> \$TMPFILE
    done

    pigz -c \$TMPFILE > asm_stats.tsv.gz
    echo "[INFO] Assembly statistics written: asm_stats.tsv.gz"
    """

    stub:
    """
    printf 'ASMID\\ttotal_length_bp\\tN50_bp\\tcontig_count\\n' | pigz -c > asm_stats.tsv.gz
    echo "[STUB] ASM_STATS"
    """
}
