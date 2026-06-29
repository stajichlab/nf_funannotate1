// Soft-mask each assembly using funannotate mask with tantan.
// storeDir caches the masked FASTA alongside the clean genome at
// input_clean_genomes/<asmid>.masked.fasta.gz, so this runs at most once per assembly.
// Downstream tools inflate it on the fly.
// Resources overridden per-profile by withName: '.*:MASKREPEAT_TANTAN_RUN' in
// conf/profile_annotate.config (short queue, 8 cpus, 16 GB, 2h).
process MASKREPEAT_TANTAN_RUN {
    label 'funannotate'
    label 'process_medium'
    tag "${meta.id}"

    storeDir "${launchDir}/input_clean_genomes"

    input:
    tuple val(meta), val(genome_fa)

    output:
    tuple val(meta), path("${meta.asmid}.masked.fasta.gz"), emit: masked

    script:
    def asmid = meta.asmid
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac
    funannotate mask -i "\$GENOME_IN" -o ${asmid}.masked.fasta -m tantan --cpus ${task.cpus}
    # Deliver the soft-masked genome gzip-compressed to save space; consumers inflate it.
    pigz -f ${asmid}.masked.fasta
    """

    stub:
    def asmid = meta.asmid
    """
    echo ">stub_${asmid}_masked" | pigz -c > ${asmid}.masked.fasta.gz
    """
}
