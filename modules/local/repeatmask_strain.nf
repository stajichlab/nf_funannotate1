// Apply the species curated TE library (from EARLGREY_BUILD_LIB) to a conspecific
// strain with RepeatMasker (-xsmall = soft-masking). storeDir caches under
// params.earlgrey_outdir/<sp_safe>/strains/ — separate from input_clean_genomes/.
// Resources: label 'repeatmask' in conf/profile_earlgrey.config (32 cpus, 32 GB, 8 h).
process REPEATMASK_STRAIN {
    label 'repeatmask'
    label 'process_medium'
    tag  "${asmid} (${species})"

    storeDir { "${params.earlgrey_outdir}/${species.replaceAll(/[^A-Za-z0-9._-]+/, '_')}/strains" }

    input:
    tuple val(asmid), val(species), path(genome), path(library)

    output:
    tuple val(asmid), path("${asmid}.masked.fasta.gz")

    script:
    """
    GENOME_IN="${genome}"
    case "${genome}" in
        *.gz) echo "[INFO] Inflating ${genome}"; gzip -dc "${genome}" > genome_input.fa; GENOME_IN=genome_input.fa ;;
    esac

    mkdir -p rmask_out
    RepeatMasker \\
        -lib ${library} \\
        -xsmall \\
        -pa ${task.cpus} \\
        -dir rmask_out \\
        "\$GENOME_IN"

    MASKED="rmask_out/\$(basename "\$GENOME_IN").masked"
    if [ -f "\$MASKED" ]; then
        gzip -c "\$MASKED" > ${asmid}.masked.fasta.gz
    else
        echo "WARN: no repeats masked for ${asmid}; using unmasked genome" >&2
        gzip -c "\$GENOME_IN" > ${asmid}.masked.fasta.gz
    fi
    """

    stub:
    """
    printf '>stub_${asmid}\\nacgtACGTacgt\\n' | gzip -c > ${asmid}.masked.fasta.gz
    """
}
