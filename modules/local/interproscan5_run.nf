process INTERPROSCAN5_RUN {
    label 'interproscan'
    label 'process_medium'
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}/annotate_misc/iprscan.xml"), emit: results
    path 'versions.yml',                                            emit: versions

    script:
    def out      = meta.id
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    mkdir -p ${out}/annotate_misc
    interproscan.sh -i ${proteins} -f XML -o ${out}/annotate_misc/iprscan.xml \\
        -dp -goterms -pa -t p -cpu ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(interproscan.sh --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """

    stub:
    def out = meta.id
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/iprscan.xml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: 5.65-97.0
    END_VERSIONS
    """
}
