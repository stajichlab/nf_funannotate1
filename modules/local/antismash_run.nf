process ANTISMASH_RUN {
    label 'antismash'
    label 'process_medium'
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}.antismash.done"), emit: results
    path 'versions.yml',                                emit: versions

    script:
    def out    = meta.id
    def gbk    = "${params.target}/${out}/predict_results/${out}.gbk"
    // Option B persistence (like FUNANNOTATE_PREDICT): write directly to the
    // persistent target dir, not a task-workdir-relative path -- annotate_genome.nf's
    // done-check and FUNANNOTATE_ANNOTATE both look for output under
    // params.target/${out}/antismash_local, which a relative `${out}/antismash_local`
    // here would never actually reach (confirmed: real antismash output sat
    // unpublished in the task work dir after a successful run).
    def outdir = "${params.target}/${out}/antismash_local"
    """
    GBK="${gbk}"
    if [ ! -f "\$GBK" ] && [ -f "${gbk}.gz" ]; then
        zcat "${gbk}.gz" > ${out}.predict.gbk
        GBK=${out}.predict.gbk
    fi
    if [ ! -f "\$GBK" ]; then
        echo "ERROR: predict GBK not found: ${gbk}[.gz]" >&2
        exit 1
    fi
    rm -rf ${outdir}
    mkdir -p ${outdir}
    antismash --taxon ${params.antismash_taxon} \\
        --databases ${params.antismash_db} \\
        --output-dir ${outdir} \\
        --genefinding-tool none \\
        --fullhmmer --clusterhmmer --cb-general --pfam2go \\
        -c ${task.cpus} \\
        \$GBK
    gzip ${outdir}/*.json
    touch ${out}.antismash.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version 2>&1 | grep -oP '(?<=antiSMASH )\\S+' || antismash --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    def out    = meta.id
    def outdir = "${params.target}/${out}/antismash_local"
    """
    mkdir -p ${outdir}
    touch ${outdir}/${out}.json.gz
    touch ${outdir}/index.html
    touch ${out}.antismash.done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: 8.0.4
    END_VERSIONS
    """
}
