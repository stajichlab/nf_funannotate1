// Run DeepTMHMM on one assembly's proteome, producing annotate_misc/TMRs.gff3
// for funannotate annotate's --tmhmm arg (overrides Phobius's TM calls).
//
// DeepTMHMM has no UCR HPCC Lmod module — it only ships as the
// container built from ~/projects/containers/deeptmhmm_container/deeptmhmm.def
// (params.container_deeptmhmm), so it always runs under -profile singularity.
// Queue/clusterOptions GPU routing (params.deeptmhmm_gpu) lives in
// conf/provision_singularity.config's withName block, not the usual
// conf/provision_ucr_hpcc.config location, since that file is skipped under
// -profile singularity.
process DEEPTMHMM_ANNOTATION {
    label 'deeptmhmm'
    label 'process_medium'
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}/annotate_misc/TMRs.gff3"), emit: results
    path 'versions.yml',                                          emit: versions

    script:
    def out      = meta.id
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    TASKDIR=\$PWD
    mkdir -p ${out}/annotate_misc

    # /opt/deeptmhmm (baked into the image at build time, see
    # conf/provision_singularity.config) holds predict.py plus its model
    # weights, loaded via paths relative to cwd -- cwd must be /opt/deeptmhmm
    # for the run itself, hence the cd. --output-dir writes results directly
    # (no copy-back dance needed, unlike the old biswasaneel image's predict.py).
    cd /opt/deeptmhmm
    python3 predict.py --fasta ${proteins} --output-dir \${TASKDIR}/${out}/annotate_misc
    cd \${TASKDIR}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptmhmm: 1.0-academic
    END_VERSIONS
    """

    stub:
    def out = meta.id
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/TMRs.gff3
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptmhmm: stub
    END_VERSIONS
    """
}
