// Completeness of funannotate predict's FINAL delivered gene set.
//
// funannotate's own logfiles/busco.log only captures INTERMEDIATE ab-initio
// training checkpoints -- BUSCO run in genome mode BEFORE any gene prediction
// (to assess assembly completeness / pick augustus training loci), and BUSCO
// run in protein mode against augustus's raw training predictions
// specifically. Neither reflects the completeness of the FINAL EVM/PASA-
// merged gene set actually delivered in predict_results/, so this reruns
// BUSCO directly against predict_results/<asmid>.proteins.fa.
//
// Always containerized (params.container_busco_completeness), regardless of
// which provisioning axis (conda/pixi/ucr_hpcc/singularity) this cell uses --
// see conf/profile_annotate.config's apptainer.enabled/withLabel block --
// so BUSCO's own version/environment is held constant across all 6 benchmark
// cells and never becomes a confound in the conda-vs-container /
// v1.8.17-vs-v1.9.0-beta10 comparison itself.
//
// storeDir-cached like MASKREPEAT_TANTAN_RUN: a pure function of the protein
// set + lineage, keyed on the delivered output path, so a `-resume` never
// redoes it. Writes straight into genome_annotation/<asmid>/busco_completeness/,
// which Funannotate_benchmarking's scripts/lib.py find_busco_summary() (a
// recursive walk for short_summary*.txt) already picks up with no changes
// needed downstream in compare_predictions.py / validate_harness.py.
//
// Resources overridden by withName: '.*:BUSCO_COMPLETENESS' in
// conf/profile_annotate.config; queue routing in conf/provision_ucr_hpcc.config.
process BUSCO_COMPLETENESS {
    label 'busco_completeness'
    label 'process_medium'
    tag "${meta.id}"

    storeDir "${launchDir}/genome_annotation/${meta.asmid}/busco_completeness"

    input:
    tuple val(meta), path(proteins_fa)

    output:
    path("${meta.asmid}")

    script:
    """
    busco -i ${proteins_fa} -l ${params.busco_lineages}/${meta.busco} \\
        -m proteins -c ${task.cpus} -o ${meta.asmid} --out_path . -f
    """

    stub:
    """
    mkdir -p ${meta.asmid}
    printf '# BUSCO version is: 6.1.0\\n# The lineage dataset is: ${meta.busco}\\n\\tC:99.0%%[S:98.0%%,D:1.0%%],F:0.5%%,M:0.5%%,n:758\\n' \\
        > "${meta.asmid}/short_summary.specific.${meta.busco}.${meta.asmid}.txt"
    """
}
