// Completeness of funannotate predict's FINAL delivered gene set.
//
// funannotate's own logfiles/busco.log only captures INTERMEDIATE ab-initio
// training checkpoints -- BUSCO run in genome mode BEFORE any gene prediction
// (to assess assembly completeness / pick augustus training loci), and BUSCO
// run in protein mode against augustus's raw training predictions
// specifically. Neither reflects the completeness of the FINAL EVM/PASA-
// merged gene set actually delivered in predict_results/, so this reruns
// BUSCO directly against predict_results/<id>.proteins.fa (id = SPECIES_STRAIN tag).
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
// redoes it. Writes straight into genome_annotation/<id>/busco_completeness/
// (the same <id> directory as predict_results/; was <asmid> before 2026-09-26),
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

    storeDir "${params.target}/${meta.id}/busco_completeness"

    input:
    tuple val(meta), path(proteins_fa)

    output:
    path("${meta.id}")

    script:
    """
    # params.busco_lineages may be a BUSCO --download_path root (lineage dirs under
    # <root>/lineages/) or the lineage dirs' own directory, as in BUSCO_GENOME. For
    # the latter, give BUSCO a local download_path whose lineages/ entry links to it.
    # (A lineages-level value passed straight to --download_path doubled the level:
    # ".../v10/lineages/lineages/eurotiomycetes_odb10 does not exist", 2026-10-01.)
    BUSCO_DL="${params.busco_lineages}"
    if [ ! -d "\$BUSCO_DL/lineages/${meta.busco}" ] && [ -d "\$BUSCO_DL/${meta.busco}" ]; then
        BUSCO_DL="\$PWD/busco_download"
        mkdir -p "\$BUSCO_DL/lineages"
        ln -sfn "${params.busco_lineages}/${meta.busco}" "\$BUSCO_DL/lineages/${meta.busco}"
    fi
    busco -i ${proteins_fa} -l ${meta.busco} \\
        -m proteins -c ${task.cpus} -o ${meta.id} --out_path . -f \\
        --offline --download_path "\$BUSCO_DL"

    # hmmer_output/ and busco_sequences/ hold ~99% of the ~3,500 small files per genome (3,567 of
    # 3,576 in a measured result) and nothing downstream reads them; the summaries, full_table.tsv
    # and missing_busco_list.tsv stay as plain files. Keep the two directories as one tarball so
    # the storeDir move is a handful of files, not thousands (a retry colliding with a half-moved
    # directory failed with "Directory not empty", 2026-10-02).
    RUN_DIR=\$(ls -d ${meta.id}/run_* 2>/dev/null | head -1 || true)
    if [ -n "\$RUN_DIR" ]; then
        BULK=""
        for d in hmmer_output busco_sequences; do [ -d "\$RUN_DIR/\$d" ] && BULK="\$BULK \$d"; done
        if [ -n "\$BULK" ]; then
            tar -C "\$RUN_DIR" -czf "\$RUN_DIR/busco_run_dirs.tar.gz" \$BULK
            for d in \$BULK; do
                find "\$RUN_DIR/\$d" -type f -delete
                find "\$RUN_DIR/\$d" -depth -type d -delete
            done
        fi
    fi
    """

    stub:
    """
    mkdir -p ${meta.id}
    printf '# BUSCO version is: 6.1.0\\n# The lineage dataset is: ${meta.busco}\\n\\tC:99.0%%[S:98.0%%,D:1.0%%],F:0.5%%,M:0.5%%,n:758\\n' \\
        > "${meta.id}/short_summary.specific.${meta.busco}.${meta.id}.txt"
    """
}
