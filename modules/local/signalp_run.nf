process SIGNALP_RUN {
    label 'signalp'
    label 'process_medium'
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}/annotate_misc/signalp.results.txt"), emit: results
    path 'versions.yml',                                                    emit: versions

    script:
    def out      = meta.id
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    TMPDIR=\${SCRATCH:-/tmp}
    MODEL_DIR_ARG=()

    if [ "${params.signalp_gpu}" = "true" ]; then
        # The baked-in signalp6-fast.sif checkpoint is saved for CPU; GPU
        # inference needs a converted copy of the weights (signalp6 reads
        # its target device off the loaded checkpoint, not from --nv). The
        # conversion itself needs a real GPU to run, so it happens here
        # rather than at image-build time. Cached at a persistent, shared
        # location so it only runs once across all pipeline runs/genomes.
        GPU_MODELS="${params.signalp_gpu_models}"
        if [ ! -f "\$GPU_MODELS/.converted" ]; then
            echo "[INFO] SIGNALP_RUN ${out}: no cached GPU model conversion at \$GPU_MODELS -- converting now"
            SIGNALP_DIR=\$(python3 -c "import signalp, os; print(os.path.dirname(signalp.__file__))")
            CONVERT_TMP="\${TMPDIR}/signalp_gpu_convert.\$\$"
            mkdir -p "\$CONVERT_TMP"
            cp -r "\$SIGNALP_DIR/model_weights/"* "\$CONVERT_TMP/"
            signalp6_convert_models gpu "\$CONVERT_TMP"
            mkdir -p "\$(dirname "\$GPU_MODELS")"
            # Race-safe against concurrent tasks converting at the same time:
            # stage under a unique name, then atomically claim the final path.
            if mv -T "\$CONVERT_TMP" "\$GPU_MODELS" 2>/dev/null; then
                touch "\$GPU_MODELS/.converted"
            else
                echo "[INFO] SIGNALP_RUN ${out}: \$GPU_MODELS already claimed by a concurrent task; discarding our copy"
                rm -rf "\$CONVERT_TMP"
            fi
        fi
        MODEL_DIR_ARG=(--model_dir "\$GPU_MODELS")
    fi

    signalp6 -od \$TMPDIR/${out}_signalp \\
        -org euk --mode fast -format txt \\
        -fasta ${proteins} \\
        --write_procs ${task.cpus} -bs 16 "\${MODEL_DIR_ARG[@]}"
    mkdir -p ${out}/annotate_misc
    cp \$TMPDIR/${out}_signalp/prediction_results.txt ${out}/annotate_misc/signalp.results.txt
    rm -rf \$TMPDIR/${out}_signalp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: \$(signalp6 --version 2>&1 | grep -oP '\\d+\\.\\d+\\S*' | head -1)
    END_VERSIONS
    """

    stub:
    def out = meta.id
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/signalp.results.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: 6.0g
    END_VERSIONS
    """
}
