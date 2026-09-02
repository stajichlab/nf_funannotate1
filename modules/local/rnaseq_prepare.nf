// Run funannotate train for the representative assembly of a species, producing the
// shared Trinity-GG transcriptome. With params.stop_after_trinity=true (funannotate
// >=1.9) uses `train --stop_after_trinity` -> scratch train stopped after Trinity.
// With stop_after_trinity=false (funannotate 1.8.x, no such flag) runs the FULL
// train, persisted into training_target (the rep is trained once; FUNANNOTATE_TRAIN
// skips it) and still archives the shared trinity.fasta. All non-representative
// strains reuse the archived FASTA in FUNANNOTATE_TRAIN via --trinity.
// storeDir-cached.
// Resources overridden by withName: '.*:RNASEQ_PREPARE' in conf/profile_annotate.config.
process RNASEQ_PREPARE {
    label 'funannotate'
    label 'process_high'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_data"

    input:
    tuple val(species_tag), val(meta), val(genome_fa), path(r1), path(r2), path(se)

    output:
    tuple val(species_tag),
            path("${species_tag}.trinity-GG.fasta"), emit: shared

    script:
    def out           = meta.id
    def species       = meta.species
    def strain        = meta.strain
    def header_length = params.header_length
    // funannotate >= 1.9 supports `train --stop_after_trinity`: assemble the shared
    // Trinity-GG transcriptome for the representative WITHOUT running PASA, letting
    // FUNANNOTATE_TRAIN do the alignment for every strain via --trinity. The 1.8.x
    // line has NO such flag (train.py rejects it as an unrecognized argument), so for
    // older releases we run the FULL train for the representative, persisting it into
    // training_target (so the rep is trained ONCE, not twice — FUNANNOTATE_TRAIN sees
    // the finished train and skips) and still archive its trinity.fasta as the shared
    // transcriptome. Computed here in Groovy (not inside the shell string) so the
    // training dir/flags below interpolate correctly at parse time.
    def stopTrinity = params.stop_after_trinity.toBoolean()
    def OUTDIR      = stopTrinity ? "\$SCRATCH/${out}" : "${params.training_target}/${out}"
    def TRINITY     = stopTrinity ? '--stop_after_trinity --no_trimmomatic' : ''
    def STAGE       = stopTrinity ? 'scratch' : 'persistent (full train, no --stop_after_trinity)'
    """
    # ── Empty-reads sentinel: no RNA-seq found by SRA_FETCH / SRA_FETCH_SE ──
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ]; then
        echo "[INFO] No RNAseq reads for ${species_tag}; writing empty shared markers"
        touch ${species_tag}.trinity-GG.fasta
        exit 0
    fi

    # ── If representative was already trained, just extract shared files ──────
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    if [ -f "\$TRAIN_GFF3" ]; then
        echo "[INFO] Training already complete for ${out}; extracting shared files to rnaseq_data"
        TRAINDIR="${params.training_target}/${out}/training"
        TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
        if [ -n "\$TRINITY_FA" ]; then
            cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
        else
            touch ${species_tag}.trinity-GG.fasta
        fi
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    # Node-local scratch may not exist / be writable if an inherited \$SCRATCH
    # points at another node's path (e.g. the driver's) — fall back to the
    # task workdir so Trinity's OUTDIR and the cleanup below land somewhere
    # real. Any child using \$TMPDIR inherits the same resolved path.
    SCRATCH=\$(printf '%s' "\${SCRATCH:-}" | tr -d '\\n\\r')
    SCRATCH=\${SCRATCH:-/tmp}
    if [ ! -d "\$SCRATCH" ] || [ ! -w "\$SCRATCH" ]; then
        SCRATCH="\$PWD"
    fi
    TMPDIR=\$SCRATCH

    echo "[INFO] RNASEQ_PREPARE: running funannotate train for representative ${out} (species: ${species_tag})"
    if [ "${params.debug.toBoolean()}" = "true" ]; then
        echo "[DEBUG] RNASEQ_PREPARE: stop_after_trinity=${stopTrinity} OUTDIR=${OUTDIR} TRINITY_FLAGS=${TRINITY}"
    fi

    # Inflate a gzipped clean genome to a local uncompressed copy.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    if [ -s "${r1}" ]; then
        echo "[INFO] RNASEQ_PREPARE: funannotate train (PE, ${STAGE}) for representative ${out}"
        funannotate train -i "\$GENOME_IN" -o ${OUTDIR} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            ${TRINITY}
    else
        echo "[INFO] RNASEQ_PREPARE: funannotate train (SE, ${STAGE}) for representative ${out}"
        funannotate train -i "\$GENOME_IN" -o ${OUTDIR} \\
            --single_norm ${se} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            ${TRINITY}
    fi

    # ── Copy shared outputs to rnaseq_data/ ──────────────────────────────────
    TRAINDIR="${OUTDIR}/training"
    TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
    if [ -n "\$TRINITY_FA" ]; then
        cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
    else
        echo "[WARN] No trinity.fasta found under \$TRAINDIR for ${out}"
        touch ${species_tag}.trinity-GG.fasta
    fi

    # Scratch (stop_after_trinity) runs clean up; full persistent runs keep the
    # training dir for FUNANNOTATE_PREDICT (via symlink) and so FUNANNOTATE_TRAIN
    # can skip the already-trained representative.
    if [ "${stopTrinity}" = "true" ]; then
        rm -rf "\$SCRATCH/${out}"
    fi
    echo "[INFO] RNASEQ_PREPARE complete for ${species_tag}"
    """

    stub:
    def out = meta.id
    """
    echo ">stub_trinity_${species_tag}" > ${species_tag}.trinity-GG.fasta
    mkdir -p ${params.training_target}/${out}/training
    touch ${params.training_target}/${out}/training/funannotate_train.pasa.gff3
    """
}
