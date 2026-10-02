// Run funannotate predict for one assembly. Writes directly into params.target/<id>/
// (Option B persistence: no publishDir copy). Emits a small marker file to carry the
// DAG edge without transferring the full predict tree through Nextflow's work/ directory.
// Resources overridden by withName: '.*:FUNANNOTATE_PREDICT' in conf/profile_annotate.config.
process FUNANNOTATE_PREDICT {
    label 'funannotate'
    label 'process_high'
    tag "${meta.id}"

    input:
    // train_fp: fingerprint of the training evidence this predict will consume
    // (FunannotateUtils.trainingFingerprint). It is NOT used by the script -- it
    // exists solely so the training evidence participates in Nextflow's task
    // hash. Without it, re-training a genome leaves the predict hash unchanged
    // and -resume serves the OLD annotation. See the helper's comment for the
    // three occurrences that motivated this.
    tuple val(meta), val(genome_fa), val(genemark_gtf), val(other_gff), val(train_fp)

    output:
    val meta, emit: metadata
    path("${meta.id}.predict.done"), emit: done

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    // Reuse siblings (FUNANNOTATE_PREDICT_SIB) also go stale when the species'
    // shared ab-initio store is refreshed; see the skip check below.
    def shared_params_json = (task.process.endsWith('_SIB') && params.gene_prediction_shared_abinitio) ?
        "${params.gene_prediction_shared_abinitio}/${species.replaceAll(/\s+/, '_')}/parameters.json" : ''
    def header_length = params.header_length
    def transl_table  = meta.transl_table
    // GeneMark GTF supplied by the standalone GENEMARK_RUN step; empty string
    // means "let funannotate run GeneMark internally (or auto-skip it)".
    // Checked by actual file size, not by Groovy truthiness on genemark_gtf
    // itself: GENEMARK_RUN's too-small-genome skip path emits a real but
    // deliberately empty ${out}.genemark.gtf (see genemark_run.nf), and a
    // non-null Path is always truthy in Groovy regardless of its size --
    // without this check every skip would still pass --genemark_gtf
    // <0-byte file> to funannotate.
    def genemark_gtf_file = genemark_gtf ? file(genemark_gtf as String) : null
    def genemark_gtf_ok   = genemark_gtf_file && genemark_gtf_file.exists() && genemark_gtf_file.size() > 0
    def genemark_cli      = genemark_gtf_ok ? "--genemark_gtf ${genemark_gtf}" : "--auto-skip-genemark"
    // -w genemark:1 MUST be passed explicitly whenever genemark_gtf is used,
    // in the SAME -w group as codingquarry:0/glimmerhmm:0 (funannotate's
    // argparse -w/--weights is nargs='+' without action='append', so a
    // second -w occurrence replaces the whole list rather than merging it --
    // ported from BFD's FUNANNOTATE_PREDICT/main.nf, which hit this).
    // Needed because none of this pipeline's funannotate provisioning modes
    // (module/pixi/singularity funannotate env -- see conf/provision_*.config)
    // put gmes_petap.pl on PATH inside this task's environment: GeneMark now
    // runs standalone in GENEMARK_RUN, so funannotate predict's own
    // genemarkcheck is always False here, and predict.py unconditionally
    // zeroes StartWeights["genemark"] when that's the case -- with NO check
    // for whether --genemark_gtf was supplied as an alternative.
    def weight_args = genemark_gtf_ok ? 'codingquarry:0 glimmerhmm:0 genemark:1' : 'codingquarry:0 glimmerhmm:0'
    // Optional external gene-model pass-through (PRODIGAL_RUN output). Passed as
    // --other_gff <gff>:<weight>; funannotate renames the source to
    // 'other_pred1', EVM-validates it, and applies the weight as
    // StartWeights['other_pred1'] (see predict.py --other_gff parsing). Empty
    // string -> omitted entirely.
    // PRODIGAL_RUN now emits gene/mRNA/CDS blocks (prodigal_gff_hier.py): EVM
    // assembles consensus models from gene/mRNA blocks, so a CDS-only source
    // is structurally inert no matter its weight (OC4 gene Sn 0.584 -> 0.839).
    def other_gff_file = other_gff ? file(other_gff as String) : null
    def other_gff_ok   = other_gff_file && other_gff_file.exists() && other_gff_file.size() > 0
    def other_gff_cli  = other_gff_ok ? "--other_gff ${other_gff}:${params.prodigal_weight}" : ''
    """
    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    # Node-local scratch may not exist / be writable if an inherited \$SCRATCH
    # points at another node's path — fall back to the task workdir.
    TMPDIR=\$(printf '%s' "\${SCRATCH:-}" | tr -d '\\n\\r')
    TMPDIR=\${TMPDIR:-/tmp}
    if [ ! -d "\$TMPDIR" ] || [ ! -w "\$TMPDIR" ]; then
        TMPDIR="\$PWD"
    fi

    PREDICTDIR="${params.target}/${out}"
    PREDICT_GBK="\$PREDICTDIR/predict_results/${out}.gbk"
    # RUNDIR is where funannotate predict actually writes (-o). Normally the
    # persistent PREDICTDIR itself. With params.predict_local_scratch it is a
    # node-local dir under \$TMPDIR instead: BUSCO/Augustus training there
    # creates thousands of small files, which crawls on a network filesystem
    # (CephFS: ~1 h blocked on metadata for one 12 Mb genome). Only the final,
    # pruned tree is copied back to PREDICTDIR (see sync_back below).
    RUNDIR="\$PREDICTDIR"
    if [ "${params.predict_local_scratch.toBoolean()}" = "true" ]; then
        RUNDIR="\$TMPDIR/funannotate_predict_${out}"
    fi
    RUN_GBK="\$RUNDIR/predict_results/${out}.gbk"
    # Local-scratch mode: copy logs back on failure (for debugging), the whole
    # pruned tree on success. No-ops when RUNDIR is PREDICTDIR.
    copy_logs_back() {
        if [ "\$RUNDIR" != "\$PREDICTDIR" ]; then
            if [ -d "\$RUNDIR/logfiles" ]; then
                mkdir -p "\$PREDICTDIR" && cp -a "\$RUNDIR/logfiles" "\$PREDICTDIR/"
            fi
            rm -rf "\$RUNDIR"
        fi
    }
    sync_back() {
        if [ "\$RUNDIR" != "\$PREDICTDIR" ]; then
            rm -rf "\$PREDICTDIR/predict_results" "\$PREDICTDIR/predict_misc"
            mkdir -p "\$PREDICTDIR" && cp -a "\$RUNDIR/." "\$PREDICTDIR/"
            rm -rf "\$RUNDIR"
        fi
    }

    if [ "${params.debug.toBoolean()}" = "true" ]; then
        echo "[DEBUG] out=${out} asmid=${asmid} species=${species} strain=${strain}"
        echo "[DEBUG] locustag=${locustag} busco=${busco_lineage} transl_table=${transl_table}"
        echo "[DEBUG] proteins=${params.proteins} genome_fa=${genome_fa}"
        echo "[DEBUG] PREDICTDIR=\$PREDICTDIR TMPDIR=\$TMPDIR pwd=\$(pwd)"
    fi

    # ── Skip vs. refresh decision ─────────────────────────────────────────────
    if [ -s "\$PREDICT_GBK" ]; then
        SPECIES_TAG=\$(printf '%s' "${species}" | sed -E 's/[[:space:]]+/_/g')
        STALE=0
        # Same evidence the channel-level checks use (FunannotateUtils.needsPredict
        # / staleSharedParams): the species' reads and Trinity, this genome's
        # PASA training output and, for reuse siblings only, the shared
        # ab-initio store (a representative's own prediction builds the store,
        # so for it the store is always newer). Without the last two, a genome
        # the pipeline had correctly flagged as stale (re-trained, or a new
        # representative's store) exited here as "current" and kept its old
        # annotation. Seen 2026-09-29 on the Bd pangenome's 10 pilot strains.
        for f in "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_R1.fastq.gz" \\
                 "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_SE.fastq.gz" \\
                 "${launchDir}/rnaseq_data/\${SPECIES_TAG}.trinity-GG.fasta" \\
                 "${params.training_target}/${out}/training/funannotate_train.pasa.gff3" \\
                 "${params.training_target}/${out}/training/funannotate_train.transcripts.gff3" \\
                 ${shared_params_json ? "\"${shared_params_json}\"" : ''}; do
            if [ -s "\$f" ] && [ "\$f" -nt "\$PREDICT_GBK" ]; then STALE=1; echo "[INFO] \$f is newer than the existing GBK"; fi
        done
        if [ "\$STALE" -eq 0 ]; then
            echo "[INFO] Prediction already complete and current for ${out}; nothing to do"
            touch ${out}.predict.done
            exit 0
        fi
        echo "[INFO] Stale prediction for ${out}: evidence newer than GBK — clearing predict outputs for a fresh run"
        rm -rf "\$PREDICTDIR/predict_results" "\$PREDICTDIR/predict_misc"
    fi

    mkdir -p "\$PREDICTDIR"

    # ── Guard against a corrupt partial from a previous attempt ───────────────
    if [ ! -d "\$PREDICTDIR/predict_misc" ] && [ -d "\$PREDICTDIR/predict_results" ]; then
        echo "[WARN] predict_results/ present without predict_misc/ for ${out}; clearing stale partial"
        rm -rf "\$PREDICTDIR/predict_results"
    fi

    # Local-scratch mode: always start from a clean local dir. A partial
    # predict_misc left in PREDICTDIR by an earlier non-local attempt is not
    # copied in: it is thousands of small files, the exact cost this avoids.
    if [ "\$RUNDIR" != "\$PREDICTDIR" ]; then
        rm -rf "\$RUNDIR" && mkdir -p "\$RUNDIR"
        echo "[INFO] funannotate predict working in node-local \$RUNDIR; results copied to \$PREDICTDIR at the end"
    fi

    # Point funannotate at the persistent training dir via symlink.
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "\$RUNDIR/training"
    fi

    TBL2ASN_PARAMS="-l paired-ends"

    # Inflate a gzipped clean/masked genome to a local uncompressed copy.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # funannotate predict rejects FASTA deflines longer than 24 chars; NCBI-style
    # headers survive the AAFTF clean verbatim (scripts/clean_genome_fa.py keeps
    # headers untouched), so rewrite each header to its accession (first
    # whitespace token). Idempotent -- safe on already-short headers too.
    # params.predict_defline_first_word=false passes deflines through unchanged,
    # as the BFD pipeline does.
    if [ "${params.predict_defline_first_word}" = "true" ]; then
        awk '/^>/{print \$1; next} {print}' "\$GENOME_IN" > "\$GENOME_IN.hdr" && mv "\$GENOME_IN.hdr" "\$GENOME_IN"
    fi

    # ── Too-small-genome pre-flight guard ────────────────────────────────────
    # Shared with GENEMARK_RUN, which needs the identical policy upstream of
    # this process (see genemark_run.nf, bin/asm_preflight_stats.py). BFD rule:
    # any verdict other than "ok" skips, except a small_fragmented genome that
    # has PRODIGAL_RUN evidence. The same pass reports the soft-masked share
    # (repeat_pct) for the repeat-aware EVM rule below.
    SKIP_REPORT="${params.target}/predict_skipped_too_small.tsv"
    read ASM_BP ASM_CTG ASM_N50 ASM_VERDICT ASM_REPEAT_PCT < <(
        python "${workflow.projectDir}/bin/asm_preflight_stats.py" "\$GENOME_IN" \\
            --min-bp ${params.predict_min_asm_bp} --max-n50 ${params.predict_frag_max_n50} \\
            --max-contigs ${params.predict_frag_max_contigs} \\
            --min-contig-len ${params.predict_min_training_contig_len} \\
            --min-training-contigs ${params.predict_min_training_contigs} \\
            --abs-min-bp ${params.predict_abs_min_asm_bp} \\
            --report-repeat-pct)
    echo "[INFO] Pre-flight assembly stats for ${out}: \${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}, \${ASM_REPEAT_PCT}% repeat-masked"
    if [ "\$ASM_VERDICT" != "ok" ] && ! { [ "\$ASM_VERDICT" = "small_fragmented" ] && [ "${other_gff_ok}" = "true" ]; }; then
        echo "[WARN] ${out} failed preflight ('\$ASM_VERDICT': \${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}); skipping predict" >&2
        mkdir -p "${params.target}"
        [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "preflight_\$ASM_VERDICT" "\$ASM_BP" "\$ASM_CTG" "\$ASM_N50" >> "\$SKIP_REPORT"
        touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
        touch ${out}.predict.done
        exit 0
    elif [ "\$ASM_VERDICT" = "small_fragmented" ]; then
        echo "[INFO] ${out} is small/fragmented but has Prodigal evidence (${other_gff}); proceeding instead of skipping" >&2
    fi

    # ── Repeat-aware EVM mode (BFD rule) ──────────────────────────────────────
    # See conf/profile_annotate.config predict_evm_repeat_pct_threshold.
    WEIGHT_ARGS=(${weight_args})
    # EVM weight refit (BFD pasa_train_performance_evaluate/DECISIONS.md D126/D127):
    # predict_evm_weights (genome has a PASA training set) or
    # predict_evm_weights_norna (no PASA set) is appended after genemark:1 so its
    # genemark value wins (funannotate applies -w entries in order, last wins), and
    # before the repeat-aware snap:0 below so that override still wins. Not applied
    # when Prodigal evidence (other_gff) is used. Empty = funannotate's own weights.
    if [ "${other_gff_ok}" != "true" ]; then
        if [ -s "${params.training_target}/${out}/training/funannotate_train.pasa.gff3" ]; then
            EVM_TUNED_WEIGHTS="${params.predict_evm_weights ?: ''}"
        else
            EVM_TUNED_WEIGHTS="${params.predict_evm_weights_norna ?: ''}"
        fi
        if [ -n "\$EVM_TUNED_WEIGHTS" ]; then
            echo "[INFO] ${out}: EVM weights \$EVM_TUNED_WEIGHTS"
            WEIGHT_ARGS+=(\$EVM_TUNED_WEIGHTS)
        fi
    fi
    EVM_REPEAT_FLAGS=()
    THRESHOLD=${params.predict_evm_repeat_pct_threshold}
    if [ "\$THRESHOLD" != "0" ] && awk -v p="\$ASM_REPEAT_PCT" -v t="\$THRESHOLD" 'BEGIN{exit !(p>=t)}'; then
        echo "[INFO] ${out}: \${ASM_REPEAT_PCT}% repeat-masked >= \${THRESHOLD}% -- repeat-aware EVM mode (--repeats2evm --evm-partition-interval ${params.predict_evm_repeat_aware_interval})"
        EVM_REPEAT_FLAGS=(--repeats2evm --evm-partition-interval ${params.predict_evm_repeat_aware_interval})
        if [ "${params.predict_evm_repeat_aware_drop_snap}" = "true" ]; then
            echo "[INFO] ${out}: repeat-aware mode also sets -w snap:0"
            WEIGHT_ARGS+=(snap:0)
        fi
    fi

    # ── PASA training-set gate and short-transcript guard ─────────────────────
    # Only a gate-aware funannotate (>= 1.9.0-rc.2) knows --min_pasa_complete_models;
    # rc.1 would reject it. `funannotate predict --help` prints funannotate's own usage
    # text, which lists the flag only from after rc.3, so check the installed module.
    GATE_ARGS=()
    if python3 -c "import inspect, sys, funannotate.predict as p; sys.exit(0 if 'min_pasa_complete_models' in inspect.getsource(p) else 1)" 2>/dev/null; then
        if [ -n "${params.predict_min_pasa_complete_models != null ? params.predict_min_pasa_complete_models : ''}" ]; then
            GATE_ARGS+=(--min_pasa_complete_models ${params.predict_min_pasa_complete_models})
        fi
        if [ -s "${params.training_target}/${out}/training/.trinity_short_transcripts" ]; then
            echo "[INFO] ${out}: Trinity assembly has short transcripts; training Augustus/SNAP from BUSCO"
            GATE_ARGS+=(--min_pasa_complete_models 1000000000)
        fi
    else
        echo "[INFO] ${out}: this funannotate has no PASA training-set gate (< 1.9.0-rc.2); gate args not passed"
    fi

    funannotate predict --name ${locustag} -i "\$GENOME_IN" --strain "${strain}" \\
        -o "\$RUNDIR" -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w "\${WEIGHT_ARGS[@]}" \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence ${params.proteins} \\
        --max_intronlen ${params.max_intronlen} --min_intronlen ${params.min_intronlen} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table} ${genemark_cli} ${other_gff_cli} \\
        "\${EVM_REPEAT_FLAGS[@]}" "\${GATE_ARGS[@]}" || true

    # ── Post-predict catch ────────────────────────────────────────────────────
    if [ ! -s "\$RUN_GBK" ]; then
        copy_logs_back
        # copy_logs_back has already moved the logs to PREDICTDIR and removed a local-scratch
        # RUNDIR, so read the copy (RUNDIR == PREDICTDIR when not using local scratch).
        PLOG="\$PREDICTDIR/logfiles/funannotate-predict.log"
        if [ -f "\$PLOG" ] && grep -q "Not enough gene models .* to train Augustus" "\$PLOG"; then
            NMODELS=\$(grep -oE "Not enough gene models [0-9]+" "\$PLOG" | grep -oE "[0-9]+" | tail -1)
            echo "[WARN] ${out}: funannotate found only \${NMODELS:-<min} training models (needs 30); too small/fragmented to annotate — skipping" >&2
            mkdir -p "${params.target}"
            [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "funannotate_too_few_models:\${NMODELS:-NA}" "" "" "" >> "\$SKIP_REPORT"
            touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
            touch ${out}.predict.done
            exit 0
        fi
        echo "ERROR: funannotate predict did not produce expected GBK: \$RUN_GBK" >&2
        exit 1
    fi
    if [ -d "\$RUNDIR/predict_misc/ab_initio_parameters" ]; then
        # Besides the ab-initio parameters and tRNAs, keep the small files needed
        # to diagnose a gene-count difference after the fact (same keep list as
        # BFD's FUNANNOTATE_PREDICT):
        #   weights.evm.txt            -- EVM weights actually used
        #   final_training_models.gff3 -- models Augustus/SNAP were trained on
        # Missing ones are skipped.
        KEEP_DIR="\$RUNDIR/.predict_misc_keep"
        rm -rf "\$KEEP_DIR"; mkdir -p "\$KEEP_DIR"
        for f in ab_initio_parameters trnascan.no-overlaps.gff3 weights.evm.txt final_training_models.gff3; do
            if [ -e "\$RUNDIR/predict_misc/\$f" ]; then mv "\$RUNDIR/predict_misc/\$f" "\$KEEP_DIR/"; fi
        done
        if [ -f "\$KEEP_DIR/final_training_models.gff3" ]; then pigz "\$KEEP_DIR/final_training_models.gff3"; fi
        rm -rf "\$RUNDIR/predict_misc"
        mv "\$KEEP_DIR" "\$RUNDIR/predict_misc"
    fi
    find "\$RUNDIR/predict_results/" -maxdepth 1 \\( -name "*.txt" -o -name "*.mrna-transcripts.fa" \\) -print0 \
        | xargs -0 --no-run-if-empty pigz
    sync_back
    [ -s "\$PREDICT_GBK" ] || { echo "ERROR: copy of predict results to \$PREDICTDIR failed" >&2; exit 1; }
    # New gene models invalidate BUSCO_COMPLETENESS's result. It is storeDir-cached
    # under \$PREDICTDIR/busco_completeness and storeDir skips on the directory
    # existing, regardless of which proteins produced it, so a re-prediction (e.g.
    # after RNA-seq/PASA evidence arrives) otherwise kept the old score.
    rm -rf "\$PREDICTDIR/busco_completeness"
    # Training intermediates are no longer needed once predict has succeeded
    # (bin/train_cleanup.sh keeps everything predict/update read).
    if [ "${params.train_cleanup}" = "true" ]; then
        bash "${workflow.projectDir}/bin/train_cleanup.sh" "${params.training_target}/${out}/training" ${params.run_update ? 1 : 0}
    fi
    # fsync only the delivered GBK. A bare `sync` flushes every dirty page on the node (other
    # tenants' writes included) and blocked in uninterruptible wait for > 3 h on a loaded CephFS
    # node after the task had finished, holding its pod slot (bfd_wave1 rc.4 pilot, 2026-10-02).
    sync "\$PREDICT_GBK" || true
    touch ${out}.predict.done
    echo "[INFO] Prediction complete for ${out} at \$PREDICTDIR"
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] Would run funannotate predict for ${out} using ${genome_fa}"
    [ -f "${genome_fa}" ] || [ -f "${genome_fa}.gz" ] || { echo "ERROR: genome not found at ${genome_fa}[.gz]" >&2; exit 1; }
    mkdir -p ${params.target}/${out}/predict_results ${params.target}/${out}/predict_misc
    echo "LOCUS stub_${out}" > ${params.target}/${out}/predict_results/${out}.gbk
    echo ">stub_${out}_p1" > ${params.target}/${out}/predict_results/${out}.proteins.fa
    touch ${out}.predict.done
    """
}
