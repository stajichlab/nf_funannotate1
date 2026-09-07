// Run funannotate train (PASA alignment) for a single assembly, using shared Trinity-GG
// transcripts produced by RNASEQ_PREPARE (representative strain) or running a full
// Trinity+PASA train when no shared Trinity is available (fallback).
// Writes output directly to params.training_target/<id>/. No publishDir — the persistent
// training directory is the primary output, accessed by FUNANNOTATE_PREDICT via symlink.
// Resources overridden by withName: '.*:FUNANNOTATE_TRAIN' in conf/profile_annotate.config.
process FUNANNOTATE_TRAIN {
    label 'funannotate'
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), val(genome_fa), path(r1), path(r2), path(se), path(trinity_fa), val(pasa_tier)

    output:
    tuple val(meta), val(genome_fa), emit: predict_input
    // Audit row emitted only when a run gracefully degrades to ab-initio-only
    // because PASA completed alignment/assignment but assigned too few loci to
    // build a training set (see the TRAIN_STATUS handling near the end).
    // Collected by TRAIN_PREDICT into one reviewable TSV. Ported from
    // BFD/Fungi_BFD_runs's modules/funannotate/predict/FUNANNOTATE_TRAIN/main.nf
    // (there generalized to every pasa_tier 2026-09-06 after a forensic sweep
    // found 207/207 real-run FUNANNOTATE_TRAIN failures shared this exact
    // signature -- PASA finished but assigned a median of 5 loci, max 249, all
    // < 500 -- across every tier, not just the hybrid-composite ones this
    // degrade path originally targeted).
    path("${meta.id}.pasa_train_failed.tsv"), optional: true, emit: pasa_failed

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def header_length = params.header_length
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no RNA-seq data at all ────────────────────────────────────────
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ] && [ ! -s "${trinity_fa}" ]; then
        echo "[INFO] No RNAseq data for ${out}, skipping funannotate train"
        exit 0
    fi

    # ── Skip if the shared Trinity-GG assembly is too thin to train on ────────
    # A too-few-transcripts trinity_fa usually means it was assembled against the
    # wrong reference strain (same species name, divergent genome) rather than a
    # real expression signal -- PASA has nothing to build a training set from and
    # funannotate train either crashes or emits junk models. See
    # train_min_trinity_transcripts in conf/profile_annotate.config.
    if [ -s "${trinity_fa}" ] && [ "${params.train_min_trinity_transcripts}" -gt 0 ]; then
        TRINITY_TX_COUNT=\$(grep -c '^>' "${trinity_fa}" || true)
        if [ "\$TRINITY_TX_COUNT" -lt "${params.train_min_trinity_transcripts}" ]; then
            echo "[WARN] ${out}: shared Trinity-GG assembly has only \$TRINITY_TX_COUNT transcripts (< ${params.train_min_trinity_transcripts}); likely assembled against the wrong reference strain. Skipping funannotate train." >&2
            exit 0
        fi
    fi

    # ── Skip if training is already resolved and evidence is not newer ────────────
    # "Resolved" is either a real published PASA GFF3, OR a durable marker recording
    # that an attempt legitimately couldn't train from transcript evidence -- PASA
    # completed alignment/assignment but assigned too few loci to build a training
    # set (see the TRAIN_STATUS handling near the end). Without this marker, a
    # gracefully-degraded strain would never publish anything and would be
    # re-attempted on every single future run/-resume forever -- the marker makes
    # that outcome durable, the same way a real GFF3 makes success durable.
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    TRAIN_FAILED_MARKER="${params.training_target}/${out}/training/.pasa_train_failed"
    PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk"
    [ -f "\$PREDICT_GBK" ] || PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk.gz"
    RESOLVED_MARKER=""
    if [ -f "\$TRAIN_GFF3" ]; then
        RESOLVED_MARKER="\$TRAIN_GFF3"
    elif [ -f "\$TRAIN_FAILED_MARKER" ]; then
        RESOLVED_MARKER="\$TRAIN_FAILED_MARKER"
    fi
    if [ -n "\$RESOLVED_MARKER" ]; then
        RETRAIN=0
        # Re-train if the trinity evidence itself is newer than the resolved marker
        # -- covers a rebuilt/rebalanced shared Trinity for a strain previously
        # marked not-trainable.
        if [ -s "${trinity_fa}" ] && [ "${trinity_fa}" -nt "\$RESOLVED_MARKER" ]; then
            echo "[INFO] Trinity evidence newer than training output for ${out}; retraining"
            RETRAIN=1
        fi
        if [ -f "\$PREDICT_GBK" ]; then
            if [ -s "${r1}" ] && [ "${r1}" -nt "\$PREDICT_GBK" ]; then
                echo "[INFO] RNAseq R1 reads newer than predict GBK for ${out}; retraining"
                RETRAIN=1
            elif [ -s "${se}" ] && [ "${se}" -nt "\$PREDICT_GBK" ]; then
                echo "[INFO] RNAseq SE reads newer than predict GBK for ${out}; retraining"
                RETRAIN=1
            fi
        fi
        if [ \$RETRAIN -eq 0 ]; then
            if [ "\$RESOLVED_MARKER" = "\$TRAIN_GFF3" ]; then
                echo "[INFO] Training already complete for ${out}; skipping"
            else
                echo "[INFO] ${out} previously determined not trainable from transcript evidence (pasa_tier=${pasa_tier}); skipping re-attempt"
            fi
            exit 0
        fi
    fi

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    # Node-local scratch may not exist / be writable if an inherited \$SCRATCH
    # points at another node's path — fall back to the task workdir.
    TMPDIR=\$(printf '%s' "\${SCRATCH:-}" | tr -d '\\n\\r')
    TMPDIR=\${TMPDIR:-/tmp}
    if [ ! -d "\$TMPDIR" ] || [ ! -w "\$TMPDIR" ]; then
        TMPDIR="\$PWD"
    fi
    export PASACONF=""
    pasa_db_arg="--pasa_db sqlite"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        # Lives under \$TMPDIR (== node-local \$SCRATCH under SLURM) rather
        # than under training_target on shared storage. This datadir is pure
        # per-task sidecar infra for PASA's MySQL backend -- nothing
        # downstream ever reads it back -- so a persistent, species-keyed
        # copy on shared storage just let a crashed prior attempt strand an
        # orphaned InnoDB tablespace file that MariaDB's own metadata didn't
        # know about; PASA's "-r" (drop-then-recreate) DROP DATABASE then
        # couldn't rmdir the stale directory (errno 39 "Directory not
        # empty"), and the following CREATE DATABASE failed with "database
        # exists". Same bug confirmed 2026-09-01 in BFD/Fungi_BFD_runs's
        # FUNANNOTATE_TRAIN (orphaned splice_variation.ibd from a crashed
        # attempt); ported here since this module has the identical pattern.
        # Putting it on \$SCRATCH instead makes that bug class impossible:
        # each SLURM job gets its own fresh, node-local scratch dir that
        # SLURM tears down when the job ends (see the FUNANNOTATE_TRAIN
        # clusterOptions in provision_ucr_hpcc.config for the
        # SCRATCH=,TMPDIR= clearing this depends on).
        MYSQL_SCRATCH=\$TMPDIR/mysql_db_${out}
        rm -rf \$MYSQL_SCRATCH
        mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
        # cp -a (not rsync): this runs inside the funannotate container image,
        # which has no rsync binary; the datadir is a one-shot bootstrap into a
        # fresh empty dir, and cp -a preserves perms/times/symlinks identically.
        cp -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
            { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
        cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
            { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        # \$MYSQL_SCRATCH (not the old, never-created "\$MYSQL_SCRATCH/mysql_db"
        # dead-bind path) is now a subdirectory of \$TMPDIR, which is already
        # bound; this SINGULARITY_BINDPATH is mostly redundant with the
        # explicit -B flags on `instance start` below but kept for parity.
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        # apptainer (not the old `singularity` module) so squashfuse is pulled
        # in automatically -- see conf/provision_singularity.config for the
        # same rationale on the main container axis. The `singularity` binary
        # used below is apptainer's own compat symlink.
        module load apptainer
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.container_mariadb} mysqldb_${asmid} /usr/bin/mysqld_safe
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # Inflate a gzipped clean genome to a local uncompressed copy.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # ── pasa_tier -> relaxed PASA alignment thresholds ────────────────────────
    # 'stringent' (PASA's own defaults, untouched) is the only value this
    # pipeline currently ever assigns -- see the module-level doc comment and
    # subworkflows/local/train_predict.nf for what's not ported yet
    # (ANI-driven representative/sibling distance, hybrid-cross composite
    # Trinity). Kept as a real branch, not dead code, so a real tier can be
    # wired in later without touching this module again.
    PASA_TIER_ARGS=""
    if [ "${pasa_tier}" = "relaxed" ]; then
        PASA_TIER_ARGS="--pasa_min_avg_per_id ${params.pasa_shared_min_avg_per_id} --pasa_min_pct_aligned ${params.pasa_shared_min_pct_aligned} --pasa_num_bp_splice ${params.pasa_shared_num_bp_splice}"
    elif [ "${pasa_tier}" = "composite" ]; then
        PASA_TIER_ARGS="--pasa_min_avg_per_id ${params.pasa_composite_min_avg_per_id} --pasa_min_pct_aligned ${params.pasa_composite_min_pct_aligned} --pasa_num_bp_splice ${params.pasa_composite_num_bp_splice}"
    elif [ "${pasa_tier}" = "composite_fallback" ]; then
        PASA_TIER_ARGS="--pasa_min_avg_per_id ${params.pasa_composite_fallback_min_avg_per_id} --pasa_min_pct_aligned ${params.pasa_composite_fallback_min_pct_aligned} --pasa_num_bp_splice ${params.pasa_composite_fallback_num_bp_splice}"
    fi

    # Whole invocation wrapped in a group + tee so a failure can be inspected
    # below without re-running anything -- \${PIPESTATUS[0]} (not plain \$?,
    # which after a pipe would report tee's exit code) captures the group's
    # real exit status, i.e. whichever funannotate train branch ran last.
    {
    # ── Use shared Trinity transcripts (PASA only) or run full train ──────────
    if [ -s "${trinity_fa}" ]; then
        if [ -s "${r1}" ]; then
            echo "[INFO] Running funannotate train (PASA+PE) for ${out} using shared Trinity (pasa_tier=${pasa_tier})"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --left_norm ${r1} --right_norm ${r2} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$PASA_TIER_ARGS \\
                \$pasa_db_arg
        elif [ -s "${se}" ]; then
            echo "[INFO] Running funannotate train (PASA+SE) for ${out} using shared Trinity (pasa_tier=${pasa_tier})"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --single_norm ${se} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$PASA_TIER_ARGS \\
                \$pasa_db_arg
        else
            # No reads at all -- r1/se are present-but-empty (0-byte) placeholders,
            # not missing paths (Nextflow path inputs can't be null). Do NOT pass
            # --left_norm/--right_norm pointed at those empty files: that hands
            # funannotate empty FASTQs to normalize instead of omitting the flags
            # entirely (same latent bug fixed in BFD/Fungi_BFD_runs's
            # FUNANNOTATE_TRAIN alongside its composite-tier work).
            echo "[INFO] Running funannotate train (PASA only, no reads) for ${out} using shared Trinity (pasa_tier=${pasa_tier})"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$PASA_TIER_ARGS \\
                \$pasa_db_arg
        fi
    elif [ -s "${r1}" ]; then
        echo "[INFO] Running funannotate train (full PE, no shared Trinity) for ${out}"
        funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    else
        echo "[INFO] Running funannotate train (full SE, no shared Trinity) for ${out}"
        funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
            --single_norm ${se} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    fi
    } 2>&1 | tee funannotate_train_capture.log
    TRAIN_STATUS=\${PIPESTATUS[0]}
    if [ "\$TRAIN_STATUS" -ne 0 ]; then
        # One chance to degrade gracefully instead of the usual hard-fail+retry,
        # PROVIDED PASA itself actually completed its alignment/assignment step
        # (its own "PASA assigned N transcripts to M loci" summary line appears
        # in the captured output). Ported from BFD/Fungi_BFD_runs's
        # FUNANNOTATE_TRAIN -- see the pasa_failed output's doc comment above
        # for the forensic basis. Absence of that line still means PASA (or
        # something before it -- MySQL, OOM, a crash) never finished, which is
        # exactly the kind of infra failure that SHOULD keep retrying with more
        # memory/resources, not get silently absorbed into "this strain isn't
        # trainable".
        if grep -qE 'PASA assigned [0-9]+ transcripts to [0-9]+ loci' funannotate_train_capture.log; then
            echo "[WARN] ${out}: funannotate train failed (exit \$TRAIN_STATUS) but PASA completed alignment/assignment for this run (pasa_tier=${pasa_tier}) -- treating as 'not enough usable transcript evidence' rather than an infra failure. Degrading to ab-initio-only; predict will proceed without PASA evidence for this strain." >&2
            mkdir -p "${params.training_target}/${out}/training"
            : > "${params.training_target}/${out}/training/.pasa_train_failed"
            printf "out\\tspecies\\tpasa_tier\\texit_code\\ttimestamp\\n%s\\t%s\\t%s\\t%s\\t%s\\n" \\
                "${out}" "${species}" "${pasa_tier}" "\$TRAIN_STATUS" "\$(date -Iseconds)" \\
                > "${out}.pasa_train_failed.tsv"
            if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
            exit 0
        fi
        echo "[ERROR] funannotate train failed for ${out} (exit \$TRAIN_STATUS)" >&2
        if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
        exit "\$TRAIN_STATUS"
    fi

    # ── Remove large intermediates not needed for predict or update ─────────────
    TRAINDIR="${params.training_target}/${out}/training"
    echo "[INFO] Removing large training intermediates in \$TRAINDIR"
    rm -rf "\$TRAINDIR/hisat2"
    rm -rf "\$TRAINDIR/trinity_gg"
    echo "[INFO] Training cleanup complete for ${out}"
    echo "mysql is ${params.pasa_mysql}"
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] FUNANNOTATE_TRAIN stub for ${out}"
    mkdir -p ${params.training_target}/${out}/training
    touch ${params.training_target}/${out}/training/funannotate_train.pasa.gff3
    """
}
