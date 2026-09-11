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
    # funannotate 1.9's trinity.py falls back to \$TMPDIR for Trinity's
    # --workdir when funannotate train isn't given one explicitly (which we
    # never do here). Trinity itself refuses to run unless "trinity" appears
    # literally in that path (its own safety check against auto-deleting the
    # wrong directory on cleanup) -- confirmed 2026-09-08 failing every
    # genome under v1.9.0-beta.11 (both perl and Rust Trinity, since they
    # share this same Python wrapper) with a bare \$SCRATCH path
    # (/scratch/<user>/<jobid>, no "trinity" in it). 1.8.17's older
    # trinity.py has no such fallback, so this never affected v1.8.17_conda.
    # Must be `export`ed: TMPDIR here was previously a plain (non-exported)
    # local var, so `funannotate train` was actually inheriting SLURM's own
    # auto-exported TMPDIR (== raw \$SCRATCH) instead of this computed value.
    TMPDIR="\$TMPDIR/trinity_work"
    mkdir -p "\$TMPDIR"
    export TMPDIR
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
        # System-tables init moved into the branches below (each has a
        # different tool available to run it) -- see there instead of a
        # ${params.mysql_datadir}/mysql template copy. A pre-built,
        # externally-staged datadir template was previously required before
        # this pipeline could even start (not self-contained: the file has to
        # already exist at a specific host path outside version control), and
        # every task paid its full copy cost (previously ~121MB/89 files via
        # cp -a) even though `mariadb-install-db`/`mysql_install_db` initializes
        # fresh system tables in a few seconds regardless.
        cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
            { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        # 127.0.0.1, not \$MYHOSTNAME: mariadbd binds loopback-only
        # (my.cnf's bind-address) and PASA runs in the same node/namespace,
        # so this connection never needs to leave loopback. Routing it
        # through the node's real hostname instead round-trips over the
        # cluster network fabric and made MariaDB's host-ACL check reject
        # the connection (confirmed 2026-09-09 against Fungi_BFD's identical
        # setup -- reverse-DNS resolved the peer to the InfiniBand FQDN,
        # which matched none of the ACL entries mariadb-install-db creates).
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=127.0.0.1:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        # Read the account PASA will actually connect as straight out of
        # \$PASACONF rather than hardcoding it here, so assets/pasa_conf/
        # conf.txt stays the single source of truth for these credentials.
        PASA_MYSQL_USER=\$(grep '^MYSQL_RW_USER=' \$PASACONF | cut -d= -f2)
        PASA_MYSQL_PASS=\$(grep '^MYSQL_RW_PASSWORD=' \$PASACONF | cut -d= -f2)
        # Point MariaDB's on-disk temp-table dir at this job's own node-local
        # \$TMPDIR (== \$SCRATCH) instead of assets/pasa_conf/my.cnf's default
        # of /tmp -- mariadbd shares the host mount namespace (no
        # --containall), so its "/tmp" is literally the compute node's real,
        # SHARED /tmp; when several FUNANNOTATE_TRAIN tasks land on the same
        # node at once (routine here -- this benchmark runs every genome in
        # a cell in parallel), their MariaDB instances all write large
        # PASA-alignment MyISAM temp tables into that same shared /tmp
        # concurrently, and once it fills, mysqld fails mid-write on a temp
        # table and then fails again trying to delete the file it never
        # finished creating (confirmed against Fungi_BFD 2026-08-29,
        # ported here 2026-09-10 -- ${out} would otherwise be exposed to
        # the identical failure mode once enough cells run concurrently).
        MYSQL_TMP="\$TMPDIR/pasa_mysql_tmp_${asmid}"
        mkdir -p "\$MYSQL_TMP"
        sed -i "s#^tmpdir[[:space:]]*=.*#tmpdir\t\t= \$MYSQL_TMP#" \$MYSQL_SCRATCH/conf/my.cnf
        # \$MYSQL_SCRATCH (not the old, never-created "\$MYSQL_SCRATCH/mysql_db"
        # dead-bind path) is now a subdirectory of \$TMPDIR, which is already
        # bound; this SINGULARITY_BINDPATH is mostly redundant with the
        # explicit -B flags on the sidecar-container branch below but kept
        # for parity.
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH
        # ── Prefer an in-image mariadbd/mysqld_safe when present ──────────────
        # Under the container/singularity provisioning profile, this whole
        # task already runs INSIDE params.container_funannotate via
        # Nextflow's own container wrap. The old unconditional path below
        # then tried to nest a SECOND singularity container (the MariaDB
        # sidecar) from inside that wrap via `module load apptainer` +
        # `singularity instance start` -- neither Lmod nor an
        # apptainer/singularity client exist inside the funannotate image,
        # so that always failed under container mode (confirmed 2026-09-08,
        # every genome, both v1.8.17_container and v1.9.0-beta11_container_rust).
        # Once funannotate-live's Dockerfile.base bundles mariadb-server
        # (see its apt-get install block), mysqld_safe/mariadbd is on PATH
        # inside the SAME container the task is already running in, so we
        # can just start it in-place -- no nesting, no module/singularity
        # dependency at all. This is a runtime capability check, not a
        # provisioning-profile branch: conda-profile tasks run bare on the
        # host (no funannotate container involved) and essentially never
        # have a local mariadbd on PATH, so they transparently keep using
        # the sidecar-container fallback below, unchanged from today's
        # working behavior.
        # Require BOTH a server daemon AND an install-db tool -- bioconda's
        # mysql-libs/mysql-common ship a mysqld_safe *wrapper script* (with
        # no actual mariadbd/mysqld behind it) but no install-db tool at all,
        # which false-positived this check into the in-image branch with
        # nothing to actually run (confirmed 2026-09-09/10, every conda
        # cell: "bundled mariadbd/mysqld_safe found but no
        # mariadb-install-db/mysql_install_db on PATH").
        if { command -v mariadbd >/dev/null 2>&1 || command -v mysqld_safe >/dev/null 2>&1; } && \\
           { command -v mariadb-install-db >/dev/null 2>&1 || command -v mysql_install_db >/dev/null 2>&1; }; then
            MYSQLD_BIN=\$(command -v mariadbd || command -v mysqld_safe)
            # Fresh in-image init -- no external datadir template needed.
            # NOTE: exact flag name/availability not yet verified against the
            # rebuilt image's actual mariadb-server package (Debian trixie);
            # confirm `mariadb-install-db --help` there and adjust if needed.
            MYSQL_INSTALL_BIN=\$(command -v mariadb-install-db || command -v mysql_install_db)
            if [ -n "\$MYSQL_INSTALL_BIN" ]; then
                echo "[INFO] Initializing fresh MariaDB system tables via \$MYSQL_INSTALL_BIN"
                "\$MYSQL_INSTALL_BIN" --datadir=\$MYSQL_SCRATCH/db/mysql \\
                    --auth-root-authentication-method=normal || \\
                    { echo "ERROR: \$MYSQL_INSTALL_BIN failed" >&2; exit 1; }
            else
                echo "ERROR: bundled mariadbd/mysqld_safe found but no mariadb-install-db/mysql_install_db on PATH" >&2
                exit 1
            fi
            echo "[INFO] Using in-image \$MYSQLD_BIN for the PASA MariaDB backend (no sidecar container needed)"
            "\$MYSQLD_BIN" --defaults-file=\$MYSQL_SCRATCH/conf/my.cnf \\
                --datadir=\$MYSQL_SCRATCH/db/mysql \\
                --socket=\$MYSQL_SCRATCH/mysqld.sock \\
                --pid-file=\$MYSQL_SCRATCH/mysqld.pid &
            MYSQLD_PID=\$!
            stop_mysqldb() { kill \$MYSQLD_PID 2>/dev/null || true; wait \$MYSQLD_PID 2>/dev/null || true; }
        else
            stop_mysqldb() { singularity instance stop mysqldb_${asmid}_\${SLURM_JOB_ID:-\$\$} 2>/dev/null || true; }
            # apptainer (not the old `singularity` module) so squashfuse is
            # pulled in automatically -- see conf/provision_singularity.config
            # for the same rationale on the main container axis. The
            # `singularity` binary used below is apptainer's own compat symlink.
            module load apptainer
            # Fresh init via the sidecar image's OWN bundled install-db tool
            # (same rationale as the in-image branch above: no external
            # datadir template needed) -- NOT yet verified that
            # params.container_mariadb actually has mariadb-install-db/
            # mysql_install_db on its PATH; confirm before relying on this.
            singularity exec -B \$MYSQL_SCRATCH/db/:/var/lib/mysql \\
                ${params.container_mariadb} sh -c \\
                'command -v mariadb-install-db || command -v mysql_install_db' \\
                > /tmp/mysql_install_bin_\$\$.txt 2>/dev/null
            MYSQL_INSTALL_BIN=\$(cat /tmp/mysql_install_bin_\$\$.txt 2>/dev/null)
            rm -f /tmp/mysql_install_bin_\$\$.txt
            if [ -n "\$MYSQL_INSTALL_BIN" ]; then
                echo "[INFO] Initializing fresh MariaDB system tables via sidecar image's \$MYSQL_INSTALL_BIN"
                # --datadir=/var/lib/mysql, NOT /var/lib/mysql/mysql: the
                # `instance start` below (and assets/pasa_conf/my.cnf's own
                # `datadir = /var/lib/mysql`) reads tables straight out of
                # this same bind point with no extra nesting -- an earlier
                # version of this line added a spurious /mysql suffix, so
                # mysqld_safe found nothing there and never actually started
                # listening (confirmed 2026-09-10: every genome failed with
                # "Can't connect to MySQL server ... (111)" despite
                # `instance start` itself reporting success).
                singularity exec -B \$MYSQL_SCRATCH/db/:/var/lib/mysql \\
                    ${params.container_mariadb} \\
                    "\$MYSQL_INSTALL_BIN" --datadir=/var/lib/mysql \\
                    --auth-root-authentication-method=normal || \\
                    { echo "ERROR: sidecar \$MYSQL_INSTALL_BIN failed" >&2; exit 1; }
            else
                echo "ERROR: no mariadb-install-db/mysql_install_db found in ${params.container_mariadb}" >&2
                exit 1
            fi
            singularity instance start --writable-tmpfs \\
                -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
                ${params.container_mariadb} mysqldb_${asmid}_\${SLURM_JOB_ID:-\$\$} /usr/bin/mysqld_safe
        fi
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        pasa_db_arg="--pasa_db mysql"
        if command -v mariadb >/dev/null 2>&1 || command -v mysql >/dev/null 2>&1; then
            MYSQL_CLIENT_BIN=\$(command -v mariadb || command -v mysql)
        else
            MYSQL_CLIENT_BIN=\$(singularity exec ${params.container_mariadb} sh -c 'command -v mariadb || command -v mysql' 2>/dev/null)
            MYSQL_CLIENT_BIN="singularity exec ${params.container_mariadb} \$MYSQL_CLIENT_BIN"
        fi
        if [ -z "\$MYSQL_CLIENT_BIN" ]; then
            echo "ERROR: no mariadb/mysql client found" >&2
            exit 1
        fi
        # Poll for the listener instead of a fixed sleep: mariadb-install-db +
        # mysqld_safe startup time varies with node load/cold page cache, and
        # a fixed short delay is exactly what produced "Can't connect to
        # MySQL server ... (111)" despite "instance started successfully"
        # (confirmed 2026-09-10, every conda cell -- the daemon just wasn't
        # listening yet). 30 x 1s covers the slowest cold-start case seen so
        # far with a wide margin; a still-failed daemon after that is treated
        # as a real infra failure, not silently waited on forever.
        MYSQL_READY=0
        for _ in \$(seq 1 30); do
            if \$MYSQL_CLIENT_BIN -uroot -h127.0.0.1 -P\${PORT} -e 'SELECT 1' >/dev/null 2>&1; then
                MYSQL_READY=1
                break
            fi
            sleep 1
        done
        if [ "\$MYSQL_READY" -ne 1 ]; then
            echo "ERROR: mariadbd on 127.0.0.1:\${PORT} did not become ready within 30s" >&2
            exit 1
        fi
        # mariadb-install-db (--auth-root-authentication-method=normal, above)
        # only creates root@localhost/127.0.0.1/::1/<hostname> with no
        # password -- it never creates the account conf.txt tells PASA to
        # connect as. Confirmed 2026-09-09 (same bug hit in Fungi_BFD): this
        # is a fresh, throwaway, loopback-only DB that lives for one task, so
        # a shared generic account (assets/pasa_conf/conf.txt) is fine.
        # Grant to '...'@'127.0.0.1', not '...@localhost': MariaDB's ACL
        # matches the literal connecting host/IP, and a TCP connection to
        # 127.0.0.1 is not treated as 'localhost' (reserved for Unix-socket
        # connections) -- same ACL-mismatch class as the MYSQLSERVER fix above.
        \$MYSQL_CLIENT_BIN -uroot -h127.0.0.1 -P\${PORT} -e \
            "CREATE USER IF NOT EXISTS '\${PASA_MYSQL_USER}'@'127.0.0.1' IDENTIFIED BY '\${PASA_MYSQL_PASS}'; GRANT ALL ON *.* TO '\${PASA_MYSQL_USER}'@'127.0.0.1'; FLUSH PRIVILEGES;" || \
            { echo "ERROR: failed to create \${PASA_MYSQL_USER} mysql user" >&2; exit 1; }
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
        # Trinity's genome-guided checkpoint files under trinity_gg/ hardcode
        # this job's node-local \$SCRATCH path (e.g. /scratch/\$USER/\$SLURM_JOB_ID).
        # A retry runs as a brand-new SLURM job with a brand-new \$SCRATCH, so
        # resuming from these checkpoints fails with "mkdir: cannot create
        # directory '/scratch/.../<dead job id>': Permission denied" and
        # silently produces "0 transcripts derived from Trinity" -- an
        # infinite failure loop across retries. hisat2/ has the same
        # cross-job staleness risk if the alignment was interrupted mid-write
        # (an empty/truncated BAM gets trusted as "existing alignments found"
        # on resume). Wipe both here so a retry always re-aligns and
        # re-assembles from scratch instead of trusting a checkpoint tied to
        # this now-dead job. Discovered 2026-09-11 in
        # BFD/Funannotate_benchmarking's conda-provisioned cells.
        TRAINDIR="${params.training_target}/${out}/training"
        rm -rf "\$TRAINDIR/hisat2"
        rm -rf "\$TRAINDIR/trinity_gg"
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
