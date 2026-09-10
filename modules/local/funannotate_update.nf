// Run funannotate update (PASA re-alignment with RNA-seq reads post-predict).
// Resources overridden by withName: '.*:FUNANNOTATE_UPDATE' in conf/profile_annotate.config.
process FUNANNOTATE_UPDATE {
    label 'funannotate'
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    val meta

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    def header_length = params.header_length
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no reads (empty marker file from SRA_FETCH) ──────────────────
    if [ ! -s "${r1}" ]; then
        echo "[INFO] No RNAseq reads for ${out}, skipping funannotate update"
        exit 0
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
        # than under training_target on shared storage -- see
        # funannotate_train.nf (confirmed 2026-09-01, ported from
        # BFD/Fungi_BFD_runs) for why: this datadir is pure per-task sidecar
        # infra that nothing downstream ever reads back, so a persistent,
        # species-keyed copy on shared storage just let a crashed prior
        # attempt strand an orphaned InnoDB tablespace file that broke
        # PASA's "-r" drop-then-recreate ("Directory not empty" / "database
        # exists"). Each SLURM job gets its own fresh node-local scratch
        # dir (see the FUNANNOTATE_UPDATE clusterOptions in
        # provision_ucr_hpcc.config), so that bug class is now impossible.
        MYSQL_SCRATCH=\$TMPDIR/mysql_db_${out}
        rm -rf \$MYSQL_SCRATCH
        mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
        # System-tables init happens in the branches below (each has a
        # different tool available) -- see funannotate_train.nf's identical
        # setup for the full rationale (no external ${params.mysql_datadir}
        # template dependency any more).
        cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
            { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        # 127.0.0.1, not \$MYHOSTNAME -- see funannotate_train.nf (confirmed
        # 2026-09-09): routing through the node's real hostname round-trips
        # over the cluster network fabric and fails MariaDB's host-ACL check.
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=127.0.0.1:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        PASA_MYSQL_USER=\$(grep '^MYSQL_RW_USER=' \$PASACONF | cut -d= -f2)
        PASA_MYSQL_PASS=\$(grep '^MYSQL_RW_PASSWORD=' \$PASACONF | cut -d= -f2)
        # \$MYSQL_SCRATCH (not the old, never-created "\$MYSQL_SCRATCH/mysql_db"
        # dead-bind path) is now a subdirectory of \$TMPDIR, which is already
        # bound; this SINGULARITY_BINDPATH is mostly redundant with the
        # explicit -B flags on the sidecar-container branch below but kept
        # for parity.
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH
        # ── Prefer an in-image mariadbd/mysqld_safe when present ──────────────
        # See funannotate_train.nf's identical branch for the full rationale.
        # Require BOTH a server daemon AND an install-db tool -- see
        # funannotate_train.nf's identical branch for why (bioconda's
        # mysql-libs/mysql-common ship a mysqld_safe wrapper script with no
        # real server or install-db tool behind it).
        if { command -v mariadbd >/dev/null 2>&1 || command -v mysqld_safe >/dev/null 2>&1; } && \\
           { command -v mariadb-install-db >/dev/null 2>&1 || command -v mysql_install_db >/dev/null 2>&1; }; then
            MYSQLD_BIN=\$(command -v mariadbd || command -v mysqld_safe)
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
            module load apptainer
            singularity exec -B \$MYSQL_SCRATCH/db/:/var/lib/mysql \\
                ${params.container_mariadb} sh -c \\
                'command -v mariadb-install-db || command -v mysql_install_db' \\
                > /tmp/mysql_install_bin_\$\$.txt 2>/dev/null
            MYSQL_INSTALL_BIN=\$(cat /tmp/mysql_install_bin_\$\$.txt 2>/dev/null)
            rm -f /tmp/mysql_install_bin_\$\$.txt
            if [ -n "\$MYSQL_INSTALL_BIN" ]; then
                echo "[INFO] Initializing fresh MariaDB system tables via sidecar image's \$MYSQL_INSTALL_BIN"
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
        sleep 5
        # See funannotate_train.nf for why this account has to be created
        # explicitly (mariadb-install-db only makes root accounts).
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
        \$MYSQL_CLIENT_BIN -uroot -h127.0.0.1 -P\${PORT} -e \
            "CREATE USER IF NOT EXISTS '\${PASA_MYSQL_USER}'@'127.0.0.1' IDENTIFIED BY '\${PASA_MYSQL_PASS}'; GRANT ALL ON *.* TO '\${PASA_MYSQL_USER}'@'127.0.0.1'; FLUSH PRIVILEGES;" || \
            { echo "ERROR: failed to create \${PASA_MYSQL_USER} mysql user" >&2; exit 1; }
    fi

    # Link training data into work dir so funannotate update finds it at the relative path it expects.
    mkdir -p ${out}
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "${out}/training"
    fi

    echo "[INFO] Running funannotate update for ${out}"
    funannotate update -i ${params.target}/${out} \\
        --left ${r1} --right ${r2} \\
        --cpus ${task.cpus} \\
        \$pasa_db_arg
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    EXPECTED="${params.target}/${out}/update_results/${out}.gbk"
    if [ ! -f "\$EXPECTED" ]; then
        echo "ERROR: funannotate update did not produce expected GBK: \$EXPECTED" >&2
        exit 1
    fi
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] FUNANNOTATE_UPDATE stub for ${out} (r1=${r1}, r2=${r2})"
    mkdir -p ${params.target}/${out}/update_results
    touch ${params.target}/${out}/update_results/${out}.tbl
    touch ${params.target}/${out}/update_results/${out}.gbk
    touch ${params.target}/${out}/update_results/${out}.gff3
    """
}
