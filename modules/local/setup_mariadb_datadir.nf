// Build a fresh MariaDB data directory once (mariadb-install-db), storeDir-cached
// at params.mysql_datadir. funannotate_train.nf / funannotate_update.nf cp -a
// just the "mysql" system-schema subfolder out of this into each task's private
// scratch datadir to seed a disposable per-task mysqld instance for PASA
// (--pasa_db mysql; see params.pasa_mysql) -- so this only needs to produce that
// one seed, not a live/growing database.
//
// Runs mariadb-install-db INSIDE params.container_mariadb (via `apptainer exec`
// / `singularity exec`) rather than requiring a host mariadb-install-db, so no
// site-specific pre-built datadir is needed -- this works under any provisioning
// profile as long as an apptainer/singularity binary is on PATH (falling back to
// the UCR HPCC Lmod module below when it isn't).
process SETUP_MARIADB_DATADIR {
    label 'setup'
    label 'process_single'

    storeDir { file(params.mysql_datadir).parent }

    output:
    path "${datadir_name}", emit: ready

    script:
    datadir_name = file(params.mysql_datadir).name
    """
    set -euo pipefail
    APPTAINER_BIN=\$(command -v apptainer || command -v singularity || true)
    if [ -z "\$APPTAINER_BIN" ]; then
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
        APPTAINER_BIN=\$(command -v apptainer || command -v singularity || true)
    fi
    [ -n "\$APPTAINER_BIN" ] || {
        echo "ERROR: no apptainer/singularity binary on PATH (needed to run mariadb-install-db from ${params.container_mariadb})" >&2
        exit 1
    }
    mkdir -p ${datadir_name}
    "\$APPTAINER_BIN" exec -B "\$PWD":"\$PWD" ${params.container_mariadb} \\
        bash -c 'DB="\$(command -v mariadb-install-db || command -v mysql_install_db)"; [ -n "\$DB" ] || DB=mysql_install_db; echo "[INFO] using \$DB"; exec "\$DB" --datadir="\$1" --auth-root-authentication-method=normal' _ \\
        "\$PWD/${datadir_name}"
    echo "[INFO] MariaDB seed datadir built at ${datadir_name}"
    """

    stub:
    datadir_name = file(params.mysql_datadir).name
    """
    mkdir -p ${datadir_name}/mysql
    """
}
