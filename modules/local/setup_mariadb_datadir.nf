// Build a fresh MariaDB data directory once (mariadb-install-db), storeDir-cached
// at params.mysql_datadir. funannotate_train.nf / funannotate_update.nf cp -a
// just the "mysql" system-schema subfolder out of this into each task's private
// scratch datadir to seed a disposable per-task mysqld instance for PASA
// (--pasa_db mysql; see params.pasa_mysql) -- so this only needs to produce that
// one seed, not a live/growing database.
//
// Runs mariadb-install-db INSIDE the MariaDB image (params.container_mariadb,
// else params.container_funannotate; see FunannotateUtils.mariadbImage) (via `apptainer exec`
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
    // PASA MariaDB image (container_mariadb, else container_funannotate), as the
    // local file in params.sif_dir -- a docker:// URI maps to Nextflow's own
    // cache file, pulled once if missing, never converted on every call.
    def mariadb_image      = FunannotateUtils.mariadbImage(params)
    def mariadb_img        = FunannotateUtils.localImageFile(mariadb_image, params.sif_dir as String)
    def ensure_mariadb_img = FunannotateUtils.ensureLocalImageScript(mariadb_image, params.sif_dir as String)
    """
    set -euo pipefail
    APPTAINER_BIN=\$(command -v apptainer || command -v singularity || true)
    if [ -z "\$APPTAINER_BIN" ]; then
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
        APPTAINER_BIN=\$(command -v apptainer || command -v singularity || true)
    fi
    [ -n "\$APPTAINER_BIN" ] || {
        echo "ERROR: no apptainer/singularity binary on PATH (needed to run mariadb-install-db from ${mariadb_img})" >&2
        exit 1
    }
    ${ensure_mariadb_img}
    mkdir -p ${datadir_name} mariadb_tmp
    # TMPDIR -> a dir under the (bound) task workdir: the host \$TMPDIR (e.g.
    # SLURM's /scratch/\$USER/<jobid>) is passed into the container but not
    # bound, and MariaDB 11.8 (funannotate image) writes InnoDB temp files
    # there -> "Read-only file system", install aborts (confirmed 2026-09-24;
    # the old 10.3.9 mariadb.sif did not hit this).
    "\$APPTAINER_BIN" exec -B "\$PWD":"\$PWD" '${mariadb_img}' \\
        bash -c 'export TMPDIR="\$2"; DB="\$(command -v mariadb-install-db || command -v mysql_install_db)"; [ -n "\$DB" ] || DB=mysql_install_db; echo "[INFO] using \$DB"; exec "\$DB" --datadir="\$1" --auth-root-authentication-method=normal' _ \\
        "\$PWD/${datadir_name}" "\$PWD/mariadb_tmp"
    echo "[INFO] MariaDB seed datadir built at ${datadir_name}"
    """

    stub:
    datadir_name = file(params.mysql_datadir).name
    """
    mkdir -p ${datadir_name}/mysql
    """
}
