// Build the funannotate databases into a local directory (params.funannotate_db).
// Two-pass: BUSCO lineage DBs first (-b all -i busco), then all remaining databases
// (-i all). Runs under the 'funannotate' label so the DBs are built with whichever
// funannotate the active provisioning profile supplies (module / pixi / singularity).
// storeDir caches the populated directory at params.funannotate_db, so this runs at
// most once across all pipeline runs; if the directory already exists (e.g. pointed
// at a prebuilt shared DB) the task is skipped entirely but still emits `db`.
process SETUP_FUNANNOTATE_DB {
    label 'funannotate'
    label 'process_low'

    // Closure defers evaluation to task runtime so a missing pipeline profile is
    // caught by the workflow guard (clear message) instead of throwing here.
    storeDir { file(params.funannotate_db).parent }

    output:
    path "${db_dir}", emit: db

    script:
    db_dir = file(params.funannotate_db).name
    """
    set -euo pipefail
    export FUNANNOTATE_DB=\$(readlink -f ${db_dir})
    funannotate setup -d ${db_dir} -b all -i busco
    funannotate setup -d ${db_dir} -i all
    echo "[INFO] funannotate database built at ${db_dir}"
    """

    stub:
    db_dir = file(params.funannotate_db).name
    """
    mkdir -p ${db_dir}
    : > ${db_dir}/funannotate-db-info.txt
    """
}
