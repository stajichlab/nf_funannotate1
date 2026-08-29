// Download antiSMASH's reference databases (ClusterBlast/MIBiG/PFAM/etc) into
// a local directory (params.antismash_db). The antismash-lite biocontainer
// ships the antismash Python package only -- no reference data baked in
// (confirmed: ANTISMASH_RUN fails with "No matching database in location
// .../antismash/databases/knownclusterblast" without this) -- so, like
// funannotate's own setup step, this has to run once before ANTISMASH_RUN can
// work. storeDir caches the populated directory at params.antismash_db, so
// this runs at most once across all pipeline runs.
process SETUP_ANTISMASH_DB {
    label 'antismash'
    label 'process_low'

    storeDir { file(params.antismash_db).parent }

    output:
    path "${db_dir}", emit: db

    script:
    db_dir = file(params.antismash_db).name
    """
    set -euo pipefail
    download-antismash-databases --database-dir ${db_dir}
    echo "[INFO] antismash database built at ${db_dir}"
    """

    stub:
    db_dir = file(params.antismash_db).name
    """
    mkdir -p ${db_dir}
    """
}
