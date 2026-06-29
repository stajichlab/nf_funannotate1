/*
 * SETUP_DBS — build/seed the three run-once databases
 *
 * Calls SETUP_TAXONDB, SETUP_FUNANNOTATE_DB, and SETUP_AUGUSTUS_CONFIG.
 * All three use storeDir so they are no-ops on any run where their target
 * directories already exist.
 *
 * Emits three gate channels consumed separately by the workflow:
 *   taxondb  — val(String path)        gates GENOME_CLEAN / GENOME_CLEAN_BATCH
 *   db       — path(funannotate_db)    \
 *   config   — path(augustus_config)   / combined to gate the predict chain
 *
 * Note: storeDir processes do not emit versions.yml — adding it as an output
 * would break the storeDir cache for any pre-existing installation where only
 * the primary outputs exist. Version info for these once-run DBs is static.
 */

include { SETUP_TAXONDB        } from './../../modules/local/setup_taxondb'
include { SETUP_FUNANNOTATE_DB } from './../../modules/local/setup_funannotate_db'
include { SETUP_AUGUSTUS_CONFIG } from './../../modules/local/setup_augustus_config'

workflow SETUP_DBS {

    main:
    SETUP_TAXONDB()
    SETUP_FUNANNOTATE_DB()
    SETUP_AUGUSTUS_CONFIG()

    emit:
    taxondb = SETUP_TAXONDB.out.ready.map { params.taxondb }
    db      = SETUP_FUNANNOTATE_DB.out.db
    config  = SETUP_AUGUSTUS_CONFIG.out.config
}
