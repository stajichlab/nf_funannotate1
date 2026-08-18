/*
 * PICK_REPRESENTATIVE_STRAIN — select one representative strain per species
 * (BUSCO completeness, N50 tiebreak, alphabetical final), restricted to strains
 * with ANI coverage, and write abinitio_reuse_assignments.csv + repr_assignments.tsv.
 *
 * Port of BFD/Fungi_BFD/nextflow/modules/ani/report/PICK_REPRESENTATIVE_STRAIN/
 * main.nf. BFD read merged Parquet tables here; this port reads the raw
 * BUSCO_GENOME *.BUSCO_summary.*.txt reports + ASM_STATS's asm_stats.tsv.gz
 * directly (see bin/pick_representative_strain.py), so no duckdb/pip install
 * is needed — plain python3 from the base provisioning env.
 *
 * Option B persistence: writes directly into ${params.ani_outdir}/_reuse_assignments/
 * (per-species + merged CSVs merge-lactically via flock), publishes copies.
 * Backfills the shared ab-initio store when the representative already has a
 * prediction on disk (see bin/backfill_abinitio_params.py).
 */
process PICK_REPRESENTATIVE_STRAIN {
    tag   "PICK_REPRESENTATIVE_STRAIN"
    label 'report'

    publishDir "${params.ani_outdir}/_reuse_assignments", mode: 'copy', pattern: 'abinitio_reuse_assignments*.csv'
    publishDir "${params.ani_outdir}/_reuse_assignments", mode: 'copy', pattern: 'repr_assignments.tsv'

    input:
    path ani_tsv
    path predict_input_tsv
    path samples_csv
    path busco_summaries
    path asm_stats

    output:
    path("abinitio_reuse_assignments.csv"), emit: outCSV
    path("repr_assignments.tsv"),           emit: outAssignments

    script:
    def threshold   = params.ani_reuse_threshold ?: 99.0
    def shared_root = params.gene_prediction_shared_abinitio ?: ''
    def target      = params.target ?: ''
    def aug_cfg     = params.augustus_config ?: ''
    """
    mkdir -p _reuse_assignments
    python "${projectDir}/bin/pick_representative_strain.py" \
        --ani-tsv "${ani_tsv}" \
        --predict-input "${predict_input_tsv}" \
        --samples "${samples_csv}" \
        --busco-summaries ${busco_summaries} \
        --asm-stats "${asm_stats}" \
        --out-dir "_reuse_assignments" \
        --ani-threshold "${threshold}" \
        ${shared_root ? "--shared-root ${shared_root}" : ''} \
        ${target ? "--target ${target}" : ''} \
        ${aug_cfg ? "--augustus-config ${aug_cfg}" : ''}
    mv _reuse_assignments/abinitio_reuse_assignments*.csv .
    mv _reuse_assignments/repr_assignments.tsv .
    """

    stub:
    """
    printf 'species\\tout\\tis_representative\\trepresentative_out\\tani_to_representative\\treuse_eligible\\n' > abinitio_reuse_assignments.csv
    printf '' > repr_assignments.tsv
    """
}
