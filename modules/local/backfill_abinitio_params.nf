/*
 * BACKFILL_ABINITIO_PARAMS — copy representative strains' trained AUGUSTUS/
 * GeneMark/SNAP ab-initio parameters into the shared per-species parameter
 * store (params.gene_prediction_shared_abinitio), so sibling strains can reuse
 * them via `funannotate predict -p parameters.json`.
 *
 * Port of BFD/Fungi_BFD/nextflow/modules/funannotate/predict/
 * BACKFILL_ABINITIO_PARAMS/main.nf. Actual logic lives in
 * bin/backfill_abinitio_params.py (staging dir + content-hash idempotency +
 * atomic swap). Batched, not one task per representative: callers group
 * candidates with .collate(100), so one job loops over up to ~100 species via
 * --manifest. Each line is independent and idempotent (content-hash
 * short-circuit), so a retry/-resume just re-confirms entries as up to date.
 *
 * Option B persistence model (same as FUNANNOTATE_PREDICT): the real output is
 * written directly into the persistent shared-params store, not the work dir.
 * No publishDir copy. The emitted marker only carries the DAG edge; callers
 * re-derive the actual path from disk via sharedParamsJsonFor() once this task
 * completes (the write happens synchronously before the marker is touched).
 */
process BACKFILL_ABINITIO_PARAMS {
    tag   "batch_$batch_id"
    label 'report'

    input:
    // items: List of [species, rep_out, genemark_mod] -- 3rd field is ''
    // when GENEMARK_RUN is off (run_genemark=false) or the caller is the
    // standalone backfill sweep. backfill_abinitio_params.py falls back to
    // deriving genemark.mod from predict_misc/ when the field is empty.
    tuple val(batch_id), val(items)

    output:
    tuple val(items), path("*.backfill.done"), emit: done

    script:
    def shared_root  = params.gene_prediction_shared_abinitio
    def target       = params.target
    def threshold    = params.ani_reuse_threshold ?: 99.0
    def aug_cfg      = params.augustus_config ?: ''
    // Manifest content is fully resolved in Groovy before the shell sees it,
    // so there's no bare `$` in this string for Groovy to misinterpolate.
    def manifest = items.collect { sp, out, mod -> "${sp}\t${out}\t${mod ?: ''}" }.join('\n')
    """
    cat > manifest.tsv <<'BACKFILL_MANIFEST_EOF'
${manifest}
BACKFILL_MANIFEST_EOF
    python "${projectDir}/bin/backfill_abinitio_params.py" \
        --manifest manifest.tsv \
        --target "${target}" \
        --shared-root "${shared_root}" \
        --ani-threshold "${threshold}" \
        ${aug_cfg ? "--augustus-config ${aug_cfg}" : ''}
    touch "batch_${batch_id}.backfill.done"
    """

    stub:
    """
    touch "batch_${batch_id}.backfill.done"
    """
}
