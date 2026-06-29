/*
 * CLEAN_GENOMES — FCS-GX contaminant purge + length-filter
 *
 * Accepts the full jobs channel and the taxondb gate channel from SETUP_DBS.
 * Handles the batching decision internally (params.clean_batch_size > 0 and
 * !params.skip_fcs → GENOME_CLEAN_BATCH; otherwise per-genome GENOME_CLEAN).
 *
 * Emits:
 *   genomes — channel: tuple(val(meta), val(genome_fa))
 *             genome_fa is an absolute-path String to the cleaned .fa/.fa.gz;
 *             samples with no cleaned file are silently dropped with a warning.
 */

include { GENOME_CLEAN       } from './../../modules/local/genome_clean'
include { GENOME_CLEAN_BATCH } from './../../modules/local/genome_clean_batch'

workflow CLEAN_GENOMES {

    take:
    jobs    // channel: tuple(val(meta), path(gz))
    taxondb // channel: val(String taxondb_path)

    main:
    // Only submit cleaning jobs for genomes not already in input_clean_genomes/.
    // Same logic as FunannotateUtils.genomeFile(): prefer .fa.gz if non-empty, else .fa.
    // GENOME_CLEAN_BATCH also re-checks per genome at runtime for partial retry.
    def jobs_to_clean = jobs.filter { meta, gz ->
        def fgz = file("${launchDir}/input_clean_genomes/${meta.asmid}.fa.gz")
        def fa  = file("${launchDir}/input_clean_genomes/${meta.asmid}.fa")
        !(fgz.exists() && fgz.size() > 0) && !fa.exists()
    }

    def clean_done_ch
    int clean_batch_size = params.clean_batch_size as int
    if (clean_batch_size > 0 && !params.skip_fcs.toBoolean()) {
        // Wrap each collated batch in a single-element list so .combine() appends
        // taxondb as the 2nd tuple element rather than spreading the batch rows.
        def clean_batches = jobs_to_clean.collate(clean_batch_size).map { batch -> [ batch ] }
        GENOME_CLEAN_BATCH(clean_batches.combine(taxondb))
        clean_done_ch = GENOME_CLEAN_BATCH.out.manifest.collect().ifEmpty([])
    } else {
        GENOME_CLEAN(jobs_to_clean.combine(taxondb))
        clean_done_ch = GENOME_CLEAN.out.genome.map { meta, gz -> gz }.collect().ifEmpty([])
    }

    // Re-attach the cleaned genome path to each sample's meta. Gate on clean_done_ch
    // so this channel only opens once all cleaning tasks have finished.
    // Use it[0] index access (not fixed-arity destructuring) because combine with
    // ifEmpty([]) produces a variable-length tuple (2 elements when the sentinel is
    // empty, 3+ when it carries collected paths).
    // Same logic as FunannotateUtils.genomeFile(): prefer .fa.gz if non-empty, else .fa.
    def ch_genomes = jobs
        .combine(clean_done_ch)
        .map { row ->
            def meta = row[0]
            def fgz  = file("${launchDir}/input_clean_genomes/${meta.asmid}.fa.gz")
            def fa   = file("${launchDir}/input_clean_genomes/${meta.asmid}.fa")
            def g    = (fgz.exists() && fgz.size() > 0) ? fgz : fa
            tuple(meta, g)
        }
        .filter { meta, g ->
            if (!g.exists()) {
                log.warn "No cleaned genome for ${meta.id} (asmid=${meta.asmid}) — skipping downstream"
                return false
            }
            return true
        }
        .map { meta, genome_fa ->
            tuple(meta, genome_fa.toAbsolutePath().toString())
        }

    emit:
    genomes = ch_genomes
}
