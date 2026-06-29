/*
 * MASK_GENOME — tantan soft-masking with run_repeatmasker fallback
 *
 * When params.run_repeatmasker is true: runs MASKREPEAT_TANTAN_RUN (storeDir-cached).
 * When false: reuses the masked genome from a prior run if it exists, otherwise
 * passes the clean genome through unmasked.
 *
 * Emits:
 *   genomes — channel: tuple(val(meta), val(genome_fa))
 *             genome_fa is an absolute-path String consumed by the predict chain.
 *
 * Note: genomeFile() is a top-level def in funannotate.nf, not accessible here.
 * The fallback path inlines the equivalent logic (prefer .masked.fasta.gz if
 * non-empty, else .masked.fasta).
 */

include { MASKREPEAT_TANTAN_RUN } from './../../modules/local/maskrepeat_tantan_run'

workflow MASK_GENOME {

    take:
    ch_clean  // channel: tuple(val(meta), val(genome_fa))

    main:
    def ch_predict
    if (params.run_repeatmasker.toBoolean()) {
        MASKREPEAT_TANTAN_RUN(ch_clean)
        ch_predict = MASKREPEAT_TANTAN_RUN.out.masked
            .map { meta, masked_fa -> tuple(meta, masked_fa.toAbsolutePath().toString()) }
    } else {
        // --run_repeatmasker false: use masked genome from a prior run if available.
        ch_predict = ch_clean
            .map { meta, genome_fa ->
                def fgz    = file("${launchDir}/input_clean_genomes/${meta.asmid}.masked.fasta.gz")
                def fa     = file("${launchDir}/input_clean_genomes/${meta.asmid}.masked.fasta")
                def masked = (fgz.exists() && fgz.size() > 0) ? fgz : fa
                def use_fa = masked.exists() ? masked.toString() : genome_fa
                if (params.debug.toBoolean()) {
                    log.info "[DEBUG] ${meta.asmid}: genome_fa=${use_fa} (masked=${masked.exists()})"
                }
                tuple(meta, use_fa)
            }
    }

    emit:
    genomes = ch_predict
}
