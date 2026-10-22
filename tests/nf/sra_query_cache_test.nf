// Test harness: run SRA_QUERY_BATCH alone for one species (params.species_tag,
// params.taxonid). Driven by tests/test_sra_query_cache.sh, which puts fake
// esearch/efetch on PATH.
include { SRA_QUERY_BATCH } from '../../modules/local/sra_query_batch'

params.species_tag      = 'Aspergillus_nidulans'
params.taxonid          = '162425'
params.enable_single_end = false
params.max_rnaseq_se_runs = 3
params.species_list     = ''

// params.species_list = 'Tag1:id+id,Tag2:id' runs several species as one batch (live test).
workflow {
    def tags = [params.species_tag]
    def ids  = [params.taxonid]
    if (params.species_list) {
        def pairs = params.species_list.toString().split(',').collect { it.split(':', 2) }
        tags = pairs.collect { it[0] }
        ids  = pairs.collect { it[1] }
    }
    SRA_QUERY_BATCH(Channel.of(tuple(tags, ids)))
}
