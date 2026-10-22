// Test harness: run SRA_QUERY_BATCH alone for one species (params.species_tag,
// params.taxonid). Driven by tests/test_sra_query_cache.sh, which puts fake
// esearch/efetch on PATH.
include { SRA_QUERY_BATCH } from '../../modules/local/sra_query_batch'

params.species_tag      = 'Aspergillus_nidulans'
params.taxonid          = '162425'
params.enable_single_end = false
params.max_rnaseq_se_runs = 3

workflow {
    SRA_QUERY_BATCH(Channel.of(tuple([params.species_tag], [params.taxonid])))
}
