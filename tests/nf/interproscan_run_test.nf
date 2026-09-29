// Test harness: run INTERPROSCAN_RUN alone, under the real pipeline configs,
// on genomes that already have predict_results/<id>.proteins.fa under
// params.target. Driven by tests/test_interproscan6.sh.
include { INTERPROSCAN_RUN } from '../../modules/local/interproscan_run'

// Same wrapper name as the real pipeline, so the task is named
// ANNOTATE_GENOME:INTERPROSCAN_RUN and the '.*:INTERPROSCAN_RUN' resource /
// queue selectors in conf/ apply exactly as in production.
workflow ANNOTATE_GENOME {
    take:
    ch_meta

    main:
    INTERPROSCAN_RUN(ch_meta)

    emit:
    versions = INTERPROSCAN_RUN.out.versions
}

workflow {
    def ids = (params.test_ids as String).tokenize(',')
    ANNOTATE_GENOME(channel.fromList(ids).map { id -> [id: id] })
    ANNOTATE_GENOME.out.versions.view { v -> v.text }
}
