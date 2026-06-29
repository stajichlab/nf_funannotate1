/*
 * ANNOTATE_GENOME — optional pre-annotate steps + funannotate update + funannotate annotate.
 *
 * Accepts all assemblies with a completed predict GBK (both prior-run and same-run).
 * Chains optional tools (ANTISMASH, INTERPRO, SIGNALP) before annotation; each tool
 * splits the channel into todo/done, processes todo, then re-mixes so FUNANNOTATE_ANNOTATE
 * only fires once all requested tools are complete.
 *
 * Takes:
 *   ch_predict_meta — val(meta) for every assembly ready for annotation
 *   ch_reads        — tuple(val(species_tag), path(r1), path(r2), path(se))
 *                     May be empty when run_sra_fetch=false; only used by UPDATE.
 *
 * Emits:
 *   versions — mixed versions.yml paths from ANTISMASH / INTERPRO / SIGNALP
 */

include { ANTISMASH_RUN    } from './../../modules/local/antismash_run'
include { INTERPROSCAN_RUN } from './../../modules/local/interproscan_run'
include { SIGNALP_RUN      } from './../../modules/local/signalp_run'
include { FUNANNOTATE_UPDATE  } from './../../modules/local/funannotate_update'
include { FUNANNOTATE_ANNOTATE } from './../../modules/local/funannotate_annotate'

// GBK may be stored compressed to save space; return the first non-empty file or null.
def _gbkResult(String dir, String id) {
    def plain = file("${dir}/${id}.gbk")
    if (plain.exists() && plain.size() > 0) return plain
    def gz = file("${dir}/${id}.gbk.gz")
    if (gz.exists() && gz.size() > 0) return gz
    return null
}

workflow ANNOTATE_GENOME {

    take:
    ch_predict_meta  // val(meta)
    ch_reads         // tuple(val(species_tag), path(r1), path(r2), path(se)) — may be empty

    main:
    def ch_versions       = Channel.empty()
    def annotate_ready_ch = ch_predict_meta

    if (params.run_antismash.toBoolean()) {
        def as_todo = annotate_ready_ch.filter { meta ->
            def asDir = file("${params.target}/${meta.id}/antismash_local")
            !(asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') })
        }
        def as_done = annotate_ready_ch.filter { meta ->
            def asDir = file("${params.target}/${meta.id}/antismash_local")
            asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') }
        }
        ANTISMASH_RUN(as_todo)
        ch_versions = ch_versions.mix(ANTISMASH_RUN.out.versions)
        annotate_ready_ch = ANTISMASH_RUN.out.results.map { meta, _files -> meta }.mix(as_done)
    }

    if (params.run_interpro.toBoolean()) {
        def ipr_todo = annotate_ready_ch.filter { meta ->
            !file("${params.target}/${meta.id}/annotate_misc/iprscan.xml").exists()
        }
        def ipr_done = annotate_ready_ch.filter { meta ->
            file("${params.target}/${meta.id}/annotate_misc/iprscan.xml").exists()
        }
        INTERPROSCAN_RUN(ipr_todo)
        ch_versions = ch_versions.mix(INTERPROSCAN_RUN.out.versions)
        annotate_ready_ch = INTERPROSCAN_RUN.out.results.map { meta, _xml -> meta }.mix(ipr_done)
    }

    if (params.run_signalp.toBoolean()) {
        def sp_todo = annotate_ready_ch.filter { meta ->
            !file("${params.target}/${meta.id}/annotate_misc/signalp.results.txt").exists()
        }
        def sp_done = annotate_ready_ch.filter { meta ->
            file("${params.target}/${meta.id}/annotate_misc/signalp.results.txt").exists()
        }
        SIGNALP_RUN(sp_todo)
        ch_versions = ch_versions.mix(SIGNALP_RUN.out.versions)
        annotate_ready_ch = SIGNALP_RUN.out.results.map { meta, _txt -> meta }.mix(sp_done)
    }

    if (params.run_update.toBoolean()) {
        if (params.run_sra_fetch.toBoolean()) {
            def upd_input = ch_predict_meta
                .map { meta ->
                    def species_tag = meta.species.replaceAll(/\s+/, '_')
                    tuple(species_tag, meta)
                }
                .combine(ch_reads, by: 0)
                .map { _st, meta, r1, r2, _se -> tuple(meta, r1, r2) }
            def upd_todo = upd_input.filter { meta, _r1, _r2 ->
                _gbkResult("${params.target}/${meta.id}/update_results", meta.id as String) == null
            }
            def upd_done_signal = upd_input
                .filter { meta, _r1, _r2 ->
                    _gbkResult("${params.target}/${meta.id}/update_results", meta.id as String) != null
                }
                .map { meta, _r1, _r2 -> tuple(meta.id, 'upd') }
            FUNANNOTATE_UPDATE(upd_todo)
            def upd_signal = FUNANNOTATE_UPDATE.out
                .map { meta -> tuple(meta.id, 'upd') }
                .mix(upd_done_signal)
            annotate_ready_ch = annotate_ready_ch
                .map { meta -> tuple(meta.id, meta) }
                .join(upd_signal)
                .map { _id, meta, _flag -> meta }
        } else {
            log.warn "run_update=true but run_sra_fetch=false; funannotate update skipped (no reads available)"
        }
    }

    if (params.run_annotate.toBoolean()) {
        FUNANNOTATE_ANNOTATE(annotate_ready_ch.filter { meta ->
            _gbkResult("${params.target}/${meta.id}/annotate_results", meta.id as String) == null
        })
    }

    emit:
    versions = ch_versions
}
