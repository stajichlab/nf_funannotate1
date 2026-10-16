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

include { ANTISMASH_RUN        } from './../../modules/local/antismash_run'
include { SETUP_ANTISMASH_DB   } from './../../modules/local/setup_antismash_db'
include { INTERPROSCAN_RUN     } from './../../modules/local/interproscan_run'
include { INTERPROSCAN5_RUN    } from './../../modules/local/interproscan5_run'
include { SIGNALP_RUN          } from './../../modules/local/signalp_run'
include { DEEPTMHMM_ANNOTATION } from './../../modules/local/deeptmhmm_annotation'
include { FUNANNOTATE_UPDATE   } from './../../modules/local/funannotate_update'
include { FUNANNOTATE_ANNOTATE } from './../../modules/local/funannotate_annotate'

// funannotate annotate gzips a plain iprscan.xml after parsing it and deletes
// the original, and INTERPROSCAN_RUN (IPS6) writes the .gz directly -- so
// either form means InterProScan is done for this assembly. A result is reused only if it is
// a real file (>= 100 bytes; the stub writes a 30-byte placeholder) and, when INTERPROSCAN_RUN
// recorded the md5 of the proteins it scanned (annotate_misc/iprscan.proteins.md5), that md5
// still matches predict_results/<id>.proteins.fa. Results without a record are reused.
def iprscanStatus(meta) {
    def misc = "${params.target}/${meta.id}/annotate_misc"
    def xml  = ['iprscan.xml', 'iprscan.xml.gz'].collect { new File("${misc}/${it}") }.find { it.exists() }
    return FunannotateUtils.reuseStatus(
        xml,
        new File("${params.target}/${meta.id}/predict_results/${meta.id}.proteins.fa"),
        new File("${misc}/iprscan.proteins.md5"),
        100L)
}
def iprscanDone(meta) { return FunannotateUtils.reusable(iprscanStatus(meta)) }

// antiSMASH: the JSON in antismash_local/, checked against the predict GBK it ran on
// (antismash_local/antismash.input.md5, written by ANTISMASH_RUN).
def antismashStatus(meta) {
    def dir = new File("${params.target}/${meta.id}/antismash_local")
    def json = dir.isDirectory() ? dir.listFiles().find { it.name.endsWith('.json') || it.name.endsWith('.json.gz') } : null
    def predictDir = "${params.target}/${meta.id}/predict_results"
    def gbk = FunannotateUtils.gbkResult(predictDir, meta.id as String)
    return FunannotateUtils.reuseStatus(json, gbk, new File(dir, 'antismash.input.md5'), 100L)
}
def antismashDone(meta) { return FunannotateUtils.reusable(antismashStatus(meta)) }

// Log why a result is being regenerated (only for results that exist but cannot be reused).
def logRerun(String what, meta, String status) {
    if (status == 'stale') log.info "${meta.id}: ${what} result is stale (input changed since it was computed); rerunning"
    else if (status == 'tiny') log.info "${meta.id}: ${what} result is empty or a placeholder; rerunning"
}

workflow ANNOTATE_GENOME {

    take:
    ch_predict_meta  // val(meta)
    ch_reads         // tuple(val(species_tag), path(r1), path(r2), path(se)) — may be empty

    main:
    def ch_versions       = Channel.empty()
    def annotate_ready_ch = ch_predict_meta

    if (params.run_antismash.toBoolean()) {
        SETUP_ANTISMASH_DB()
        def as_todo = annotate_ready_ch.filter { meta ->
            def st = antismashStatus(meta)
            logRerun('antiSMASH', meta, st)
            !FunannotateUtils.reusable(st)
        }
        def as_done = annotate_ready_ch.filter { meta -> antismashDone(meta) }
        // Gate on the DB being ready (storeDir means this is a no-op after the
        // first run) -- .combine() forces ANTISMASH_RUN to wait rather than
        // racing SETUP_ANTISMASH_DB on a cold cache.
        def as_todo_gated = as_todo.combine(SETUP_ANTISMASH_DB.out.db).map { meta, _db -> meta }
        ANTISMASH_RUN(as_todo_gated)
        ch_versions = ch_versions.mix(ANTISMASH_RUN.out.versions)
        annotate_ready_ch = ANTISMASH_RUN.out.results.map { meta, _files -> meta }.mix(as_done)
    }

    if (params.run_interpro.toBoolean()) {
        def ipr_todo = annotate_ready_ch.filter { meta ->
            def st = iprscanStatus(meta)
            logRerun('InterProScan', meta, st)
            !FunannotateUtils.reusable(st)
        }
        def ipr_done = annotate_ready_ch.filter { meta -> iprscanDone(meta) }
        if (params.interproscan_engine == 'ips5') {
            INTERPROSCAN5_RUN(ipr_todo)
            ch_versions = ch_versions.mix(INTERPROSCAN5_RUN.out.versions)
            annotate_ready_ch = INTERPROSCAN5_RUN.out.results.map { meta, _xml -> meta }.mix(ipr_done)
        } else {
            INTERPROSCAN_RUN(ipr_todo)
            ch_versions = ch_versions.mix(INTERPROSCAN_RUN.out.versions)
            annotate_ready_ch = INTERPROSCAN_RUN.out.results.map { meta, _done -> meta }.mix(ipr_done)
        }
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

    if (params.run_deeptmhmm.toBoolean()) {
        def tmhmm_todo = annotate_ready_ch.filter { meta ->
            !file("${params.target}/${meta.id}/annotate_misc/TMRs.gff3").exists()
        }
        def tmhmm_done = annotate_ready_ch.filter { meta ->
            file("${params.target}/${meta.id}/annotate_misc/TMRs.gff3").exists()
        }
        DEEPTMHMM_ANNOTATION(tmhmm_todo)
        ch_versions = ch_versions.mix(DEEPTMHMM_ANNOTATION.out.versions)
        annotate_ready_ch = DEEPTMHMM_ANNOTATION.out.results.map { meta, _gff3 -> meta }.mix(tmhmm_done)
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
                FunannotateUtils.gbkResult("${params.target}/${meta.id}/update_results", meta.id as String) == null
            }
            def upd_done_signal = upd_input
                .filter { meta, _r1, _r2 ->
                    FunannotateUtils.gbkResult("${params.target}/${meta.id}/update_results", meta.id as String) != null
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
        // SBT per genome: the samples.csv `SBT` column (a path, or a file name under
        // --sbt_dir), else --sbt_template, else none (funannotate's own template).
        // --sbt is an alias for --sbt_template (--sbt wins when both are given).
        // An unusable SBT entry or a BioSample mismatch only warns, since an NCBI
        // deposit is not required for these runs; --sbt_strict turns both into errors.
        def strict      = params.sbt_strict.toBoolean()
        def sbt_file    = (params.sbt ?: params.sbt_template) as String
        def sbt_default = sbt_file && file(sbt_file).exists() ? file(sbt_file).toRealPath() : null
        if (!sbt_default) {
            log.warn "SBT template '${sbt_file}' not found; genomes without a usable SBT column entry use funannotate's default template"
        }
        def sbt_by_id = channel.fromPath(params.samples)
            .splitCsv(header: true)
            .map { row -> tuple(SampleUtils.makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: ''),
                                row.SBT?.trim() ?: '', row.BIOSAMPLE?.trim() ?: '') }
        def annotate_in = annotate_ready_ch
            .filter { meta ->
                FunannotateUtils.gbkResult("${params.target}/${meta.id}/annotate_results", meta.id as String) == null
            }
            .map { meta -> tuple(meta.id, meta) }
            .join(sbt_by_id)
            .map { _id, meta, sbt_name, biosample ->
                def sbt = SampleUtils.resolveSbt(sbt_name, params.sbt_dir as String, launchDir.toString())
                if (sbt_name && !sbt) {
                    def msg = "${meta.id}: SBT '${sbt_name}' not found (sbt_dir=${params.sbt_dir ?: launchDir})"
                    if (strict) error "${msg}; --sbt_strict is set"
                    log.warn "${msg}; using ${sbt_default ?: 'funannotate default template'}"
                }
                def chosen = sbt ?: sbt_default
                // The SBT must name the same BioSample as the sheet's BIOSAMPLE column.
                // A template that names none (the generic default) is not checked.
                if (chosen && biosample) {
                    def named = SampleUtils.sbtBiosamples(chosen)
                    if (named && !named.contains(biosample)) {
                        def msg = "${meta.id}: SBT ${chosen.name} names BioSample ${named.join(',')} but the sheet has ${biosample}"
                        if (strict) error "${msg}; --sbt_strict is set"
                        log.warn msg
                    }
                }
                tuple(meta, chosen ?: file("${projectDir}/assets/NO_FILE"))
            }
        FUNANNOTATE_ANNOTATE(annotate_in)
    }

    emit:
    versions = ch_versions
}
