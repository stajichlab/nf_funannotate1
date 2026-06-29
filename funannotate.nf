#!/usr/bin/env nextflow

/*
 * SOURCE: ../../../1KFG/common_annotate/pipeline/nextflow/funannotate.nf
 * Last synced: 2026-05-23
 * Changes vs source: removed nextflow.enable.dsl=2; params block moved to
 *                    conf/profile_annotate.config.
 *
 * Usage (from project root — a pipeline profile is REQUIRED; without it
 * params.taxondb / params.funannotate_db are null and parsing fails):
 *   sbatch nextflow/run_annotate.sh
 *   nextflow run nextflow/funannotate.nf -c nextflow/nextflow.config \
 *       -profile annotate,slurm,ucr_hpcc -resume
 */

// Data contract: every channel element is `tuple val(meta), val/path(genome)`.
// meta is a Map built by SampleUtils.makeMeta(row) — see lib/SampleUtils.groovy.
//   meta.id is the ONLY field used for tag{} and file naming.
//   meta.asmid, meta.species, meta.strain, meta.locustag, meta.busco,
//   meta.transl_table, meta.taxonid carry payload used inside process scripts.
//   header_length is NOT in meta — it comes from params.header_length (default 24).



// staleRnaseq wraps FunannotateUtils.staleRnaseq to add a log.info (requires DSL scope).
def staleRnaseq(String out, String species) {
    if (FunannotateUtils.staleRnaseq(out, species, params.target as String, launchDir.toString())) {
        log.info "stale prediction for ${out}: rnaseq/trinity newer than GBK — scheduling retrain+repredict"
        return true
    }
    return false
}

include { validateParameters; paramsSummaryLog; paramsHelp } from 'plugin/nf-schema'
include { ASM_STATS }        from './modules/local/asm_stats'
include { INPUT_CHECK }      from './subworkflows/local/input_check'
include { SETUP_DBS }        from './subworkflows/local/setup_dbs'
include { CLEAN_GENOMES }    from './subworkflows/local/clean_genomes'
include { MASK_GENOME }      from './subworkflows/local/mask_genome'
include { FETCH_RNASEQ }     from './subworkflows/local/fetch_rnaseq'
include { TRAIN_PREDICT }    from './subworkflows/local/train_predict'
include { ANNOTATE_GENOME }  from './subworkflows/local/annotate_genome'

workflow {
    // `--help` prints schema-driven parameter help (grouped, with types/defaults) and exits.
    if (params.help) {
        log.info paramsHelp()
        exit 0
    }
    // Type-check params against nextflow_schema.json and log the resolved set.
    // (Unrecognised params warn rather than fail — see nextflow.config.)
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Fail fast with an actionable message when a pipeline profile was not selected
    // (these params come from conf/profile_annotate.config). Without it, downstream
    // file(params.funannotate_db) calls throw a cryptic "file() ... cannot be null".
    if( !params.taxondb || !params.funannotate_db )
        error "Missing params.taxondb / params.funannotate_db — add a pipeline profile, e.g. -profile annotate,slurm,module (or use: sbatch nextflow/run_annotate.sh)"

    // ── Samplesheet ingestion (INPUT_CHECK) ──────────────────────────────────
    // Parses samples CSV, applies taxon/asmid/suppress/n_test filters, builds
    // meta maps, and resolves genome paths. Two outputs:
    //   jobs        — tuple(meta, gz)  with genome existence filter (cleaning path)
    //   postpredict — meta only        no genome filter (annotate/update paths)
    INPUT_CHECK()
    def jobs = INPUT_CHECK.out.genomes

    def ch_versions = Channel.empty()
    if (params.debug.toBoolean()) {
        jobs.view { meta, gz -> "[CHANNEL] Submitting: out=${meta.id}, asmid=${meta.asmid}, transl_table=${meta.transl_table}, gz=${gz}" }
    }

    // Build/seed the three run-once databases. All use storeDir so they are no-ops
    // on any run where their target directories already exist.
    SETUP_DBS()
    def taxondb_ch = SETUP_DBS.out.taxondb

    CLEAN_GENOMES(jobs, taxondb_ch)

    if (!params.only_clean.toBoolean()) {
        def clean_genome_ch = CLEAN_GENOMES.out.genomes

        // ── Generate assembly statistics (for earlgrey_mask.nf SELECT_REPS) ────────
        // Generate asm_stats.tsv if --gen_asm_stats is true and the file doesn't exist.
        // This is used by earlgrey_mask.nf to select representative genomes per species.
        if (params.gen_asm_stats.toBoolean()) {
            def asm_stats_path = file(params.tables_dir).toAbsolutePath()
            def asm_stats_gz = file("${asm_stats_path}/asm_stats.tsv.gz")
            if (!asm_stats_gz.exists()) {
                log.info "Generating assembly statistics: ${asm_stats_gz}"
                ASM_STATS(
                    file(params.samples),
                    file(params.genome_dir)
                )
                ch_versions = ch_versions.mix(ASM_STATS.out.versions)
            } else {
                log.info "Assembly statistics already exist: ${asm_stats_gz}"
            }
        }

        // ── Repeat masking ────────────────────────────────────────────────────────
        // predict_genome_ch: tantan soft-masked (default) or clean/prior-masked genome.
        // MASK_GENOME handles the run_repeatmasker if/else and storeDir-cached masking.
        MASK_GENOME(clean_genome_ch)
        def predict_genome_ch = MASK_GENOME.out.genomes

        // Gate the predict chain on funannotate DB + augustus config being ready.
        // SETUP_DBS was already called above; its storeDir-cached outputs are free
        // on resumed runs. Gating here threads the dependency through the entire
        // downstream funannotate subgraph (train, predict, update, annotate).
        // (MASKREPEAT uses `funannotate mask`, which needs neither, so it is intentionally
        // left ungated and can run in parallel with these setup steps.)
        predict_genome_ch = predict_genome_ch
            .combine(SETUP_DBS.out.db)
            .combine(SETUP_DBS.out.config)
            .map { row -> row[0..-3] }

        // SRA read fetching (FETCH_RNASEQ) + training (RNASEQ_PREPARE + FUNANNOTATE_TRAIN)
        // + prediction (FUNANNOTATE_PREDICT). All three are composed in TRAIN_PREDICT.
        // When run_sra_fetch=false, reads channel is empty and genomes go straight to predict.
        def reads_ch = Channel.empty()
        if (params.run_sra_fetch.toBoolean()) {
            def sra_input = predict_genome_ch
                .map { meta, _genome_fa ->
                    def species_tag = meta.species.replaceAll(/\s+/, '_')
                    tuple(species_tag, meta.taxonid)
                }
                .groupTuple(by: 0)
                .map { species_tag, taxonids -> tuple(species_tag, taxonids[0]) }

            FETCH_RNASEQ(sra_input)
            reads_ch = FETCH_RNASEQ.out.reads
        }

        if (!params.run_sra_fetch.toBoolean() || (!params.stop_after_sra_fetch.toBoolean() && !params.stop_after_sra_query.toBoolean())) {
        TRAIN_PREDICT(predict_genome_ch, reads_ch)

        // Genomes predicted in a PRIOR run (complement of what TRAIN_PREDICT runs this run).
        def postpredict = INPUT_CHECK.out.samples
            .filter { meta ->
                FunannotateUtils.gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) != null &&
                !staleRnaseq(meta.id as String, meta.species as String)
            }
        // Merge prior-run and same-run predictions; feed into annotation chain.
        def predict_meta = postpredict.mix(TRAIN_PREDICT.out.metadata)

        ANNOTATE_GENOME(predict_meta, reads_ch)
        ch_versions = ch_versions.mix(ANNOTATE_GENOME.out.versions)
        } // end if (!params.stop_after_sra_fetch || !params.run_sra_fetch)
    }

    // Collect software versions from all processes that emit versions.yml.
    // Written to logs/software_versions.yml alongside the trace file.
    ch_versions
        .unique()
        .collectFile(
            name:     'software_versions.yml',
            storeDir: "${launchDir}/logs/nextflow",
            newLine:  true
        )
}

