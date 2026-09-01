/*
 * PREDICT_REUSE — representative-gated FUNANNOTATE_PREDICT/GENEMARK_RUN with
 * per-species ab-initio reuse (port of BFD/Fungi_BFD/nextflow/subworkflows/
 * local/FUNANNOTATE_PREDICTION.nf, expressed with this repo's (meta, genome_fa)
 * channel contract).
 *
 * Only invoked from TRAIN_PREDICT when --run_ani_reuse=true (see the validation
 * in funannotate.nf: run_ani_reuse requires --gene_prediction_shared_abinitio,
 * the shared per-species ab-initio store). The reuse map is an
 * abinitio_reuse_assignments.csv produced in-run by ANI_REUSE; its `out` column
 * == meta.id. Strains with no row (singleton species, or ANI below threshold)
 * are independents and train ungated, exactly as before this feature existed.
 *
 * Branching (BFD logic, unchanged):
 *   representative   — is_representative=true. Never reuses a shared .mod (they
 *                      DEFINE the store): fresh GeneMark train + predict, then
 *                      BACKFILL_ABINITIO_PARAMS copies their freshly-trained
 *                      .mod/parameters into the shared store once predict lands.
 *   eligible_sibling — a non-rep strain whose ANI-to-rep >= threshold. Only
 *                      proceeds once its species' shared params are available —
 *                      pre-existing on disk, or freshly backfilled by THIS run's
 *                      own representative. GENEMARK_RUN_SIB fast-reuses the
 *                      shared species .mod (--predict_with) when present.
 *   independent      — everything else: fresh train, no shared anything.
 *
 * For eligible siblings whose species has no shared store by the time the run
 * needs them (blocked): written to
 * params.target/predict_blocked_awaiting_representative.tsv and excluded —
 * unless --allow_independent_fallback, in which case they train independently
 * anyway (still logged).
 */

include { FUNANNOTATE_PREDICT }                 from './../../modules/local/funannotate_predict'
include { FUNANNOTATE_PREDICT as FUNANNOTATE_PREDICT_SIB } from './../../modules/local/funannotate_predict'
include { GENEMARK_RUN }                        from './../../modules/local/genemark_run'
include { GENEMARK_RUN as GENEMARK_RUN_SIB }    from './../../modules/local/genemark_run'
include { PRODIGAL_RUN }                        from './../../modules/local/prodigal_run'
include { PRODIGAL_RUN as PRODIGAL_RUN_SIB }    from './../../modules/local/prodigal_run'
include { BACKFILL_ABINITIO_PARAMS }            from './../../modules/local/backfill_abinitio_params'

workflow PREDICT_REUSE {

    take:
    predict_input_ch   // tuple(val(meta), path(genome_fa)) — assemblies to potentially predict
    reuse_map_csv      // val(path) — ANI_REUSE's abinitio_reuse_assignments.csv

    main:
    def allowFallback = (params.allow_independent_fallback ?: false).toString().toBoolean()
    def runGenemark   = (params.run_genemark ?: true).toString().toBoolean()
    def genemarkMode  = ((params.genemark_mode ?: 'ES') as String).toUpperCase()
    def sharedRoot    = params.gene_prediction_shared_abinitio ?: ''

    // BFD's forceIndependentSet / forceIndependentGenemarkSet are per-species /
    // per-out lists; this repo has a single params.force_independent applying to
    // every genome. forceAll => no species eligible AND every GeneMark retrains.
    def forceAll = (params.force_independent ?: false).toString().toBoolean()

    // ── Classify each strain once, up front ──────────────────────────────────
    // Load the reuse map OFFLINE (plain synchronous file read, per BFD) into a
    // single-emission channel, then broadcast it to every predict row via
    // .combine() — which also makes prediction WAIT for the in-run ANI
    // front-end (PICK_REPRESENTATIVE_STRAIN) to land the file, so the read
    // below is never raced. An empty/missing map is NOT silently treated as
    // "everything independent": BFD semantics are to fail early.
    def reuseMapPerItem = reuse_map_csv
        .map { csvPath ->
            def p = csvPath as String
            if (!p)
                error "PREDICT_REUSE: empty reuse-map channel — run_ani_reuse is on but no abinitio_reuse_assignments.csv was produced"
            def csv = file(p, glob: false)
            if (!csv.exists() || csv.size() == 0)
                error "PREDICT_REUSE: abinitio_reuse_assignments.csv not found at ${p} — the ANI front-end (ANI_REUSE) must run before prediction when --run_ani_reuse=true"
            FunannotateUtils.loadAbinitioReuseMap(p)
        }

    def classified = predict_input_ch
        .combine(reuseMapPerItem)
        .map { meta, gfa, reuseMap ->
            def assignment = reuseMap[meta.id]
            def is_rep     = assignment?.is_representative ?: false
            def eligible   = (assignment?.reuse_eligible ?: false) && !forceAll
            tuple(meta, gfa, is_rep, eligible)
        }

    def branched = classified.branch {
        representative:   it[2]
        eligible_sibling: it[3]
        independent:      true
    }

    // Representative out-ids (for backfill filtering), as a single-value Set.
    def repOutSet = branched.representative
        .map { meta, _gfa, _is_rep, _elig -> meta.id.toString() }
        .collect()
        .map { it as Set }
        .ifEmpty([] as Set)

    // ── Representatives + independents: shared fresh GeneMark + predict path ──
    // Both always pass shared_mod='' (reps define the store, independents were
    // never eligible to reuse it), so ONE GENEMARK_RUN call on the mixed channel
    // is correct — join-by-meta doesn't care which branch a row came from.
    def rep_and_indep = branched.representative.mix(branched.independent)
        .map { meta, gfa, _is_rep, _elig -> tuple(meta, gfa) }

    def fresh_todo = rep_and_indep.filter { meta, _gfa ->
        FunannotateUtils.gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) == null ||
            FunannotateUtils.staleRnaseq(meta.id as String, meta.species as String, params.target as String, launchDir.toString()) ||
            FunannotateUtils.staleGenome(meta.id as String, meta.asmid as String, params.source as String, params.target as String)
    }

    def fresh_with_gtf
    if (runGenemark) {
        def genemark_in = fresh_todo.map { meta, gfa ->
            def training_bam = FunannotateUtils.trainingTranscriptBamFor(meta.id as String, params.training_target as String)
            def mode = genemarkMode == 'AUTO' ? (training_bam ? 'ET' : 'ES') : genemarkMode
            tuple(meta, gfa, mode, training_bam, forceAll ? "true" : "false", '')
        }
        GENEMARK_RUN(genemark_in)
        fresh_with_gtf = fresh_todo.join(GENEMARK_RUN.out.gtf, by: 0)
            .map { meta, gfa, gtf -> tuple(meta, gfa, gtf) }
    } else {
        fresh_with_gtf = fresh_todo.map { meta, gfa -> tuple(meta, gfa, '') }
    }
    // Optional prodigal pass-through (--other_gff <gff>:<weight>), same as the
    // TRAIN_PREDICT non-reuse path. PRODIGAL_RUN runs on the same fresh_todo
    // that FUNANNOTATE_PREDICT consumes, so joins are 1:1. Contingency: see
    // prodigal_lineages — enforced inside PRODIGAL_RUN (non-matching lineages
    // emit an empty GFF that FUNANNOTATE_PREDICT drops).
    def runProdigal = (params.run_prodigal ?: false).toString().toBoolean()
    def fresh_final
    if (runProdigal) {
        PRODIGAL_RUN(fresh_todo)
        fresh_final = fresh_with_gtf.join(PRODIGAL_RUN.out.gff3, by: 0)
            .map { meta, gfa, gtf, other_gff -> tuple(meta, gfa, gtf, other_gff) }
    } else {
        fresh_final = fresh_with_gtf
            .map { meta, gfa, gtf -> tuple(meta, gfa, gtf, '') }
    }
    FUNANNOTATE_PREDICT(fresh_final)

    // ── Backfill every representative predicted this run into the shared store.
    // Batched (~100/batch); join by `out` with GENEMARK_RUN's fresh .mod (1:1,
    // reps never reuse). Reps-with-preexisting-current-gbk get filtered out of
    // fresh_todo above, so only reps actually (re)trained this run backfill.
    def freshBackfillInput
    if (runGenemark) {
        freshBackfillInput = FUNANNOTATE_PREDICT.out.metadata
            .map { meta -> tuple(meta.id.toString(), meta.species.toString()) }
            .combine(repOutSet)
            .filter { out, _sp, repSet -> repSet.contains(out) }
            .map { out, sp, _repSet -> tuple(out, sp) }
            .join(GENEMARK_RUN.out.mod.map { meta, mod -> tuple(meta.id.toString(), mod.toString()) })
            .map { out, sp, mod -> tuple(sp, out, mod) }
            .collate(100)
            .map { batch -> tuple(batch.hashCode(), batch) }
    } else {
        freshBackfillInput = FUNANNOTATE_PREDICT.out.metadata
            .map { meta -> tuple(meta.id.toString(), meta.species.toString()) }
            .combine(repOutSet)
            .filter { out, _sp, repSet -> repSet.contains(out) }
            .map { out, sp, _repSet -> tuple(sp, out, '') }
            .collate(100)
            .map { batch -> tuple(batch.hashCode(), batch) }
    }

    // ── Species with shared params available by the time siblings need them ──
    // pre-existing (on disk now) mixed with THIS run's freshly-backfilled ones
    // (real channel dependency on rep predict+backfill), collected to a Set.
    // combine() (not join(by:0, remainder:true)): a species key is 1-to-many
    // (many eligible siblings share it) — join silently satisfies only one left
    // row per key; combine-with-a-scalar-Set broadcasts to every sibling.
    def preexistingSpeciesCh = branched.eligible_sibling
        .map { meta, _gfa, _is_rep, _elig ->
            def sp   = meta.species.toString()
            def spTag = sp.replaceAll(/\s+/, '_')
            def json = file("${sharedRoot}/${spTag}/parameters.json", glob: false)
            tuple(sp, (json.exists() && json.size() > 0) ? 'yes' : 'no')
        }
        .filter { sp, has -> has == 'yes' }
        .map { sp, _has -> sp }

    if (sharedRoot) {
        BACKFILL_ABINITIO_PARAMS(freshBackfillInput)
    } else {
        // run_ani_reuse without a shared store is a config error caught in
        // funannotate.nf; keep the channel graph well-formed regardless.
        freshBackfillInput.collect().ifEmpty([]).map { true }
    }

    def freshSpeciesCh = sharedRoot
        ? BACKFILL_ABINITIO_PARAMS.out.done
            .flatMap { items, _marker -> items.collect { sp, _rep_out, _mod -> sp.toString() } }
        : Channel.empty()

    def availableSpeciesSet = preexistingSpeciesCh
        .mix(freshSpeciesCh)
        .unique()
        .collect()
        .map { it as Set }
        .ifEmpty([] as Set)

    // ── Eligible siblings: gate on availability ──────────────────────────────
    def eligibleSiblingKeyed = branched.eligible_sibling
        .map { meta, gfa, _is_rep, _elig -> tuple(meta.species.toString(), tuple(meta, gfa)) }

    def gated = eligibleSiblingKeyed.combine(availableSpeciesSet)

    def readyRows = gated
        .filter { sp, _row, availSet -> availSet.contains(sp) }
        .map { sp, row, _availSet -> tuple(row[0], row[1], sp) }

    def blockedRows = gated
        .filter { sp, _row, availSet -> !availSet.contains(sp) }
        .map { _sp, row, _availSet -> row }

    def sibling_todo
    if (allowFallback) {
        def fallbackRows = blockedRows.map { meta, gfa ->
            log.warn "predict: ${meta.id} (${meta.species}) — shared ab-initio parameters not available; training independently (--allow_independent_fallback)"
            tuple(meta, gfa, meta.species.toString())
        }
        sibling_todo = readyRows.mix(fallbackRows)
    } else {
        blockedRows
            .map { meta, _gfa ->
                "${meta.species}\t${meta.id}\t${meta.asmid}\tshared_ab_initio_params_not_available"
            }
            .collectFile(name: 'predict_blocked_awaiting_representative.tsv',
                         storeDir: params.target, newLine: true, sort: true,
                         seed: 'species\tout\tasmid\treason')
        sibling_todo = readyRows
    }

    // A sibling's predict is also stale when the shared parameters.json it was
    // built on has since been refreshed (rep re-annotated / assignment changed).
    def sibling_predict_todo = sibling_todo
        .filter { meta, _gfa, sp ->
            def sharedJson = FunannotateUtils.sharedParamsJsonFor(sp as String, sharedRoot as String)
            FunannotateUtils.gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) == null ||
                FunannotateUtils.staleRnaseq(meta.id as String, meta.species as String, params.target as String, launchDir.toString()) ||
                FunannotateUtils.staleGenome(meta.id as String, meta.asmid as String, params.source as String, params.target as String) ||
                FunannotateUtils.staleSharedParams(meta.id as String, sharedJson ? sharedJson.toString() : '', params.target as String)
        }

    // Eligible siblings: GENEMARK_RUN_SIB reads post-gating (readyRows), so
    // shared_mod resolves only once availability is known. A forced-independent
    // sibling's own fresh .mod is deliberately never wired into backfill —
    // only representatives define the shared per-species store.
    if (runGenemark) {
        def sibling_genemark_in = sibling_predict_todo.map { meta, gfa, sp ->
            def training_bam = FunannotateUtils.trainingTranscriptBamFor(meta.id as String, params.training_target as String)
            def mode = genemarkMode == 'AUTO' ? (training_bam ? 'ET' : 'ES') : genemarkMode
            def sharedMod = FunannotateUtils.sharedGenemarkModFor(sp as String, sharedRoot as String)
            tuple(meta, gfa, mode, training_bam, forceAll ? "true" : "false", sharedMod ? sharedMod.toString() : '')
        }
        GENEMARK_RUN_SIB(sibling_genemark_in)

        sibling_predict_todo = sibling_predict_todo
            .join(GENEMARK_RUN_SIB.out.gtf, by: 0)
            .map { meta, gfa, sp, gtf -> tuple(meta, gfa, gtf) }
    } else {
        sibling_predict_todo = sibling_predict_todo.map { meta, gfa, sp -> tuple(meta, gfa, '') }
    }

    // Optional prodigal pass-through for reusing siblings, mirroring the fresh
    // path above (PRODIGAL_RUN_SIB on the same gated sibling_predict_todo).
    // Contingency: enforced inside the module (non-matching BUSCO lineages
    // emit an empty GFF -> no --other_gff).
    def runProdigalSib = (params.run_prodigal ?: false).toString().toBoolean()
    def sibling_final
    if (runProdigalSib) {
        PRODIGAL_RUN_SIB(sibling_predict_todo.map { meta, gfa, sp -> tuple(meta, gfa) })
        sibling_final = sibling_predict_todo.join(PRODIGAL_RUN_SIB.out.gff3, by: 0)
            .map { meta, gfa, gtf, other_gff -> tuple(meta, gfa, gtf, other_gff) }
    } else {
        sibling_final = sibling_predict_todo
            .map { meta, gfa, gtf -> tuple(meta, gfa, gtf, '') }
    }

    FUNANNOTATE_PREDICT_SIB(sibling_final)

    emit:
    metadata = FUNANNOTATE_PREDICT.out.metadata.mix(FUNANNOTATE_PREDICT_SIB.out.metadata)
    done     = FUNANNOTATE_PREDICT.out.done.mix(FUNANNOTATE_PREDICT_SIB.out.done)
}
