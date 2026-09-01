/*
 * TRAIN_PREDICT — RNA-seq training + gene prediction for one batch of assemblies.
 *
 * Takes:
 *   ch_genomes   — channel: tuple(val(meta), val(genome_fa))  (from MASK_GENOME)
 *   ch_reads     — channel: tuple(val(species_tag), path(r1), path(r2), path(se))
 *                  May be empty when run_sra_fetch=false or stop_after_sra_fetch=true.
 *   reuse_map_csv — path to abinitio_reuse_assignments.csv (ANI_REUSE's output).
 *                  '' (default) when run_ani_reuse=false; drives PREDICT_REUSE gating.
 *
 * Steps:
 *   1. When reads are available:
 *      a. RNASEQ_PREPARE: run once per species on the representative assembly (storeDir).
 *      b. FUNANNOTATE_TRAIN: run PASA alignment for every assembly using the shared Trinity.
 *         Species with no reads bypass training entirely.
 *   2. Prediction routing:
 *      a. run_ani_reuse=true  → PREDICT_REUSE (BFD FUNANNOTATE_PREDICTION port): reps
 *         train fresh + BACKFILL_ABINITIO_PARAMS into the shared store; eligible siblings
 *         fast-reuse the species' shared GeneMark .mod via GENEMARK_RUN_SIB (--predict_with)
 *         once it's available; independents train ungated; blocked siblings reported to
 *         predict_blocked_awaiting_representative.tsv (unless allow_independent_fallback).
 *      b. run_ani_reuse=false → original inline path below (everyone independent).
 *
 * Emits:
 *   metadata — val(meta) for each assembly that completed predict
 *   done     — path to *.predict.done marker file
 */

include { RNASEQ_PREPARE    } from './../../modules/local/rnaseq_prepare'
include { FUNANNOTATE_TRAIN } from './../../modules/local/funannotate_train'
include { FUNANNOTATE_PREDICT } from './../../modules/local/funannotate_predict'
include { GENEMARK_RUN      } from './../../modules/local/genemark_run'
include { PRODIGAL_RUN      } from './../../modules/local/prodigal_run'
include { BACKFILL_ABINITIO_PARAMS } from './../../modules/local/backfill_abinitio_params'
include { PREDICT_REUSE     } from './predict_reuse'

workflow TRAIN_PREDICT {

    take:
    ch_genomes  // tuple(val(meta), val(genome_fa))
    ch_reads    // tuple(val(species_tag), path(r1), path(r2), path(se))  — may be empty
    reuse_map_csv  // path — ANI_REUSE's abinitio_reuse_assignments.csv ('' = disabled)

    main:
    def predict_input_ch

    if (params.run_sra_fetch.toBoolean() && !params.stop_after_sra_fetch.toBoolean()) {
        // ── Build per-assembly channel keyed by species_tag with SRA reads joined ──
        // reads_ch is a 4-tuple: (species_tag, r1, r2, se)
        def assembly_with_reads = ch_genomes
            .map { meta, genome_fa ->
                def species_tag = meta.species.replaceAll(/\s+/, '_')
                tuple(species_tag, meta, genome_fa)
            }
            .combine(ch_reads, by: 0)
        // assembly_with_reads: (species_tag, meta, genome_fa, r1, r2, se)

        // ── RNASEQ_PREPARE: once per species on the representative (first) assembly ──
        // storeDir-cached; subsequent runs reuse the Trinity-GG FASTA.
        def repr_ch = assembly_with_reads
            .groupTuple(by: 0)
            .map { species_tag, metas, genomes, r1s, r2s, ses ->
                tuple(species_tag, metas[0], genomes[0], r1s[0], r2s[0], ses[0])
            }

        def repr_branched = repr_ch.branch {
            has_reads: it[3].size() > 0 || it[5].size() > 0   // r1=[3] or se=[5]
            no_reads:  true
        }

        RNASEQ_PREPARE(repr_branched.has_reads)

        // For species with no RNA-seq, write an empty Trinity FASTA in the driver process
        // (no SLURM job) and emit it directly as a shared channel item.
        def empty_shared_ch = repr_branched.no_reads
            .map { species_tag, _meta, _gfa, _r1, _r2, _se ->
                def empty_fa = file("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
                if (!empty_fa.exists()) {
                    empty_fa.parent.mkdirs()
                    empty_fa.text = ''
                }
                tuple(species_tag, empty_fa)
            }

        def shared_ch = RNASEQ_PREPARE.out.shared.mix(empty_shared_ch)

        // ── Join shared Trinity back to every assembly for FUNANNOTATE_TRAIN ──
        def train_input = assembly_with_reads
            .combine(shared_ch, by: 0)
            .map { species_tag, meta, genome_fa, r1, r2, se, trinity_fa ->
                tuple(meta, genome_fa, r1, r2, se, trinity_fa)
            }
        // train_input: (meta, genome_fa, r1, r2, se, trinity_fa)

        // Assemblies with no RNA-seq bypass FUNANNOTATE_TRAIN entirely.
        def branched = train_input.branch {
            has_rnaseq: it[2].size() > 0 || it[4].size() > 0 || it[5].size() > 0
            no_rnaseq:  true
        }
        def predict_no_rnaseq = branched.no_rnaseq
            .map { meta, genome_fa, _r1, _r2, _se, _tf -> tuple(meta, genome_fa) }

        // Skip TRAIN when pasa.gff3 already exists and is not stale relative to reads.
        def train_todo = branched.has_rnaseq.filter { meta, _gfa, _r1, _r2, _se, _tf ->
            def gff3 = file("${params.training_target}/${meta.id}/training/funannotate_train.pasa.gff3")
            !gff3.exists() || gff3.size() == 0 || FunannotateUtils.staleRnaseq(meta.id as String, meta.species as String, params.target as String, launchDir.toString())
        }
        def train_done = branched.has_rnaseq
            .filter { meta, _gfa, _r1, _r2, _se, _tf ->
                def gff3 = file("${params.training_target}/${meta.id}/training/funannotate_train.pasa.gff3")
                gff3.exists() && gff3.size() > 0 && !FunannotateUtils.staleRnaseq(meta.id as String, meta.species as String, params.target as String, launchDir.toString())
            }
            .map { meta, genome_fa, _r1, _r2, _se, _tf -> tuple(meta, genome_fa) }

        FUNANNOTATE_TRAIN(train_todo)
        predict_input_ch = FUNANNOTATE_TRAIN.out.mix(train_done).mix(predict_no_rnaseq)
    } else {
        // No RNA-seq: pass genomes straight to predict.
        predict_input_ch = ch_genomes
    }

    // ── FUNANNOTATE_PREDICT / prediction routing ─────────────────────────────
    def metadata_out
    def done_out

    if (params.run_ani_reuse.toBoolean()) {
        // Representative-gated prediction (BFD FUNANNOTATE_PREDICTION port).
        // PREDICT_REUSE does its own not-yet-predicted/stale filtering internally
        // (fresh_todo / sibling_predict_todo), mirroring the filter below, so it
        // takes the raw post-train channel and gates + reuses inside.
        PREDICT_REUSE(predict_input_ch, reuse_map_csv)
        metadata_out = PREDICT_REUSE.out.metadata
        done_out     = PREDICT_REUSE.out.done
    } else {
    def predict_ch = predict_input_ch
        .filter { meta, _gfa ->
            FunannotateUtils.gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) == null ||
            FunannotateUtils.staleRnaseq(meta.id as String, meta.species as String, params.target as String, launchDir.toString()) ||
            FunannotateUtils.staleGenome(meta.id as String, meta.asmid as String, params.source as String, params.target as String)
        }

    // ── GENEMARK_RUN (standalone, host-side) ────────────────────────────────
    // Runs only for assemblies actually being predicted, so it emits exactly one
    // GTF per predict task. mode resolves params.genemark_mode for each assembly:
    // 'auto' -> ET when its FUNANNOTATE_TRAIN BAM exists, else ES. shared_mod
    // fast-reuses the species' shared GeneMark model when present (unless
    // force_independent); with no shared-root configured it's null -> fresh training.
    def genemark_input = predict_ch.map { meta, genome_fa ->
        def training_bam = FunannotateUtils.trainingTranscriptBamFor(meta.id as String, params.training_target as String)
        def mode = params.genemark_mode == 'auto' ? (training_bam ? 'ET' : 'ES') : params.genemark_mode
        def shared_root = params.gene_prediction_shared_abinitio
        def shared_mod  = shared_root ? FunannotateUtils.sharedGenemarkModFor(meta.species as String, shared_root as String) : null
        tuple(meta, genome_fa, mode, training_bam, params.force_independent, shared_mod ? shared_mod.toString() : '')
    }
    GENEMARK_RUN(genemark_input)

    // Join the stand-alone GTF back to each predict task. FUNANNOTATE_PREDICT
    // passes it as --genemark_gtf (empty string -> --auto-skip-genemark).
    // PRODIGAL_RUN (optional, params.run_prodigal) runs on the same predict_ch
    // and its GFF is passed as --other_gff <gff>:<weight> (empty when disabled).
    // Contingency: params.prodigal_lineages (list of BUSCO lineage names, e.g.
    // ['microsporidia_odb10']) is enforced inside PRODIGAL_RUN — genomes whose
    // lineage is NOT in the list emit an empty GFF, which FUNANNOTATE_PREDICT
    // drops (other_gff_ok requires size>0), so they get no --other_gff.
    def runProdigal = (params.run_prodigal ?: false).toString().toBoolean()
    def gtf_ch = predict_ch.join(GENEMARK_RUN.out.gtf, by: 0)
    def predict_final
    if (runProdigal) {
        PRODIGAL_RUN(predict_ch)
        predict_final = gtf_ch.join(PRODIGAL_RUN.out.gff3, by: 0)
            .map { meta, genome_fa, genemark_gtf, other_gff -> tuple(meta, genome_fa, genemark_gtf, other_gff) }
    } else {
        predict_final = gtf_ch
            .map { meta, genome_fa, genemark_gtf -> tuple(meta, genome_fa, genemark_gtf, '') }
    }

    FUNANNOTATE_PREDICT(predict_final)

    // ── BACKFILL_ABINITIO_PARAMS (fresh .mod -> shared per-species store) ─────
    // Only when a shared ab-initio store is configured. Every species freshly
    // trained this run (GENEMARK_RUN emitted a new .mod because shared_mod was
    // absent or force_independent) and completed predict is backfilled into the
    // store; the inner join with GENEMARK_RUN.out.mod drops fast-reuse rows
    // (which produce no new .mod). Batched via .collate(100) so one job loops
    // over up to ~100 species; each line is idempotent (content-hash
    // short-circuit in bin/backfill_abinitio_params.py). Representative-vs-
    // sibling filtering is future work: the compare_ANI leg that computes
    // is_representative isn't ported yet, so for now every fresh-trained species
    // defines the shared store for itself.
    if (params.gene_prediction_shared_abinitio) {
        def backfill_input = FUNANNOTATE_PREDICT.out.metadata
            .map { meta -> tuple(meta.id.toString(), meta.species.toString()) }
            .join(GENEMARK_RUN.out.mod.map { meta, mod -> tuple(meta.id.toString(), mod.toString()) })
            .map { out, sp, mod -> tuple(sp, out, mod) }
            .collate(100)
            .map { batch -> tuple(batch.hashCode(), batch) }
        BACKFILL_ABINITIO_PARAMS(backfill_input)
    }

    metadata_out = FUNANNOTATE_PREDICT.out.metadata
    done_out     = FUNANNOTATE_PREDICT.out.done
    } // end run_ani_reuse=false inline branch

    emit:
    metadata = metadata_out
    done     = done_out
}
