/*
 * TRAIN_PREDICT — RNA-seq training + gene prediction for one batch of assemblies.
 *
 * Takes:
 *   ch_genomes  — channel: tuple(val(meta), val(genome_fa))  (from MASK_GENOME)
 *   ch_reads    — channel: tuple(val(species_tag), path(r1), path(r2), path(se))
 *                 May be empty when run_sra_fetch=false or stop_after_sra_fetch=true.
 *
 * Steps:
 *   1. When reads are available:
 *      a. RNASEQ_PREPARE: run once per species on the representative assembly (storeDir).
 *      b. FUNANNOTATE_TRAIN: run PASA alignment for every assembly using the shared Trinity.
 *         Species with no reads bypass training entirely.
 *   2. FUNANNOTATE_PREDICT: run for all assemblies not yet predicted (or stale).
 *
 * Emits:
 *   metadata — val(meta) for each assembly that completed predict
 *   done     — path to *.predict.done marker file
 */

include { RNASEQ_PREPARE    } from './../../modules/local/rnaseq_prepare'
include { FUNANNOTATE_TRAIN } from './../../modules/local/funannotate_train'
include { FUNANNOTATE_PREDICT } from './../../modules/local/funannotate_predict'

// GBK may be stored compressed to save space; return the first existing non-empty file.
def _gbkResult(String dir, String id) {
    def plain = file("${dir}/${id}.gbk")
    if (plain.exists() && plain.size() > 0) return plain
    def gz = file("${dir}/${id}.gbk.gz")
    if (gz.exists() && gz.size() > 0) return gz
    return null
}

// True when rnaseq reads or shared Trinity are newer than the existing prediction GBK.
def _staleRnaseq(String id, String species) {
    def species_tag = species.replaceAll(/\s+/, '_')
    def gbk = _gbkResult("${params.target}/${id}/predict_results", id)
    if (gbk == null) return false
    def gbkMod = gbk.lastModified()
    def r1      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_R1.fastq.gz")
    def se      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_SE.fastq.gz")
    def trinity = file("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
    return (r1.exists() && r1.size() > 0 && r1.lastModified() > gbkMod) ||
           (se.exists() && se.size() > 0 && se.lastModified() > gbkMod) ||
           (trinity.exists() && trinity.size() > 0 && trinity.lastModified() > gbkMod)
}

workflow TRAIN_PREDICT {

    take:
    ch_genomes  // tuple(val(meta), val(genome_fa))
    ch_reads    // tuple(val(species_tag), path(r1), path(r2), path(se))  — may be empty

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
            !gff3.exists() || gff3.size() == 0 || _staleRnaseq(meta.id as String, meta.species as String)
        }
        def train_done = branched.has_rnaseq
            .filter { meta, _gfa, _r1, _r2, _se, _tf ->
                def gff3 = file("${params.training_target}/${meta.id}/training/funannotate_train.pasa.gff3")
                gff3.exists() && gff3.size() > 0 && !_staleRnaseq(meta.id as String, meta.species as String)
            }
            .map { meta, genome_fa, _r1, _r2, _se, _tf -> tuple(meta, genome_fa) }

        FUNANNOTATE_TRAIN(train_todo)
        predict_input_ch = FUNANNOTATE_TRAIN.out.mix(train_done).mix(predict_no_rnaseq)
    } else {
        // No RNA-seq: pass genomes straight to predict.
        predict_input_ch = ch_genomes
    }

    // ── FUNANNOTATE_PREDICT ─────────────────────────────────────────────────
    def predict_ch = predict_input_ch
        .filter { meta, _gfa ->
            _gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) == null ||
            _staleRnaseq(meta.id as String, meta.species as String)
        }
    FUNANNOTATE_PREDICT(predict_ch)

    emit:
    metadata = FUNANNOTATE_PREDICT.out.metadata
    done     = FUNANNOTATE_PREDICT.out.done
}
