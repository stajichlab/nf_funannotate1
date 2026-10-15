// GENEMARK_RUN: drop-in replacement for the old GENEMARK_RUN process (same take, same emits gtf and mod), now a
// subworkflow. Callers keep calling GENEMARK_RUN(...) and GENEMARK_RUN.out.{gtf,mod}.
//
//   GENEMARK_TRAIN --+-- monolithic / skip route ---------------------------------> gtf, mod
//                    '-- chunked route (params.genemark_chunked, >= N contigs, fresh ES, container mode)
//                          chunks > GENEMARK_PREDICT_CHUNK (one per chunk, --predict_with the model)
//                                    > groupTuple(groupKey(id, n_chunks)) > GENEMARK_MERGE ----> gtf, mod
//
// With params.genemark_chunked = false (default) this behaves exactly like the old single process.
// See modules/local/genemark_train.nf for what is and is not chunked.

include { GENEMARK_TRAIN }         from '../../modules/local/genemark_train'
include { GENEMARK_PREDICT_CHUNK } from '../../modules/local/genemark_predict_chunk'
include { GENEMARK_MERGE }         from '../../modules/local/genemark_merge'

workflow GENEMARK_RUN {
    take:
    rows   // tuple(meta, genome_fa, mode, training_bam, force_independent, shared_mod, shared_fp)

    main:
    GENEMARK_TRAIN(rows)

    // per-genome facts needed later, keyed by meta.id
    meta_ch = rows.map { meta, genome_fa, _mode, _bam, _fi, _sm, _fp -> tuple(meta.id, meta, genome_fa) }
    model   = GENEMARK_TRAIN.out.train_mod.map { meta, mod -> tuple(meta.id, mod) }

    // scatter: one row per chunk file, in chunk order
    chunk_rows = GENEMARK_TRAIN.out.chunks
        .flatMap { meta, files ->
            // sort on the integer in chunk_<NNNN>.fa.gz, not on the name, so ordering never depends on zero padding
            def fl = (files instanceof List ? files : [files]).sort { (it.name =~ /chunk_(\d+)\.fa\.gz/)[0][1] as int }
            fl.withIndex().collect { f, i -> tuple(meta.id, String.format('%04d', i), fl.size(), f) }
        }
        .combine(model, by: 0)
        .combine(meta_ch, by: 0)
        .map { _id, idx, n, chunk, mod, meta, _genome -> tuple(meta, idx, n, chunk, mod) }

    GENEMARK_PREDICT_CHUNK(chunk_rows)

    // gather: wait for all n chunks of a genome (groupKey carries the expected count)
    gathered = GENEMARK_PREDICT_CHUNK.out.gtf
        .map { meta, _idx, n, gtf -> tuple(groupKey(meta.id, n as int), gtf) }
        .groupTuple()
        .map { key, gtfs -> tuple(key.getGroupTarget(), key.getGroupSize(), gtfs) }

    // failOnDuplicate: two rows with the same id would cross-pair chunks and models. (Not failOnMismatch: one
    // permanently failed chunk must only drop its own genome, not stop the run.)
    merge_in = gathered
        .join(model, failOnDuplicate: true)
        .join(meta_ch, failOnDuplicate: true)
        .map { _id, n, gtfs, mod, meta, genome -> tuple(meta, n, genome, gtfs, mod) }

    GENEMARK_MERGE(merge_in)

    emit:
    gtf = GENEMARK_TRAIN.out.gtf.mix(GENEMARK_MERGE.out.gtf)
    mod = GENEMARK_TRAIN.out.mod.mix(GENEMARK_MERGE.out.mod)
}
