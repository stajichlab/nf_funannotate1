# Chunked GeneMark (train once, scatter predict, gather)

Applies to `genemark_sidecar.nf` only. `GENEMARK_RUN` (`subworkflows/local/genemark_run.nf`) keeps the old take and emits (`gtf`, `mod`).
`train_predict.nf` and `predict_reuse.nf` still include the old process in `modules/local/genemark_run.nf`.

## Route
- `genemark_chunked = false` (default): every genome takes the old single-task route (`GENEMARK_TRAIN`, the old process body).
- `genemark_chunked = true`: a genome with at least `genemark_chunk_min_contigs` (25,000) contigs that would train FRESH in ES
  mode, in container mode, stops after `gmes_petap.pl --ES --training_only` and writes `chunks/chunk_NNNN.fa.gz` (about
  `genemark_chunk_contigs` = 20,000 contigs each, cut by contig count, genome order). `GENEMARK_PREDICT_CHUNK` runs
  `gmes_petap.pl --predict_with <model> --max_intron --soft_mask 2000 --fungus` per chunk. `GENEMARK_MERGE` merges in chunk
  order, renumbers IDs, validates against the whole genome, and emits the GTF and the model (the model is emitted only after the merge).
- Not chunked: fast-reuse rows (this repo's reuse call omits `--max_intron` and `--soft_mask`, so a chunked reuse would not match it),
  fresh ET rows with a training BAM (ET equivalence is untested), and host (non-container) mode.

## Queues (UCR HPCC site configs)
Chunks: `short` 2 h first, then `preempt` 6 h, then the default queue 12 h; `maxForks` 32 (8 CPUs each). A chunk or merge that fails
every attempt is ignored: that genome gets no GTF. After a run:
`bin/genemark_chunked_missing.py --chunked-tsv <target>/genemark_chunked.tsv --store <published dir>` lists chunked genomes with no result.

## Checks
- The process shell has no `set -e`. Chunk and merge scripts set strict mode; `GENEMARK_TRAIN` checks every exit 0 for the outputs of
  exactly one route, and fails when `asm_preflight_stats.py` produces no verdict (an empty verdict used to skip the genome silently).
- A chunk passes only with one result file per contig piece, no gmhmme3 warning, a GTF, and contig names that belong to the chunk.

## Equivalence
ES training is not deterministic, so compare with the same published model: run `gmes_petap.pl --predict_with <published .mod>` on
the whole genome and compare with the chunked GTF (`tests/genemark_chunked/README.md`). On a 1,759-contig genome (4 chunks) the
GTFs are identical (67,486 records). The BFD pipeline's equivalent change was also checked on a 241,057-contig genome.

## Not tested here
NRP/k8s pods (no nested apptainer, no `$SCRATCH`, 4 CPU / 16 GB), genetic-code tables 6/26, ET, and thresholds below 25,000 contigs.
Lower thresholds and non-fungal genomes with many contigs are a later performance test.
