/*
 * ANI_REUSE — in-run ANI front-end for representative-strain ab-initio reuse.
 *
 * Spliced into funannotate.nf right after CLEAN_GENOMES (consumes the cleaned,
 * pre-mask genomes — same inputs BFD's compare_ani/busco_genome pipelines read)
 * when --run_ani_reuse=true. Everything runs inside the same funannotate
 * invocation as prediction; no separate compare_ani or busco_genome entrypoint
 * is needed.
 *
 * Steps (all storeDir/publishDir-cached, so -resume reuses prior results):
 *   - BUSCO_GENOME   per assembly (genome mode; written to
 *                    params.genome_stats_outdir/BUSCO_genome as BFD does).
 *   - ASM_STATS      per-assembly length/N50/contig table (params.tables_dir).
 *   - SKANI_COMPARE  symmetric all-vs-all skani ANI within each species group
 *                    → <species>.full.ani.tsv.
 *   - CONCAT_ANI_TSVS merge → all_pairs_merged.tsv + asmid manifest.
 *   - PICK_REPRESENTATIVE_STRAIN → abinitio_reuse_assignments.csv (+
 *                    repr_assignments.tsv), the map TRAIN_PREDICT/PREDICT_REUSE
 *                    consume to gate sibling reuse of the shared per-species
 *                    ab-initio store.
 *
 * Port of BFD/Fungi_BFD/nextflow workfows/compare_ANI.nf's --run_ani_reuse leg,
 * compressed into a front-end so representative picking is deterministic IN the
 * same run that gates on it.
 */

include { BUSCO_GENOME }        from './../../modules/local/busco_genome'
include { ASM_STATS }           from './../../modules/local/asm_stats'
include { SKANI_COMPARE }       from './../../modules/local/skani_compare'
include { CONCAT_ANI_TSVS }     from './../../modules/local/concat_ani_tsvs'
include { PICK_REPRESENTATIVE_STRAIN } from './../../modules/local/pick_representative_strain'

workflow ANI_REUSE {

    take:
    ch_genomes    // tuple(val(meta), path(genome_fa)) — cleaned, pre-mask genomes
    samples_csv   // path to samples.csv

    main:
    // ── Per-assembly BUSCO completeness (genome mode) — represents the
    // completeness axis of representative picking (mirrors BFD, which reads the
    // same <asmid>.BUSCO_summary.<lineage>.txt files from this same storeDir).
    BUSCO_GENOME(ch_genomes)
    def busco_summaries = BUSCO_GENOME.out.summary

    // ── Assembly stats (total bp / N50 / contig count) — the N50 tiebreak axis.
    // params.genome_dir points at launchDir/input_clean_genomes (the cleaned
    // genomes), matching what BUSCO_GENOME/SKANI_COMPARE consume above.
    ASM_STATS(samples_csv, file(params.genome_dir))
    def asm_stats = ASM_STATS.out.stats

    // ── Group by species (sanitized tag, same convention as species_tag used
    // everywhere else in this repo) and run symmetric all-vs-all skani ANI.
    def grouped = ch_genomes
        .map { meta, gfa ->
            def species_tag = meta.species.replaceAll(/\s+/, '_')
            tuple(species_tag, meta, gfa)
        }
        .groupTuple(by: 0)
        .map { species_tag, metas, gfxs ->
            tuple(species_tag, metas.size(), gfxs.collect { g -> file(g, glob: false) })
        }
    SKANI_COMPARE(grouped)
    def ani_tsv_ch = SKANI_COMPARE.out

    // ── Manifest of every group's ANI TSV (feed collectFile the channel of
    // paths directly, one per group — never a single collected List: collectFile
    // would fold a real Path into the collector by filename instead of writing
    // its path as a line; see BFD .living/learnings.md CONCAT reuse-0 bug).
    def ani_manifest = ani_tsv_ch
        .map { _group, tsv -> tsv.toString() }
        .collectFile(name: 'tsv_manifest.txt', newLine: true, sort: true)

    // ── Companion asmid manifest (deduped, sorted) for CONCAT's staleness file.
    def asmid_manifest = ch_genomes
        .map { meta, _gfa -> meta.asmid.toString() }
        .unique()
        .collectFile(name: 'asmid_manifest.txt', newLine: true, sort: true)

    CONCAT_ANI_TSVS(ani_manifest, asmid_manifest)

    // ── predict_input TSV: one tab-separated row per strain, same columns as
    // BFD's WRITE_ANI_PREDICT_INPUT (out, asmid, species, strain, locustag,
    // busco, header_length, transl_table, genome). out == meta.id so the CSV
    // keys (by out) join straight onto TRAIN_PREDICT's (meta, genome) rows.
    // collectFile without sort so the header seed stays the FIRST line that
    // csv.DictReader(fh, delimiter='\t') reads as the column list.
    def predict_input_tsv = ch_genomes
        .map { meta, gfa ->
            [ meta.id, meta.asmid, meta.species, meta.strain, meta.locustag,
              meta.busco, params.header_length, meta.transl_table, gfa.toString() ].join('\t')
        }
        .collectFile(name: 'predict_input_for_ani.tsv', newLine: true,
                     seed: 'out\tasmid\tspecies\tstrain\tlocustag\tbusco\thlen\ttable\tgenome_fa')

    PICK_REPRESENTATIVE_STRAIN(
        CONCAT_ANI_TSVS.out.out.ifEmpty(file('/dev/null')),
        predict_input_tsv,
        samples_csv,
        busco_summaries,
        asm_stats
    )

    emit:
    reuse_map    = PICK_REPRESENTATIVE_STRAIN.out.outCSV
    assignments  = PICK_REPRESENTATIVE_STRAIN.out.outAssignments
}
