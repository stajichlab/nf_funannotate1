/*
 * FETCH_RNASEQ — NCBI SRA query + download + normalization
 *
 * Steps:
 *   1. Query SRA per species (SRA_QUERY_BATCH or reuse cached CSVs via skip_sra_query).
 *   2. Merge all per-species CSVs into {stem}.rnaseq_sra.csv (COLLECT_SRA_QUERY).
 *   3. Route each species to SRA_FETCH (PE), SRA_FETCH_SE (SE_trinity / SINGLE),
 *      or WRITE_EMPTY_READS (no data). Skipped when stop_after_sra_query=true.
 *
 * Emits:
 *   reads        — channel: tuple(val(species_tag), path(r1), path(r2), path(se))
 *                  Empty channel when stop_after_sra_query=true.
 *   sra_manifest — path to {stem}.rnaseq_sra.csv (always emitted after step 2).
 */

include { SRA_QUERY_BATCH   } from './../../modules/local/sra_query_batch'
include { COLLECT_SRA_QUERY } from './../../modules/local/collect_sra_query'
include { WRITE_EMPTY_READS } from './../../modules/local/write_empty_reads'
include { SRA_FETCH         } from './../../modules/local/sra_fetch'
include { SRA_FETCH_SE      } from './../../modules/local/sra_fetch_se'

workflow FETCH_RNASEQ {

    take:
    ch_species  // channel: tuple(val(species_tag), val(taxonid))

    main:
    // ── Step 1: SRA query ────────────────────────────────────────────────────
    def sra_query_results
    if (params.skip_sra_query.toBoolean()) {
        // Reuse pre-existing per-species CSVs; skip all NCBI network calls.
        sra_query_results = ch_species
            .map { species_tag, _taxonid ->
                def csv = file("${launchDir}/rnaseq_reads/sra_query/${species_tag}.sra_query.csv")
                if (!csv.exists()) {
                    log.warn "skip_sra_query: no cached CSV for ${species_tag} — skipping this species"
                    return null
                }
                tuple(species_tag, csv)
            }
            .filter { it != null }
    } else {
        def sra_batched = ch_species
            .collate(params.sra_query_batch_size)
            .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }) }
        SRA_QUERY_BATCH(sra_batched)
        sra_query_results = SRA_QUERY_BATCH.out.query_results
            .flatten()
            .map { csv -> tuple(csv.baseName.replaceAll(/\.sra_query$/, ''), csv) }
    }

    // ── Step 2: collect manifest ─────────────────────────────────────────────
    def stem = file(params.samples).baseName
    COLLECT_SRA_QUERY(
        sra_query_results.map { _stag, csv -> csv }.collect(),
        stem
    )

    // ── Step 3: fetch reads (skipped when stop_after_sra_query=true) ─────────
    def ch_reads
    if (params.stop_after_sra_query.toBoolean()) {
        log.info "[FETCH_RNASEQ] stop_after_sra_query=true; skipping SRA download"
        ch_reads = Channel.empty()
    } else {
        // Build the blacklist map once (O(1) per-accession lookup in the closures).
        def blPath = file("${launchDir}/rnaseq_blacklist.csv")
        def blMap = blPath.exists()
            ? blPath.readLines().drop(1)
                  .findAll { it.trim() && !it.startsWith('#') }
                  .collectEntries { line ->
                      def cols = line.split(',')
                      cols.size() >= 4 ? [(cols[0].trim()): cols[3].trim()] : [:]
                  }
            : [:]

        // csvHasPE: at least one PAIRED accession not blocked or overridden to SE.
        def csvHasPE = { csv ->
            csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                def cols = line.split(',')
                if (cols.size() < 3) return false
                def layout = cols.size() > 5 ? cols[5].trim() : 'PAIRED'
                def action = blMap.get(cols[2].trim(), '')
                layout == 'PAIRED' && action != 'skip' && action != 'SE_trinity'
            }
        }

        // csvHasSEtrinity: at least one PAIRED accession overridden to SE via blacklist.
        def csvHasSEtrinity = { csv ->
            csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                def cols = line.split(',')
                cols.size() >= 3 && blMap.get(cols[2].trim(), '') == 'SE_trinity'
            }
        }

        // csvHasSingleLayout: at least one genuine SINGLE-layout accession (when enable_single_end=true).
        def csvHasSingleLayout = { csv ->
            csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                def cols = line.split(',')
                cols.size() > 5 && cols[5].trim() == 'SINGLE' && blMap.get(cols[2].trim(), '') != 'skip'
            }
        }

        // Three-way branch:
        //   has_pe  → SRA_FETCH   (PE wins; SE_trinity entries are handled inside SRA_FETCH)
        //   has_se  → SRA_FETCH_SE (SE_trinity always; SINGLE layout only if enable_single_end)
        //   no_data → WRITE_EMPTY_READS
        def branched_sra = sra_query_results
            .branch {
                has_pe: csvHasPE.call(it[1])
                has_se: csvHasSEtrinity.call(it[1]) ||
                        (params.enable_single_end.toBoolean() && csvHasSingleLayout.call(it[1]))
                no_data: true
            }

        SRA_FETCH(branched_sra.has_pe)
        SRA_FETCH_SE(branched_sra.has_se)
        WRITE_EMPTY_READS(branched_sra.no_data.map { stag, _csv -> stag })

        // Merge candidate files for human review (not consumed by downstream processes).
        SRA_FETCH.out.se_candidates
            .collectFile(name: 'rnaseq_se_candidates.csv', storeDir: launchDir, newLine: false)
        SRA_FETCH.out.blacklist_candidates
            .collectFile(name: 'rnaseq_blacklist_candidates.csv', storeDir: launchDir, newLine: false)

        ch_reads = SRA_FETCH.out.reads
            .mix(SRA_FETCH_SE.out.reads)
            .mix(WRITE_EMPTY_READS.out.reads)
    }

    emit:
    reads        = ch_reads
    sra_manifest = COLLECT_SRA_QUERY.out.manifest
}
