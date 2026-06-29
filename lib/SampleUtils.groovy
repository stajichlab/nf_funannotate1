/**
 * SampleUtils — shared helpers for constructing filesystem-safe sample identifiers
 * from samples.csv SPECIES and STRAIN columns.
 *
 * Placed in nextflow/lib/ so it is automatically on the Groovy classpath for all
 * workflow scripts in nextflow/*.nf.
 *
 * The Python equivalent is collect_busco_stats.py::build_basename_map() and
 * fix_low_trinity.py::species_key_from_row().  Keep them in sync.
 */
class SampleUtils {

    /**
     * Build a filesystem-safe "{species}_{strain}" tag from raw samples.csv values.
     *
     * Canonicalises:
     *   - strips leading/trailing whitespace and quote characters (', ") from both fields
     *   - takes only the first semicolon-delimited token of strain (some rows list
     *     multiple synonymous strains separated by ';')
     *   - replaces colons with spaces in strain (colons appear as ' colon ' separators)
     *   - collapses runs of whitespace, /, #, [, ], *, ?, {, } into single underscores
     *
     * Examples:
     *   makeSampleTag("Saccharomyces cerevisiae", "CBS 1171")    → "Saccharomyces_cerevisiae_CBS_1171"
     *   makeSampleTag("Aspergillus fumigatus", "Af293; CBS 101")  → "Aspergillus_fumigatus_Af293"
     *   makeSampleTag("Fusarium oxysporum", "")                   → "Fusarium_oxysporum"
     */
    static String makeSampleTag(String rawSpecies, String rawStrain) {
        def sp = (rawSpecies ?: '').trim().replaceAll(/['"]/, '')
        def st = (rawStrain  ?: '').trim()
                    .replaceAll(/['"]/, '')
                    .split(';')[0]
                    .trim()
                    .replace(':', ' ')
        return [sp, st].findAll { it }
                       .join('_')
                       .replaceAll(/[\s\/\#\[\]\*\?\{\}]+/, '_')
    }

    /**
     * Build the canonical per-sample `meta` map from a raw samples.csv row (the
     * map produced by splitCsv(header: true)).
     *
     * This is the data contract for the DSL2 modularization (REFACTORING_PLAN.md
     * Principle 0): channels carry `tuple val(meta), val(genome)` and `meta.id`
     * is the ONLY field used for tag{}/naming. It reproduces, field-for-field,
     * the cleaning the funannotate.nf `jobs` channel does today, so wiring it in
     * is a behaviour-preserving swap.
     *
     * NOTE: `header_length` is intentionally NOT a meta field — it is a constant
     * (params.header_length, default 24), not per-sample payload. The genome path
     * also travels separately as the 2nd tuple element, not inside meta.
     *
     * Not yet wired into the workflow; introduced ahead of the atomic channel
     * conversion so the contract has one authoritative definition.
     */
    static Map makeMeta(Map row) {
        def species  = (row.SPECIES?.trim() ?: '').replaceAll(/['"]/, '')
        def strain   = (row.STRAIN?.trim() ?: '').replaceAll(/['"]/, '')
                          .replaceAll(/;.*$/, '').trim().replace(':', ' ')
        return [
            id          : makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: ''),
            asmid       : row.ASMID?.trim(),
            species     : species,
            strain      : strain,
            locustag    : row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim(),
            busco       : row.BUSCO_LINEAGE?.trim(),
            transl_table: row.TRANSL_TABLE?.trim() ?: '1',
            taxonid     : row.NCBI_TAXONID?.trim(),
        ]
    }
}
