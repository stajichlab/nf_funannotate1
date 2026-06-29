#!/usr/bin/env nextflow

/*
 * earlgrey_mask — curated repeat masking for large fungal genomes
 *
 * Runs EarlGrey once per species on the single best representative genome
 * (>cutoff_mb), then applies the resulting curated TE library with RepeatMasker to
 * every conspecific strain. Output is a soft-masked input_clean_genomes/<asmid>.masked.fasta
 * per strain — exactly the file funannotate.nf consumes (its tantan MASKREPEAT_TANTAN_RUN
 * uses storeDir on input_clean_genomes, so it is skipped wherever this file already exists).
 *
 * Runs independently of funannotate.nf. Submit from the project root.
 *
 * Usage:
 *   nextflow run nextflow/earlgrey_mask.nf \
 *       -c nextflow/nextflow.config -profile earlgrey -resume
 *
 *   # limit to N species for testing
 *   nextflow run nextflow/earlgrey_mask.nf \
 *       -c nextflow/nextflow.config -profile earlgrey --n_test 1 -resume
 *
 * Stub/dry-run:
 *   nextflow run nextflow/earlgrey_mask.nf \
 *       -c nextflow/nextflow.config -profile earlgrey -stub-run --n_test 2
 */

// ════════════════════════════════════════════════════════════════════════════
// PARAMETER DEFAULTS  (override via --param value)
// ════════════════════════════════════════════════════════════════════════════

params.genome_dir          = "${launchDir}/input_clean_genomes"
params.genome_suffix       = '.fa'                       // clean (unmasked) genome suffix
params.tables_dir          = "${launchDir}/tables"       // where asm_stats.tsv.gz lives
params.asm_stats           = "${params.tables_dir}/asm_stats.tsv.gz"
params.gen_asm_stats       = true                        // generate asm_stats if missing
params.skip_select_reps    = false                       // skip SELECT_REPS step (just do EarlGrey on all)
params.cutoff_mb           = 200                         // species qualifies if rep > this
params.repeat_taxon        = 'fungi'                     // EarlGrey -r RepeatMasker search term
params.earlgrey_version    = '7.2.6'
params.repeatmasker_version = '4.1.8'
params.outdir              = "${launchDir}/results/repeatlibrary"
params.earlgrey_outdir     = "${params.outdir}"                    // alias used by the shared EarlGrey modules
params.masked_dir          = "${launchDir}/input_clean_genomes"    // where <asmid>.masked.fasta.gz land
params.earlgrey_workdir    = "${launchDir}/work/earlgrey_persist"  // persistent per-species EarlGrey output (enables resume)

// genomeFile: see lib/FunannotateUtils.groovy (shared with funannotate.nf and subworkflows).
def genomeFile(String base) { FunannotateUtils.genomeFile(base) }

// ════════════════════════════════════════════════════════════════════════════
// INCLUDES
// ════════════════════════════════════════════════════════════════════════════

include { ASM_STATS          } from './modules/local/asm_stats'
include { INPUT_CHECK        } from './subworkflows/local/input_check'
include { EARLGREY_BUILD_LIB } from './modules/local/earlgrey_build_lib'
include { REPEATMASK_STRAIN  } from './modules/local/repeatmask_strain'
include { DELIVER_MASK       } from './modules/local/deliver_mask'

// ════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ════════════════════════════════════════════════════════════════════════════

// ── SELECT_REPS ──────────────────────────────────────────────────────────────
// Joins samples.csv + asm_stats to pick the best representative per species and
// keep only species whose representative exceeds the size cutoff. storeDir caches
// the CSV in misc/ so re-runs skip this step (delete the file to refresh).
process SELECT_REPS {
    tag    'select'
    label  'select'

    storeDir "${launchDir}/misc"

    input:
        path samples
        path asm_stats

    output:
        path 'repeat_representatives.csv'

    script:
    """
    python ${launchDir}/scripts/select_repeat_representatives.py \\
        --samples ${samples} \\
        --asm-stats ${asm_stats} \\
        --genome-dir ${params.genome_dir} \\
        --genome-suffix ${params.genome_suffix} \\
        --cutoff-mb ${params.cutoff_mb} \\
        --output repeat_representatives.csv
    """

    stub:
    """
    printf 'SPECIES,REP_ASMID,REP_SIZE_MB,N_MEMBERS,MEMBER_ASMIDS\\n' > repeat_representatives.csv
    printf 'Stub species one,STUBASM_REP1,250.0,1,STUBASM_MEM1\\n'   >> repeat_representatives.csv
    printf 'Stub species two,STUBASM_REP2,300.0,0,\\n'              >> repeat_representatives.csv
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {

    // ── Generate assembly statistics if needed ────────────────────────────────
    // ASM_STATS generates asm_stats.tsv.gz from clean genomes (used by SELECT_REPS).
    // Only runs if gen_asm_stats=true and the file doesn't already exist.
    if (params.gen_asm_stats.toBoolean()) {
        def asm_stats_path = file(params.tables_dir).toAbsolutePath()
        def asm_stats_gz = file("${asm_stats_path}/asm_stats.tsv.gz")
        if (!asm_stats_gz.exists()) {
            log.info "Generating assembly statistics: ${asm_stats_gz}"
            ASM_STATS(
                file(params.samples, glob: false),
                file(params.genome_dir, glob: false)
            )
        } else {
            log.info "Assembly statistics already exist: ${asm_stats_gz}"
        }
    }

    // ── Select representatives (skip with --skip_select_reps) ────────────────
    // When skip_select_reps=true, all genomes are processed for EarlGrey
    // without the size/N50 filtering applied by SELECT_REPS.
    def reps
    if (params.skip_select_reps.toBoolean()) {
        // INPUT_CHECK applies the same taxon/asmid/suppress/n_test filters as funannotate.nf.
        // We use samples (meta-only, no genome-existence check) because EarlGrey reads genomes
        // from input_clean_genomes/, not from the raw params.source path.
        log.info "Skipping SELECT_REPS; processing all filtered genomes for EarlGrey"
        INPUT_CHECK()
        reps = INPUT_CHECK.out.samples
            .map { meta ->
                "SPECIES,REP_ASMID,REP_SIZE_MB,N_MEMBERS,MEMBER_ASMIDS\n${meta.species},${meta.asmid},0.0,0,"
            }
            .collectFile(name: "${launchDir}/misc/repeat_representatives.csv", newLine: false)
    } else {
        reps = SELECT_REPS(
            file(params.samples,   glob: false),
            file(params.asm_stats, glob: false),
        )
    }

    // ── Per-species records (n_test limits *species*) ─────────────────────────
    def records = reps
        .splitCsv(header: true)
        .take(params.n_test > 0 ? params.n_test as int : -1)
        .map { row ->
            tuple(row.SPECIES?.trim(), row.REP_ASMID?.trim(), (row.MEMBER_ASMIDS ?: '').trim())
        }

    // ── EarlGrey library build on the representative ──────────────────────────
    def build_in = records.map { species, rep_asmid, _members ->
        def g = genomeFile("${params.genome_dir}/${rep_asmid}${params.genome_suffix}")
        if (!g.exists() && !workflow.stubRun) {
            log.warn "Skipping ${species}: representative genome not found at ${g}"
            return null
        }
        tuple(species, rep_asmid, g)
    }.filter { it != null }

    EARLGREY_BUILD_LIB(build_in)

    // ── Expand conspecific members and pair each with its species library ─────
    def members_ch = records.flatMap { species, _rep, members ->
        (members ? members.split(';') : [])
            .findAll { it?.trim() }
            .collect { asm -> tuple(species, asm.trim()) }
    }

    def lib_by_species = EARLGREY_BUILD_LIB.out.lib
        .map { species, _rep, lib -> tuple(species, lib) }

    def mask_in = members_ch
        .combine(lib_by_species, by: 0)     // tuple(species, asmid, library)
        .map { species, asmid, library ->
            def g = genomeFile("${params.genome_dir}/${asmid}${params.genome_suffix}")
            if (!g.exists() && !workflow.stubRun) {
                log.warn "Skipping member ${asmid} (${species}): genome not found at ${g}"
                return null
            }
            tuple(asmid, species, g, library)
        }
        .filter { it != null }

    REPEATMASK_STRAIN(mask_in)

    // ── Deliver all soft-masked genomes into input_clean_genomes ──────────────
    // Both the EarlGrey-masked representative and the RepeatMasked members are
    // copied in (overwriting any older tantan mask). Each is tuple(asmid, masked).
    def masked_ch = EARLGREY_BUILD_LIB.out.masked
        .mix(REPEATMASK_STRAIN.out)

    DELIVER_MASK(masked_ch)
}
