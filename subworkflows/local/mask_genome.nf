/*
 * MASK_GENOME — genome soft-masking with three mutually-exclusive paths:
 *
 *   run_earlgrey=true    → EARLGREY_BUILD_LIB (one rep per species) →
 *                          REPEATMASK_STRAIN (all other strains) →
 *                          DELIVER_MASK (copy into params.masked_dir)
 *
 *   run_repeatmasker=true (default) → MASKREPEAT_TANTAN_RUN (storeDir-cached)
 *
 *   neither             → reuse pre-existing masked genome from params.masked_dir
 *                         if it exists, else fall through to the unmasked clean genome.
 *
 * EarlGrey path notes:
 *   • First assembly per species (by channel arrival order) is the representative;
 *     all others run REPEATMASK_STRAIN with the curated library.
 *   • Single-assembly species still run EARLGREY_BUILD_LIB and skip REPEATMASK_STRAIN.
 *   • run_earlgrey takes precedence over run_repeatmasker.
 *
 * Emits:
 *   genomes — channel: tuple(val(meta), val(genome_fa))
 *             genome_fa is an absolute-path String consumed by TRAIN_PREDICT.
 */

include { MASKREPEAT_TANTAN_RUN } from './../../modules/local/maskrepeat_tantan_run'
include { EARLGREY_BUILD_LIB    } from './../../modules/local/earlgrey_build_lib'
include { REPEATMASK_STRAIN     } from './../../modules/local/repeatmask_strain'
include { DELIVER_MASK          } from './../../modules/local/deliver_mask'

workflow MASK_GENOME {

    take:
    ch_clean  // channel: tuple(val(meta), val(genome_fa_string))

    main:
    def ch_predict

    if (params.run_earlgrey.toBoolean()) {
        // ── Key each assembly by species (filesystem-safe name) ──────────────
        def keyed = ch_clean.map { meta, gfa ->
            tuple(meta.species.replaceAll(/[^A-Za-z0-9._-]+/, '_'), meta, gfa)
        }

        // ── Group per species; first assembly is the representative ──────────
        def by_species = keyed.groupTuple(by: 0)

        // ── EARLGREY_BUILD_LIB: one per species on the representative genome ─
        def rep_in = by_species.map { sp, metas, gfas ->
            tuple(sp, metas[0].asmid, file(gfas[0]))
        }
        EARLGREY_BUILD_LIB(rep_in)

        // ── REPEATMASK_STRAIN: all non-representative assemblies ─────────────
        // lib_ch: (sp, families.fa)
        def lib_ch = EARLGREY_BUILD_LIB.out.lib
            .map { sp, _rep_asmid, lib_fa -> tuple(sp, lib_fa) }

        // Flatten members (index >= 1) from the grouped tuples
        def members_in = by_species
            .flatMap { sp, metas, gfas ->
                (1..<metas.size()).collect { i ->
                    tuple(sp, metas[i].asmid, file(gfas[i]))
                }
            }
            // join on species so each member gets the right curated library
            .combine(lib_ch, by: 0)
            .map { sp, asmid, genome, lib_fa -> tuple(asmid, sp, genome, lib_fa) }

        REPEATMASK_STRAIN(members_in)

        // ── DELIVER_MASK: copy rep + member masks into params.masked_dir ─────
        def all_masked = EARLGREY_BUILD_LIB.out.masked
            .mix(REPEATMASK_STRAIN.out)
        DELIVER_MASK(all_masked)

        // ── Reconstruct meta → emit as masked path strings ───────────────────
        // Use the actual work-dir path from DELIVER_MASK's output (not params.masked_dir);
        // publishDir copies the file to masked_dir as a side-effect for future runs.
        def asmid_to_meta = ch_clean.map { meta, _gfa -> tuple(meta.asmid, meta) }

        ch_predict = DELIVER_MASK.out.delivered
            .join(asmid_to_meta)
            .map { asmid, masked_path, meta ->
                tuple(meta, masked_path.toAbsolutePath().toString())
            }

    } else if (params.run_repeatmasker.toBoolean()) {
        MASKREPEAT_TANTAN_RUN(ch_clean)
        ch_predict = MASKREPEAT_TANTAN_RUN.out.masked
            .map { meta, masked_fa -> tuple(meta, masked_fa.toAbsolutePath().toString()) }

    } else {
        // --run_repeatmasker false / --run_earlgrey false:
        // reuse masked genome from a prior run if present in params.masked_dir.
        ch_predict = ch_clean
            .map { meta, genome_fa ->
                def fgz    = file("${params.masked_dir}/${meta.asmid}.masked.fasta.gz")
                def fa     = file("${params.masked_dir}/${meta.asmid}.masked.fasta")
                def masked = (fgz.exists() && fgz.size() > 0) ? fgz : fa
                def use_fa = masked.exists() ? masked.toString() : genome_fa
                if (params.debug.toBoolean()) {
                    log.info "[DEBUG] ${meta.asmid}: genome_fa=${use_fa} (masked=${masked.exists()})"
                }
                tuple(meta, use_fa)
            }
    }

    emit:
    genomes = ch_predict
}
