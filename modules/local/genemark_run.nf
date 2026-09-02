// Standalone GeneMark step, split out of funannotate predict's internal call so
// predict itself can move onto a container. Port of BFD's modules/funannotate/
// predict/GENEMARK_RUN — see BFD nextflow/docs/GENEMARK_RUN_DESIGN.md for the
// full design/validation.
//
// Two provisioning paths, selected by params.genemark_container_mode (set true
// by conf/provision_singularity.config, false/unset under ucr_hpcc/pixi):
//   - container mode: GeneMark/Augustus-aux commands (gmes_petap.pl, bam2hints,
//     join_mult_hints.pl) run via a manually-built `apptainer exec` against
//     params.container_genemark (public teambraker/braker3 image). Built
//     manually rather than via a process `container =` directive because this
//     script also calls python/gzip on the bare host shell (not guaranteed
//     present in that image) -- see the withLabel: 'genemark' comment in
//     conf/provision_singularity.config. GeneMark's own license still requires
//     a user-obtained ~/.gm_key even from this public image (the image only
//     ships the redistributable binary, not a license) -- resolved and bound
//     defensively below, since on this cluster ~/.gm_key can be a symlink
//     landing outside $HOME/$PWD (e.g. under /opt/linux), which apptainer's
//     plain $HOME automount would not follow.
//   - host mode (ucr_hpcc/pixi): GeneMark is located in this priority order:
//       1. params.genemark_path (user-local licensed install dir, e.g. pixi/conda users)
//       2. $GENEMARK_PATH (set by `module load genemarkESET` on UCR HPCC)
//       3. `command -v gmes_petap.pl` (on PATH from the provisioning env)
//     ET-mode auxiliary tools (bam2hints / join_mult_hints.pl) come from the
//     funannotate/augustus provisioning env (module load funannotate, or the
//     pixi funannotate env) — both already on PATH there.
//
// Reuse mode (checked first, applies regardless of ES/ET — gmes_petap.pl's
// --predict_with is its own mutually-exclusive run mode, not a variant of
// --ES/--ET): shared_mod set and force_independent != 'true' ->
// `gmes_petap.pl --predict_with <mod>` — genuinely training-free prediction.
// Faster AND more correct than funannotate predict's `-p parameters.json` reuse
// (which only seeds --ini_mod into a full ES retrain).
//
// Otherwise, fresh training, mode selected by mode + training_bam:
//   - mode == 'ET' and training_bam is a non-empty existing file:
//     `gmes_petap.pl --ET`, seeded by RNA-seq-informed intron hints derived from
//     FUNANNOTATE_TRAIN's transcript-to-genome alignment BAM
//     (training/transcript.alignments.bam). Raw bam2hints output is unstranded —
//     it must run through vendored filterIntronsFindStrand.pl (Artistic License,
//     from BRAKER) to assign strand + drop non-canonical splice sites, or
//     GeneMark's branch-point step (bp_seq_select.pl) dies with "hash is empty".
//   - otherwise (mode == 'ES', or mode == 'ET' with no training BAM): fresh
//     `gmes_petap.pl --ES` self-training, identical to predict's internal
//     RunGeneMarkES().
// Both fresh-training branches produce a new .mod (candidate for backfilling to
// the shared store by BACKFILL_ABINITIO_PARAMS — this process never writes to
// the shared store itself).
//
// --genemark_gtf consumes GeneMark's raw native GTF (genemark.gtf, produced
// directly by gmes_petap.pl in its working dir) — predict.py runs
// genemark_gtf2gff3.pl on it internally, so no conversion is done here.
process GENEMARK_RUN {
    label 'genemark'
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), val(genome_fa), val(mode), val(training_bam), val(force_independent), val(shared_mod)

    output:
    tuple val(meta), path("${meta.id}.genemark.gtf"), emit: gtf
    // Keyed tuple, not a bare path: a batch can mix fresh-train rows (which
    // emit this) and fast-reuse rows (which don't — --predict_with produces no
    // new .mod). Keying by (meta, out) keeps downstream joins correct even when
    // some rows never emit this output at all.
    tuple val(meta), path("${meta.id}.genemark.mod"), optional: true, emit: mod

    script:
    def out            = meta.id
    def transl_table   = meta.transl_table
    def filterIntronsFindStrand = "${workflow.projectDir}/bin/vendor/filterIntronsFindStrand.pl"
    // params.genemark_path defaults to null; Groovy string interpolation would
    // render that as the literal text "null" below, which is non-empty and so
    // defeats the GENEMARK_PATH/command -v fallback chain entirely -- coalesce
    // to '' first so `[ -z "$GM" ]` actually trips when unset.
    def genemarkPath   = params.genemark_path ?: ''
    def containerMode  = params.genemark_container_mode ?: false
    """
    if [ "${containerMode}" = "true" ]; then
        # ── Container mode: manually build an apptainer invocation for just the
        # GeneMark/Augustus-aux commands below -- python/gzip in this same
        # script still run on the bare host shell (not guaranteed present in
        # ${params.container_genemark}), so this can't be a process
        # `container =` full-script wrap. See conf/provision_singularity.config's
        # withLabel: 'genemark' comment.
        # Node-local scratch may not exist / be writable if an inherited
        # \$SCRATCH points at another node's path — fall back to the task
        # workdir, which is already bound into the container below.
        export TMPDIR=\$(printf '%s' "\${SCRATCH:-}" | tr -d '\\n\\r')
        TMPDIR=\${TMPDIR:-/tmp}
        if [ ! -d "\$TMPDIR" ] || [ ! -w "\$TMPDIR" ]; then
            TMPDIR="\$PWD"
        fi
        export TMPDIR
        # GeneMark's license resolution (gmes_petap.pl reads ~/.gm_key) can hop
        # through host-specific symlinks that land outside \$HOME/\$PWD (e.g. a
        # module-installed license landing under this cluster's /opt/linux
        # package tree) -- resolve the real target dynamically and bind its
        # parent dir instead of hardcoding a host-specific prefix. A no-op
        # (binds \$HOME again) at any site where ~/.gm_key isn't a symlink, or
        # doesn't exist yet.
        GM_KEY_REAL=\$(readlink -f "\$HOME/.gm_key" 2>/dev/null || echo "\$HOME/.gm_key")
        GM_KEY_DIR=\$(dirname "\$GM_KEY_REAL")
        SING_BINDS="--bind \$PWD:\$PWD,${workflow.projectDir}:${workflow.projectDir},\$GM_KEY_DIR:\$GM_KEY_DIR,\$TMPDIR:\$TMPDIR"
        # training_bam (ET mode only) lives outside the task workdir -- it's
        # read directly by the containerized bam2hints call below, not staged
        # into \$PWD first, so its directory needs its own bind.
        TRAINING_BAM="${training_bam}"
        if [ -n "\$TRAINING_BAM" ]; then
            SING_BINDS="\$SING_BINDS,\$(dirname "\$TRAINING_BAM"):\$(dirname "\$TRAINING_BAM")"
        fi
        SING="apptainer exec \$SING_BINDS ${params.container_genemark}"
        echo "[INFO] GENEMARK_RUN ${out}: running GeneMark from container ${params.container_genemark}"
    else
        SING=""
        # ── Locate GeneMark: params.genemark_path > \$GENEMARK_PATH > command -v ────
        GM="${genemarkPath}"
        if [ -z "\$GM" ] && [ -n "\$GENEMARK_PATH" ] && [ -x "\$GENEMARK_PATH/gmes_petap.pl" ]; then
            GM="\$GENEMARK_PATH"
        fi
        if [ -z "\$GM" ] && command -v gmes_petap.pl >/dev/null 2>&1; then
            GM="\$(dirname "\$(command -v gmes_petap.pl)")"
        fi
        if [ -z "\$GM" ] || [ ! -x "\$GM/gmes_petap.pl" ]; then
            echo "ERROR: GeneMark-ES/ET not found. Set params.genemark_path to your licensed install, or ensure a host module/PATH provides gmes_petap.pl (e.g. 'module load genemarkESET' on UCR HPCC)." >&2
            echo "       Resolved: GENEMARK_PATH=\$GENEMARK_PATH, GENEMARK_PATH param='${genemarkPath}', gm=\$GM" >&2
            exit 1
        fi
        export GENEMARK_PATH="\$GM"

        # gmes_petap.pl needs Perl's YAML.pm. UCR HPCC's funannotate/genemarkESET
        # host modules run it under the system perl (/usr/bin/perl via #!/usr/bin/env
        # perl) and don't ship or wire in YAML.pm (no PERL5LIB in either modulefile),
        # so a bare host-module GENEMARK_RUN dies with "Can't locate YAML.pm in @INC"
        # before even printing --help. The conda envs built for the funannotate
        # provisioning axis happen to bundle a working YAML.pm; reuse it as a
        # PERL5LIB fallback, but only when the host perl doesn't already have one
        # (keeps this a no-op anywhere YAML.pm is already resolvable/installed).
        if ! perl -MYAML -e1 >/dev/null 2>&1; then
            _yaml_site="${params.conda_envs_root ?: '/bigdata/stajichlab/shared/condaenv'}/funannotate-1.8.17/lib/perl5/site_perl"
            if [ -f "\$_yaml_site/YAML.pm" ]; then
                export PERL5LIB="\$_yaml_site\${PERL5LIB:+:\$PERL5LIB}"
                echo "[INFO] GENEMARK_RUN ${out}: host perl lacks YAML.pm; PERL5LIB+=\$_yaml_site"
            fi
        fi
        echo "[INFO] GENEMARK_RUN ${out}: GENEMARK_PATH=\$GENEMARK_PATH"
    fi

    # Dispatch helpers: \$SING is empty in host mode (bare host commands);
    # in container mode they run through the apptainer exec built above.
    run_gmes() {
        if [ -n "\$SING" ]; then \$SING gmes_petap.pl "\$@"; else "\$GENEMARK_PATH/gmes_petap.pl" "\$@"; fi
    }
    run_bam2hints() {
        if [ -n "\$SING" ]; then \$SING bam2hints "\$@"; else bam2hints "\$@"; fi
    }
    run_join_hints() {
        if [ -n "\$SING" ]; then \$SING join_mult_hints.pl "\$@"; else join_mult_hints.pl "\$@"; fi
    }
    run_filter_introns() {
        if [ -n "\$SING" ]; then \$SING perl "\$@"; else perl "\$@"; fi
    }

    GENOME_GZ="${genome_fa}"
    case "\$GENOME_GZ" in
        *.gz) gzip -dc "\$GENOME_GZ" > genome.fa ;;
        *)    cp "\$GENOME_GZ" genome.fa ;;
    esac

    # ── Too-small/fragmented-genome pre-flight guard ─────────────────────────
    # Mirrors FUNANNOTATE_PREDICT's own guard (same params.predict_min_asm_bp/
    # predict_frag_max_n50/predict_frag_max_contigs, same
    # bin/asm_preflight_stats.py -- shared, not duplicated) -- GENEMARK_RUN
    # sits upstream of that guard in the DAG, so without a check here, a
    # genome predict would go on to skip anyway (small AND fragmented)
    # instead burns a full gmes_petap.pl --ES/--ET attempt and hard-fails
    # once GeneMark's own internal contig selection (--min_contig 10000,
    # after excluding soft-masked repeat sequence) leaves too little usable
    # training sequence, e.g. "error, input sequence size is too small
    # data/training.fna: 32259" for a small, fragmented assembly. Skip
    # cleanly here (empty GTF, no .mod) so predict's own preflight guard is
    # the one that actually flags/records the skip in predict_skipped_too_small.tsv.
    read ASM_BP ASM_CTG ASM_N50 ASM_VERDICT < <(
        python "${workflow.projectDir}/bin/asm_preflight_stats.py" genome.fa \\
            --min-bp ${params.predict_min_asm_bp} --max-n50 ${params.predict_frag_max_n50} \\
            --max-contigs ${params.predict_frag_max_contigs})
    if [ "\$ASM_VERDICT" = "small_fragmented" ]; then
        echo "[WARN] GENEMARK_RUN ${out}: too small/fragmented (\$ASM_BP bp, \$ASM_CTG contigs, N50 \$ASM_N50); skipping GeneMark -- predict's own preflight guard will flag/report this genome" >&2
        touch "${out}.genemark.gtf"
        rm -f genome.fa
        exit 0
    fi

    # ── --gcode support probe (mirrors funannotate's _genemark_supports_gcode()) ──
    GCODE_ARGS=()
    if [ "${transl_table}" = "6" ] || [ "${transl_table}" = "26" ]; then
        if run_gmes 2>&1 | grep -q -- '--gcode'; then
            GCODE_ARGS=(--gcode "${transl_table}")
        else
            echo "[WARN] GeneMark does not support --gcode in this version; running with default code 1 -- gene calls in CUG/alt-table genomes may be unreliable" >&2
        fi
    elif [ "${transl_table}" != "1" ]; then
        echo "[WARN] GeneMark only supports --gcode 6 or 26; ignoring table ${transl_table} and running with default code 1" >&2
    fi

    # gmes_petap.pl's own "too small" abort (e.g. "error, input sequence size
    # is too small data/training.fna: 32259") prints to stdout, not a log
    # file, and gmes_petap.pl still exits 0 in that case -- only the missing
    # output/gmhmm.mod reveals the failure. Capture stdout alongside letting
    # it stream to the task log (tee), so the branches below can distinguish
    # "GeneMark declined because the input was too small/fragmented after its
    # own internal contig selection" (graceful skip, same outcome as the
    # pre-flight guard above) from a genuine unexpected failure (hard error).
    GMES_LOG="gmes_stdout.log"
    too_small_skip() {
        grep -qi "input sequence size is too small" "\$GMES_LOG" 2>/dev/null
    }

    if [ -n "${shared_mod}" ] && [ "${force_independent}" != "true" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fast-reuse (--predict_with) against shared model ${shared_mod}"
        cp "${shared_mod}" genemark-shared.mod
        run_gmes --predict_with genemark-shared.mod \\
            --sequence genome.fa --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" 2>&1 | tee "\$GMES_LOG"
    elif [ "${mode}" = "ET" ] && [ -n "${training_bam}" ] && [ -s "${training_bam}" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fresh ET self-training seeded by RNA-seq intron hints from ${training_bam}"
        run_bam2hints --intronsonly --in="${training_bam}" --out=raw_hints.gff
        if [ ! -s raw_hints.gff ]; then
            echo "ERROR: bam2hints produced no intron hints from ${training_bam}" >&2
            exit 1
        fi
        run_filter_introns "${filterIntronsFindStrand}" genome.fa raw_hints.gff --score > stranded_hints.gff
        if [ ! -s stranded_hints.gff ]; then
            echo "ERROR: filterIntronsFindStrand.pl dropped all hints (no canonical splice sites found) -- ${training_bam} may be too noisy for --ET; consider genemark_mode=ES" >&2
            exit 1
        fi
        sort -n -k4,4 stranded_hints.gff | sort -s -n -k5,5 | sort -s -k3,3 | sort -s -k1,1 \\
            | run_join_hints > genemark.intron-hints.gff
        run_gmes --ET genemark.intron-hints.gff --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" 2>&1 | tee "\$GMES_LOG"
        if [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        elif too_small_skip; then
            echo "[WARN] GENEMARK_RUN ${out}: GeneMark-ET declined -- not enough usable (unmasked, >=10kb) training sequence after masking; skipping GeneMark for this genome" >&2
            touch "${out}.genemark.gtf"
            rm -f genome.fa
            exit 0
        else
            echo "ERROR: GeneMark-ET did not produce output/gmhmm.mod" >&2
            exit 1
        fi
    else
        if [ "${mode}" = "ET" ]; then
            echo "[INFO] GENEMARK_RUN ${out}: mode=ET requested but no training_bam available (no RNA-seq for this genome) -- falling back to ES self-training"
        fi
        echo "[INFO] GENEMARK_RUN ${out}: fresh ES self-training (force_independent=${force_independent}, shared_mod='${shared_mod}')"
        run_gmes --ES --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" 2>&1 | tee "\$GMES_LOG"
        if [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        elif too_small_skip; then
            echo "[WARN] GENEMARK_RUN ${out}: GeneMark-ES declined -- not enough usable (unmasked, >=10kb) training sequence after masking; skipping GeneMark for this genome" >&2
            touch "${out}.genemark.gtf"
            rm -f genome.fa
            exit 0
        else
            echo "ERROR: GeneMark-ES did not produce output/gmhmm.mod" >&2
            exit 1
        fi
    fi

    if [ ! -s genemark.gtf ]; then
        if too_small_skip; then
            echo "[WARN] GENEMARK_RUN ${out}: GeneMark declined -- not enough usable training sequence; skipping GeneMark for this genome" >&2
            touch "${out}.genemark.gtf"
            rm -f genome.fa
            exit 0
        fi
        echo "ERROR: GeneMark did not produce genemark.gtf" >&2
        exit 1
    fi
    cp genemark.gtf "${out}.genemark.gtf"
    rm -f genome.fa

    # ── Also persist a fresh .mod to predict_misc/, like funannotate's own
    # internal GeneMark call used to. Option B persistence, same pattern as
    # FUNANNOTATE_PREDICT writing directly into params.target/${out} -- not a
    # publishDir copy. Consumers of GENEMARK_RUN.out.mod (the wired
    # BACKFILL_ABINITIO_PARAMS join) are unaffected and keep taking priority
    # via their explicit channel value; this only backstops any consumer that
    # reads genemark.mod off disk instead (a future standalone backfill
    # sweep, or bin/backfill_abinitio_params.py's predict_misc/ fallback when
    # genemark_mod isn't passed explicitly) so it doesn't silently miss
    # genemark for anything predicted through this process. See BFD's
    # nextflow/docs/GENEMARK_RUN_DESIGN.md "Known gap" section for the
    # original report of this.
    if [ -f "${out}.genemark.mod" ]; then
        LOWER_OUT=\$(printf '%s' "${out}" | tr '[:upper:]' '[:lower:]')
        AB_INITIO_DIR="${params.target}/${out}/predict_misc/ab_initio_parameters"
        mkdir -p "\$AB_INITIO_DIR"
        cp "${out}.genemark.mod" "\$AB_INITIO_DIR/\${LOWER_OUT}.genemark.mod"
    fi

    echo "[INFO] GENEMARK_RUN complete for ${out}"
    """

    stub:
    def out = meta.id
    """
    touch "${out}.genemark.gtf"
    if [ -z "${shared_mod}" ] || [ "${force_independent}" = "true" ]; then
        touch "${out}.genemark.mod"
    fi
    """
}
