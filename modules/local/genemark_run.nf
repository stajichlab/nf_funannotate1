// Standalone GeneMark step, split out of funannotate predict's internal call so
// predict itself can move onto a container (GeneMark's license forbids
// redistribution, so it can never be containerized and must stay on a host
// GeneMark install). Port of BFD's modules/funannotate/predict/GENEMARK_RUN —
// see BFD nextflow/docs/GENEMARK_RUN_DESIGN.md for the full design/validation.
//
// GeneMark is located in this priority order:
//   1. params.genemark_path (user-local licensed install dir, e.g. pixi/conda users)
//   2. $GENEMARK_PATH (set by `module load genemarkESET` on UCR HPCC)
//   3. `command -v gmes_petap.pl` (on PATH from the provisioning env)
// ET-mode auxiliary tools (bam2hints / join_mult_hints.pl) come from the
// funannotate/augustus provisioning env (module load funannotate, or the pixi
// funannotate env) — both already on PATH there.
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
    """
    # ── Locate GeneMark: params.genemark_path > \$GENEMARK_PATH > command -v ────
    GM="${params.genemark_path}"
    if [ -z "\$GM" ] && [ -n "\$GENEMARK_PATH" ] && [ -x "\$GENEMARK_PATH/gmes_petap.pl" ]; then
        GM="\$GENEMARK_PATH"
    fi
    if [ -z "\$GM" ] && command -v gmes_petap.pl >/dev/null 2>&1; then
        GM="\$(dirname "\$(command -v gmes_petap.pl)")"
    fi
    if [ -z "\$GM" ] || [ ! -x "\$GM/gmes_petap.pl" ]; then
        echo "ERROR: GeneMark-ES/ET not found. Set params.genemark_path to your licensed install, or ensure a host module/PATH provides gmes_petap.pl (e.g. 'module load genemarkESET' on UCR HPCC)." >&2
        echo "       Resolved: GENEMARK_PATH=\$GENEMARK_PATH, GENEMARK_PATH param='${params.genemark_path}', gm=\$GM" >&2
        exit 1
    fi
    export GENEMARK_PATH="\$GM"
    echo "[INFO] GENEMARK_RUN ${out}: GENEMARK_PATH=\$GENEMARK_PATH"

    GENOME_GZ="${genome_fa}"
    case "\$GENOME_GZ" in
        *.gz) pigz -dc "\$GENOME_GZ" > genome.fa ;;
        *)    cp "\$GENOME_GZ" genome.fa ;;
    esac

    # ── --gcode support probe (mirrors funannotate's _genemark_supports_gcode()) ──
    GCODE_ARGS=()
    if [ "${transl_table}" = "6" ] || [ "${transl_table}" = "26" ]; then
        if "\$GENEMARK_PATH/gmes_petap.pl" 2>&1 | grep -q -- '--gcode'; then
            GCODE_ARGS=(--gcode "${transl_table}")
        else
            echo "[WARN] GeneMark does not support --gcode in this version; running with default code 1 -- gene calls in CUG/alt-table genomes may be unreliable" >&2
        fi
    elif [ "${transl_table}" != "1" ]; then
        echo "[WARN] GeneMark only supports --gcode 6 or 26; ignoring table ${transl_table} and running with default code 1" >&2
    fi

    if [ -n "${shared_mod}" ] && [ "${force_independent}" != "true" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fast-reuse (--predict_with) against shared model ${shared_mod}"
        cp "${shared_mod}" genemark-shared.mod
        "\$GENEMARK_PATH/gmes_petap.pl" --predict_with genemark-shared.mod \\
            --sequence genome.fa --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}"
    elif [ "${mode}" = "ET" ] && [ -n "${training_bam}" ] && [ -s "${training_bam}" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fresh ET self-training seeded by RNA-seq intron hints from ${training_bam}"
        bam2hints --intronsonly --in="${training_bam}" --out=raw_hints.gff
        if [ ! -s raw_hints.gff ]; then
            echo "ERROR: bam2hints produced no intron hints from ${training_bam}" >&2
            exit 1
        fi
        perl "${filterIntronsFindStrand}" genome.fa raw_hints.gff --score > stranded_hints.gff
        if [ ! -s stranded_hints.gff ]; then
            echo "ERROR: filterIntronsFindStrand.pl dropped all hints (no canonical splice sites found) -- ${training_bam} may be too noisy for --ET; consider genemark_mode=ES" >&2
            exit 1
        fi
        sort -n -k4,4 stranded_hints.gff | sort -s -n -k5,5 | sort -s -k3,3 | sort -s -k1,1 \\
            | join_mult_hints.pl > genemark.intron-hints.gff
        "\$GENEMARK_PATH/gmes_petap.pl" --ET genemark.intron-hints.gff --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}"
        if [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        else
            echo "ERROR: GeneMark-ET did not produce output/gmhmm.mod" >&2
            exit 1
        fi
    else
        if [ "${mode}" = "ET" ]; then
            echo "[INFO] GENEMARK_RUN ${out}: mode=ET requested but no training_bam available (no RNA-seq for this genome) -- falling back to ES self-training"
        fi
        echo "[INFO] GENEMARK_RUN ${out}: fresh ES self-training (force_independent=${force_independent}, shared_mod='${shared_mod}')"
        "\$GENEMARK_PATH/gmes_petap.pl" --ES --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}"
        if [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        else
            echo "ERROR: GeneMark-ES did not produce output/gmhmm.mod" >&2
            exit 1
        fi
    fi

    if [ ! -s genemark.gtf ]; then
        echo "ERROR: GeneMark did not produce genemark.gtf" >&2
        exit 1
    fi
    cp genemark.gtf "${out}.genemark.gtf"
    rm -f genome.fa
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
