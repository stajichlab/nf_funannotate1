// Run funannotate predict for one assembly. Writes directly into params.target/<id>/
// (Option B persistence: no publishDir copy). Emits a small marker file to carry the
// DAG edge without transferring the full predict tree through Nextflow's work/ directory.
// Resources overridden by withName: '.*:FUNANNOTATE_PREDICT' in conf/profile_annotate.config.
process FUNANNOTATE_PREDICT {
    label 'funannotate'
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), val(genome_fa), val(genemark_gtf), val(other_gff)

    output:
    val meta, emit: metadata
    path("${meta.id}.predict.done"), emit: done

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    def header_length = params.header_length
    def transl_table  = meta.transl_table
    // GeneMark GTF supplied by the standalone GENEMARK_RUN step; empty string
    // means "let funannotate run GeneMark internally (or auto-skip it)".
    // Checked by actual file size, not by Groovy truthiness on genemark_gtf
    // itself: GENEMARK_RUN's too-small-genome skip path emits a real but
    // deliberately empty ${out}.genemark.gtf (see genemark_run.nf), and a
    // non-null Path is always truthy in Groovy regardless of its size --
    // without this check every skip would still pass --genemark_gtf
    // <0-byte file> to funannotate.
    def genemark_gtf_file = genemark_gtf ? file(genemark_gtf as String) : null
    def genemark_gtf_ok   = genemark_gtf_file && genemark_gtf_file.exists() && genemark_gtf_file.size() > 0
    def genemark_cli      = genemark_gtf_ok ? "--genemark_gtf ${genemark_gtf}" : "--auto-skip-genemark"
    // -w genemark:1 MUST be passed explicitly whenever genemark_gtf is used,
    // in the SAME -w group as codingquarry:0/glimmerhmm:0 (funannotate's
    // argparse -w/--weights is nargs='+' without action='append', so a
    // second -w occurrence replaces the whole list rather than merging it --
    // ported from BFD's FUNANNOTATE_PREDICT/main.nf, which hit this).
    // Needed because none of this pipeline's funannotate provisioning modes
    // (module/pixi/singularity funannotate env -- see conf/provision_*.config)
    // put gmes_petap.pl on PATH inside this task's environment: GeneMark now
    // runs standalone in GENEMARK_RUN, so funannotate predict's own
    // genemarkcheck is always False here, and predict.py unconditionally
    // zeroes StartWeights["genemark"] when that's the case -- with NO check
    // for whether --genemark_gtf was supplied as an alternative.
    def weight_args = genemark_gtf_ok ? 'codingquarry:0 glimmerhmm:0 genemark:1' : 'codingquarry:0 glimmerhmm:0'
    // Optional external gene-model pass-through (PRODIGAL_RUN output). Passed as
    // --other_gff <gff>:<weight>; funannotate renames the source to
    // 'other_pred1', EVM-validates it, and applies the weight as
    // StartWeights['other_pred1'] (see predict.py --other_gff parsing). Empty
    // string -> omitted entirely.
    // PRODIGAL_RUN now emits gene/mRNA/CDS blocks (prodigal_gff_hier.py): EVM
    // assembles consensus models from gene/mRNA blocks, so a CDS-only source
    // is structurally inert no matter its weight (OC4 gene Sn 0.584 -> 0.839).
    def other_gff_file = other_gff ? file(other_gff as String) : null
    def other_gff_ok   = other_gff_file && other_gff_file.exists() && other_gff_file.size() > 0
    def other_gff_cli  = other_gff_ok ? "--other_gff ${other_gff}:${params.prodigal_weight}" : ''
    """
    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    # Node-local scratch may not exist / be writable if an inherited \$SCRATCH
    # points at another node's path — fall back to the task workdir.
    TMPDIR=\$(printf '%s' "\${SCRATCH:-}" | tr -d '\\n\\r')
    TMPDIR=\${TMPDIR:-/tmp}
    if [ ! -d "\$TMPDIR" ] || [ ! -w "\$TMPDIR" ]; then
        TMPDIR="\$PWD"
    fi

    PREDICTDIR="${params.target}/${out}"
    PREDICT_GBK="\$PREDICTDIR/predict_results/${out}.gbk"

    if [ "${params.debug.toBoolean()}" = "true" ]; then
        echo "[DEBUG] out=${out} asmid=${asmid} species=${species} strain=${strain}"
        echo "[DEBUG] locustag=${locustag} busco=${busco_lineage} transl_table=${transl_table}"
        echo "[DEBUG] proteins=${params.proteins} genome_fa=${genome_fa}"
        echo "[DEBUG] PREDICTDIR=\$PREDICTDIR TMPDIR=\$TMPDIR pwd=\$(pwd)"
    fi

    # ── Skip vs. refresh decision ─────────────────────────────────────────────
    if [ -s "\$PREDICT_GBK" ]; then
        SPECIES_TAG=\$(printf '%s' "${species}" | sed -E 's/[[:space:]]+/_/g')
        STALE=0
        for f in "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_R1.fastq.gz" \\
                 "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_SE.fastq.gz" \\
                 "${launchDir}/rnaseq_data/\${SPECIES_TAG}.trinity-GG.fasta"; do
            if [ -s "\$f" ] && [ "\$f" -nt "\$PREDICT_GBK" ]; then STALE=1; fi
        done
        if [ "\$STALE" -eq 0 ]; then
            echo "[INFO] Prediction already complete and current for ${out}; nothing to do"
            touch ${out}.predict.done
            exit 0
        fi
        echo "[INFO] Stale prediction for ${out}: rnaseq/trinity newer than GBK — clearing predict outputs for a fresh run"
        rm -rf "\$PREDICTDIR/predict_results" "\$PREDICTDIR/predict_misc"
    fi

    mkdir -p "\$PREDICTDIR"

    # ── Guard against a corrupt partial from a previous attempt ───────────────
    if [ ! -d "\$PREDICTDIR/predict_misc" ] && [ -d "\$PREDICTDIR/predict_results" ]; then
        echo "[WARN] predict_results/ present without predict_misc/ for ${out}; clearing stale partial"
        rm -rf "\$PREDICTDIR/predict_results"
    fi

    # Point funannotate at the persistent training dir via symlink.
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "\$PREDICTDIR/training"
    fi

    TBL2ASN_PARAMS="-l paired-ends"

    # Inflate a gzipped clean/masked genome to a local uncompressed copy.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # funannotate predict rejects FASTA deflines longer than 24 chars; NCBI-style
    # headers survive the AAFTF clean verbatim (scripts/clean_genome_fa.py keeps
    # headers untouched), so rewrite each header to its accession (first
    # whitespace token). Idempotent -- safe on already-short headers too.
    awk '/^>/{print \$1; next} {print}' "\$GENOME_IN" > "\$GENOME_IN.hdr" && mv "\$GENOME_IN.hdr" "\$GENOME_IN"

    # ── Too-small-genome pre-flight guard ────────────────────────────────────
    # Shared with GENEMARK_RUN, which needs the identical policy upstream of
    # this process (see genemark_run.nf, bin/asm_preflight_stats.py).
    SKIP_REPORT="${params.target}/predict_skipped_too_small.tsv"
    read ASM_BP ASM_CTG ASM_N50 ASM_VERDICT < <(
        python "${workflow.projectDir}/bin/asm_preflight_stats.py" "\$GENOME_IN" \\
            --min-bp ${params.predict_min_asm_bp} --max-n50 ${params.predict_frag_max_n50} \\
            --max-contigs ${params.predict_frag_max_contigs})
    echo "[INFO] Pre-flight assembly stats for ${out}: \${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}"
    if [ "\$ASM_VERDICT" = "small_fragmented" ]; then
        echo "[WARN] ${out} is too small/fragmented for funannotate training (\${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}); skipping predict" >&2
        mkdir -p "${params.target}"
        [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "preflight_small_fragmented" "\$ASM_BP" "\$ASM_CTG" "\$ASM_N50" >> "\$SKIP_REPORT"
        touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
        touch ${out}.predict.done
        exit 0
    fi

    funannotate predict --name ${locustag} -i "\$GENOME_IN" --strain "${strain}" \\
        -o "\$PREDICTDIR" -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w ${weight_args} \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence ${params.proteins} \\
        --max_intronlen ${params.max_intronlen} --min_intronlen ${params.min_intronlen} \\
        --tbl2asn "\$TBL2ASN_PARAMS" ${genemark_cli} ${other_gff_cli} || true

    # ── Post-predict catch ────────────────────────────────────────────────────
    if [ ! -s "\$PREDICT_GBK" ]; then
        PLOG="\$PREDICTDIR/logfiles/funannotate-predict.log"
        if [ -f "\$PLOG" ] && grep -q "Not enough gene models .* to train Augustus" "\$PLOG"; then
            NMODELS=\$(grep -oE "Not enough gene models [0-9]+" "\$PLOG" | grep -oE "[0-9]+" | tail -1)
            echo "[WARN] ${out}: funannotate found only \${NMODELS:-<min} training models (needs 30); too small/fragmented to annotate — skipping" >&2
            mkdir -p "${params.target}"
            [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "funannotate_too_few_models:\${NMODELS:-NA}" "" "" "" >> "\$SKIP_REPORT"
            touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
            touch ${out}.predict.done
            exit 0
        fi
        echo "ERROR: funannotate predict did not produce expected GBK: \$PREDICT_GBK" >&2
        exit 1
    fi
    if [ -d "\$PREDICTDIR/predict_misc/ab_initio_parameters" ]; then
        mv "\$PREDICTDIR/predict_misc/ab_initio_parameters" "\$PREDICTDIR"
        mv "\$PREDICTDIR/predict_misc/trnascan.no-overlaps.gff3" "\$PREDICTDIR"
        rm -rf "\$PREDICTDIR/predict_misc"
        mkdir -p "\$PREDICTDIR/predict_misc"
        mv "\$PREDICTDIR/ab_initio_parameters" "\$PREDICTDIR/trnascan.no-overlaps.gff3" "\$PREDICTDIR/predict_misc"
    fi
    find "\$PREDICTDIR/predict_results/" -maxdepth 1 \\( -name "*.txt" -o -name "*.mrna-transcripts.fa" \\) -print0 \
        | xargs -0 --no-run-if-empty pigz
    sync
    touch ${out}.predict.done
    echo "[INFO] Prediction complete for ${out} at \$PREDICTDIR"
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] Would run funannotate predict for ${out} using ${genome_fa}"
    [ -f "${genome_fa}" ] || [ -f "${genome_fa}.gz" ] || { echo "ERROR: genome not found at ${genome_fa}[.gz]" >&2; exit 1; }
    mkdir -p ${params.target}/${out}/predict_results ${params.target}/${out}/predict_misc
    echo "LOCUS stub_${out}" > ${params.target}/${out}/predict_results/${out}.gbk
    echo ">stub_${out}_p1" > ${params.target}/${out}/predict_results/${out}.proteins.fa
    touch ${out}.predict.done
    """
}
