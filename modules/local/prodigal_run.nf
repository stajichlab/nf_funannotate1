// Prodigal ab-initio gene prediction, run standalone ahead of FUNANNOTATE_PREDICT
// (structural analogue of GENEMARK_RUN) so its gene-model GFF can be handed to
// funannotate predict as an additional EVM pass-through source via
// --other_gff <gff>:<weight>.
//
// Use case (2026-08-30): microsporidia -- tiny, intron-poor, very short-ORF
// genomes where neither Augustus nor (often) GeneMark recover enough training
// genes for funannotate's 30-model floor. Prodigal's single-genome mode is
// strong at short-ORF discovery on such compact genomes; its models then feed
// funannotate's EVM integration as source 'other_pred1' with their own weight.
// Follow-up candidate: Saccharomycotina (dense, intron-poor but larger ORFs)
// weighing prodigal models against Augustus/RNA-seq.
//
// Hierarchy fix (2026-08-31): Prodigal's native -f gff output is CDS-only,
// and EVM assembles models from gene/mRNA blocks -- a CDS-only pass-through is
// structurally inert no matter its weight. bin/prodigal_gff_hier.py re-emits
// each CDS as a gene/mRNA/CDS block (blank-line separated, no header), the
// same transform validated on the OC4 held-out gold genome: gene Sn/Sp went
// 0.584/0.698 (CDS-only, consensus was a pure GeneMark echo) -> 0.839/0.817.
// => prodigal_weight default is now 5 (was 1), calibrated from that run.
//
// Contingency (2026-09-01): params.prodigal_lineages (list of BUSCO lineage
// names, e.g. ['microsporidia_odb10']) restricts which genomes get prodigal.
// Genomes whose lineage is NOT in the list emit a 0-byte GFF/FAA so
// FUNANNOTATE_PREDICT's other_gff_ok (size>0) drops the source -- no
// --other_gff, they stay on Augustus/GeneMark. Empty list = all genomes.
//
// Skeleton status / TODOs before this is production-validated:
//   TODO(masked): decide the genome handed in. This process is fed the same
//     masked `genome_fa` as GENEMARK_RUN/FUNANNOTATE_PREDICT; prodigal will NOT
//     call genes in soft-masked (lowercase) sequence, so for repeat-poor
//     microsporidia masking is a minor concern, but for a repeat-heavy genome
//     the pass-through would silently under-call. Options: (a) run on the
//     masked genome as-is and accept under-calling in repeats (simplest,
//     EVM-compatible), or (b) run on the unmasked raw, keeping only prodigal
//     models that land outside masked intervals (needs a min-max filter vs the
//     repeatmasker BED funannotate itself writes). Start with (a).
//   TODO(contigs): prodigal -p single assumes complete genomes; for
//     short/fragmented contigs it can miss short ORFs at contig ends. -p meta
//     prefers ORF recovery over training quality on fragmented inputs --
//     expose via params.prodigal_mode and default per params.predict_frag_max_*
//     pre-flight verdict if needed.
//   TODO(validate): prodigal GFF contig names come straight from the FASTA
//     headers, so they match as long as the same genome goes in; funannotate
//     itself validates contigs + EVM-formats the file (predict_misc/
//     other1_predictions.gff3) and hard-errors on mismatches. No duplicate
//     check needed here.
//
// Outputs:
//   - ${out}.prodigal.gff3 : Prodigal GFF re-emitted as gene/mRNA/CDS blocks
//     (source will become 'other_pred1' inside funannotate, EVM-validated
//     there).
//   - ${out}.prodigal.faa  : predicted proteins (informational / future
//     validation; NOT currently wired to predict).
//
// Runs under the 'prodigal' provisioning label (prodigal is conda-solvable --
// bioconda `prodigal` -- so conda/pixi/module all supply it; no licensed or
// container-only constraints like GeneMark). Resources overridden per-profile
// by withName: '.*:PRODIGAL_RUN' in conf/profile_annotate.config.
process PRODIGAL_RUN {
    label 'prodigal'
    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), val(genome_fa)

    output:
    tuple val(meta), path("${meta.id}.prodigal.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}.prodigal.faa"), optional: true, emit: proteins

    script:
    def out           = meta.id
    def transl_table  = meta.transl_table
    def mode          = params.prodigal_mode ?: 'single'
    def prodigalLineages = (params.prodigal_lineages ?: []).collect { it.toString() }
    def skip = !prodigalLineages.isEmpty() && !prodigalLineages.contains(meta.busco?.toString())
    """
    if [ "${skip.toBoolean()}" = "true" ]; then
        # Contingency: params.prodigal_lineages (e.g. ['microsporidia_odb10'])
        # restricts which BUSCO lineages get prodigal. Non-matching genomes emit
        # a 0-byte GFF -- FUNANNOTATE_PREDICT's other_gff_ok requires size>0, so
        # it drops the source and these genomes get NO --other_gff (they stay on
        # Augustus/GeneMark). The 1:1 join in the subworkflow is preserved.
        echo "[INFO] PRODIGAL_RUN ${out}: lineage ${meta.busco} not in prodigal_lineages -- skipping (empty GFF)"
        : > "${out}.prodigal.gff3"
        : > "${out}.prodigal.faa"
        exit 0
    fi

    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome"; gzip -dc "\$GENOME_FA" > genome.fa ;;
        *)    cp "\$GENOME_FA" genome.fa ;;
    esac

    # -g uses the same meta.transl_table that drives funannotate predict
    # --table, so model translation tables always agree (microsporidia CUG
    # reassignment / alt tables must match or EVM evidence is meaningless).
    # -f gff keeps the native 'gene'/'CDS' features funannotate's
    # renameGFF/EVM validation expects. Note: prodigal emits CDS rows only --
    # bin/prodigal_gff_hier.py re-emits them as gene/mRNA/CDS blocks so EVM
    # can use the source structurally (see header comment "Hierarchy fix").
    echo "[INFO] PRODIGAL_RUN ${out}: prodigal -p ${mode} -g ${transl_table} on \$(grep -c '>' genome.fa) contigs"
    prodigal -i genome.fa -f gff -p ${mode} -g ${transl_table} \\
        -o "${out}.prodigal.raw.gff3" -a "${out}.prodigal.faa"

    if [ ! -s "${out}.prodigal.raw.gff3" ]; then
        echo "ERROR: prodigal produced no GFF for ${out}" >&2
        exit 1
    fi
    python "${workflow.projectDir}/bin/prodigal_gff_hier.py" "${out}.prodigal.raw.gff3" \\
        -o "${out}.prodigal.gff3"
    rm -f genome.fa "${out}.prodigal.raw.gff3"
    echo "[INFO] PRODIGAL_RUN ${out}: \$(grep -c \$'\\tgene\\t' "${out}.prodigal.gff3") genes predicted"
    """

    stub:
    def out = meta.id
    """
    printf '%s\\tprodigal\\tgene\\t1\\t100\\t.\\t+\\t.\\tID=prodigal_g1\\n' "${out}" > "${out}.prodigal.gff3"
    printf '%s\\tprodigal\\tmRNA\\t1\\t100\\t.\\t+\\t.\\tID=prodigal_m1;Parent=prodigal_g1\\n' "${out}" >> "${out}.prodigal.gff3"
    printf '%s\\tprodigal\\tCDS\\t1\\t100\\t58.7\\t+\\t0\\tID=prodigal_m1.cds;Parent=prodigal_m1\\n' "${out}" >> "${out}.prodigal.gff3"
    printf '>stub_1\\nMSTUB\\n' > "${out}.prodigal.faa"
    echo "[STUB] PRODIGAL_RUN ${out}"
    """
}
