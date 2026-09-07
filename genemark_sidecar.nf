#!/usr/bin/env nextflow

/*
 * genemark_sidecar — run GeneMark exactly once per genome, out of band, for
 * consumption by multiple funannotate comparison cells.
 *
 * Motivation (Funannotate_benchmarking/DESIGN.md "GeneMark sidecar"): when
 * comparing funannotate versions/environments (funannotate.nf run repeatedly
 * against the same genome set with different --conda_env / --container_funannotate
 * / --evm_backend), GeneMark itself is not part of what's being compared and
 * should be held constant across every cell rather than retrained
 * independently N times per genome (redundant compute, and a source of
 * incidental variance/inconsistency between cells if their GeneMark path
 * differs — e.g. a host module vs. a container image).
 *
 * This script runs GENEMARK_RUN (the same module funannotate.nf's
 * TRAIN_PREDICT uses internally) once per genome, always in container mode
 * (params.genemark_container_mode is forced true below — same
 * teambraker/braker3 image every genome, independent of which provisioning
 * axis a downstream cell uses). Results are published as
 * <target>/<out>.genemark.gtf (+ .mod), matching the naming convention
 * funannotate.nf's TRAIN_PREDICT reads back via --genemark_sidecar_dir.
 *
 * ES vs. ET: funannotate.nf's own TRAIN_PREDICT feeds GeneMark --ET RNA-seq
 * intron hints from training/transcript.alignments.bam, an artifact of that
 * CELL's own `funannotate train` (PASA) run -- which is exactly what's under
 * comparison (perl vs. rust PASA/Trinity/EVM) for some of the 6 cells, so
 * reusing any one cell's BAM here would make the "shared, held-constant"
 * GeneMark result silently depend on that one cell's version/backend, and
 * every other cell's GeneMark-ET result would inherit that dependency too.
 * Instead, when `--rnaseq_reads_dir` is given, ALIGN_RNASEQ_HINTS below does
 * its own plain hisat2 spliced alignment of that species' already-fetched,
 * already-normalized reads against the masked genome -- independent of any
 * cell's train/PASA/EVM backend -- and that BAM feeds GENEMARK_RUN --ET.
 * Genomes/species with no reads in --rnaseq_reads_dir (or when the param is
 * left unset) fall back to --ES, same as before.
 *
 * Usage:
 *   nextflow run genemark_sidecar.nf -c nextflow.config \
 *       -profile genemark_sidecar,slurm,singularity \
 *       --samples samples.csv --target runs/genemark_sidecar/output \
 *       --rnaseq_reads_dir runs/v1.8.17_conda/rnaseq_reads -resume
 *
 * Then point every comparison cell at the result:
 *   nextflow run funannotate.nf ... --genemark_sidecar_dir runs/genemark_sidecar/output
 *
 * `--samples` only needs GENOME (already-masked FASTA) resolved per row —
 * the same samples.csv any funannotate.nf cell already uses is sufficient;
 * this script does not re-mask or re-clean, it only reads the GENOME column
 * via INPUT_CHECK. `--rnaseq_reads_dir` (optional) points at any ONE cell's
 * already-fetched rnaseq_reads/ dir (<species_tag>_norm_R1/_R2/_SE.fastq.gz,
 * FETCH_RNASEQ:SRA_FETCH's output naming) -- the raw/normalized reads
 * themselves don't depend on funannotate version, only what's later done
 * with them does, so reusing one cell's fetch here is safe.
 */

params.target = "${launchDir}/output"
params.rnaseq_reads_dir = null

include { INPUT_CHECK } from './subworkflows/local/input_check'
include { GENEMARK_RUN } from './modules/local/genemark_run'

// Independent RNA-seq-to-genome spliced alignment for GeneMark --ET intron
// hints -- plain hisat2 + samtools, deliberately NOT routed through any
// cell's funannotate train/PASA (see header comment). Host modules (not a
// container) since this has nothing to do with which funannotate
// version/environment is under test; self-contained beforeScript rather than
// a shared provisioning label since this process exists only here.
process ALIGN_RNASEQ_HINTS {
    tag "${species_tag}"
    label 'process_medium'

    input:
    tuple val(species_tag), val(genome_fa), path(r1), path(r2), path(se)

    output:
    tuple val(species_tag), path("${species_tag}.training.bam"), emit: bam

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load hisat2 samtools

    GENOME_GZ="${genome_fa}"
    case "\$GENOME_GZ" in
        *.gz) gzip -dc "\$GENOME_GZ" > genome.fa ;;
        *)    cp "\$GENOME_GZ" genome.fa ;;
    esac
    hisat2-build -p ${task.cpus} genome.fa idx >hisat2-build.log 2>&1

    READ_ARGS=()
    if [ -s "${r1}" ] && [ -s "${r2}" ]; then
        READ_ARGS+=(-1 "${r1}" -2 "${r2}")
    fi
    if [ -s "${se}" ]; then
        READ_ARGS+=(-U "${se}")
    fi
    if [ "\${#READ_ARGS[@]}" -eq 0 ]; then
        echo "ERROR: no non-empty reads for ${species_tag} (r1=${r1} r2=${r2} se=${se})" >&2
        exit 1
    fi

    hisat2 -p ${task.cpus} --dta -x idx "\${READ_ARGS[@]}" \\
        | samtools sort -@ ${task.cpus} -o "${species_tag}.training.bam" -
    rm -f genome.fa idx.*.ht2
    """

    stub:
    """
    touch "${species_tag}.training.bam"
    """
}

// ── Publish GENEMARK_RUN's per-genome outputs into params.target ────────────
// GENEMARK_RUN itself has no publishDir (it's normally consumed directly out
// of the Nextflow work dir by TRAIN_PREDICT in the same run) — this sidecar
// needs the GTF to persist at a stable path outside work/ so other pipeline
// invocations (other cells, other runs) can find it via --genemark_sidecar_dir.
process PUBLISH_GENEMARK {
    tag "${meta.id}"
    label 'process_single'
    publishDir params.target, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(gtf)

    output:
    path(gtf)

    script:
    """
    # no-op: publishDir does the copy; this process just exists to attach one.
    """

    stub:
    """
    """
}

// .mod is optional (GENEMARK_RUN.out.mod, itself declared optional: true —
// too-small/fragmented genomes skip GeneMark entirely and never emit one), so
// it gets its own publish call over only the rows that have one, rather than
// forcing a path qualifier through a channel that can carry nulls.
process PUBLISH_GENEMARK_MOD {
    tag "${meta.id}"
    label 'process_single'
    publishDir params.target, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(mod)

    output:
    path(mod)

    script:
    """
    # no-op: publishDir does the copy; this process just exists to attach one.
    """

    stub:
    """
    """
}

workflow {
    // params.genemark_container_mode / params.container_genemark come from
    // the singularity provisioning axis (conf/provision_singularity.config);
    // force container mode regardless of which provisioning profile is
    // loaded, so the sidecar result never depends on host-module GeneMark
    // being installed/licensed at the site running it.
    params.genemark_container_mode = true

    INPUT_CHECK()

    // ── Branch by RNA-seq availability (per species, not per genome) ────────
    def readsDir = params.rnaseq_reads_dir
    // species_tag convention matches funannotate.nf/train_predict.nf/
    // ani_reuse.nf/annotate_genome.nf: meta.species with whitespace -> '_'
    // (this is also FETCH_RNASEQ:SRA_FETCH's species_tag, so it's what the
    // <species_tag>_norm_R{1,2,SE}.fastq.gz filenames under
    // --rnaseq_reads_dir actually use).
    def genomes_with_species = INPUT_CHECK.out.genomes.map { meta, genome_fa ->
        tuple(meta.species.replaceAll(/\s+/, '_') as String, meta, genome_fa.toString())
    }

    def has_reads
    if (readsDir) {
        has_reads = genomes_with_species.map { species_tag, meta, genome_fa ->
            def r1 = file("${readsDir}/${species_tag}_norm_R1.fastq.gz")
            def r2 = file("${readsDir}/${species_tag}_norm_R2.fastq.gz")
            def se = file("${readsDir}/${species_tag}_norm_SE.fastq.gz")
            def paired = r1.exists() && r1.size() > 0 && r2.exists() && r2.size() > 0
            def single = se.exists() && se.size() > 0
            tuple(species_tag, meta, genome_fa, paired || single, r1, r2, se)
        }
    } else {
        has_reads = genomes_with_species.map { species_tag, meta, genome_fa ->
            tuple(species_tag, meta, genome_fa, false, file('NO_FILE'), file('NO_FILE'), file('NO_FILE'))
        }
    }

    def branched = has_reads.branch {
        with_reads: it[3]
        no_reads:   true
    }

    // One alignment per SPECIES (not per genome/strain) -- .unique() on
    // species_tag, since GeneMark's genome argument is per-genome but the
    // intron hints only need one representative alignment per species'
    // read set. Genomes sharing a species reuse the same training BAM.
    def align_input = branched.with_reads
        .map { species_tag, meta, genome_fa, _has, r1, r2, se -> tuple(species_tag, genome_fa, r1, r2, se) }
        .unique { it[0] }

    ALIGN_RNASEQ_HINTS(align_input)

    def et_input = branched.with_reads
        .map { species_tag, meta, genome_fa, _has, r1, r2, se -> tuple(species_tag, meta, genome_fa) }
        .combine(ALIGN_RNASEQ_HINTS.out.bam, by: 0)
        .map { species_tag, meta, genome_fa, bam ->
            tuple(meta, genome_fa, 'ET', bam.toString(), true, '')
        }

    def es_input = branched.no_reads
        .map { species_tag, meta, genome_fa, _has, _r1, _r2, _se ->
            // mode='ES' (no RNA-seq hints), no shared_mod, force_independent=true:
            // always a fresh, independent, genome-only GeneMark-ES run.
            tuple(meta, genome_fa, 'ES', '', true, '')
        }

    GENEMARK_RUN(et_input.mix(es_input))

    PUBLISH_GENEMARK(GENEMARK_RUN.out.gtf)
    PUBLISH_GENEMARK_MOD(GENEMARK_RUN.out.mod)
}
