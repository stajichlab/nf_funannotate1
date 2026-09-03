#!/usr/bin/env nextflow

/*
 * genemark_sidecar — run GeneMark-ES exactly once per genome, out of band,
 * for consumption by multiple funannotate comparison cells.
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
 * axis a downstream cell uses), always fresh ES self-training (no RNA-seq
 * hints — training_bam is deliberately not wired here, so the sidecar result
 * does not depend on any cell's own RNA-seq training, keeping it identical
 * across every cell for a genome). Results are published as
 * <target>/<out>.genemark.gtf (+ .mod), matching the naming convention
 * funannotate.nf's TRAIN_PREDICT reads back via --genemark_sidecar_dir.
 *
 * Usage:
 *   nextflow run genemark_sidecar.nf -c nextflow.config \
 *       -profile genemark_sidecar,slurm,singularity \
 *       --samples samples.csv --target runs/genemark_sidecar/output -resume
 *
 * Then point every comparison cell at the result:
 *   nextflow run funannotate.nf ... --genemark_sidecar_dir runs/genemark_sidecar/output
 *
 * `--samples` only needs GENOME (already-masked FASTA) resolved per row —
 * the same samples.csv any funannotate.nf cell already uses is sufficient;
 * this script does not re-mask or re-clean, it only reads the GENOME column
 * via INPUT_CHECK.
 */

params.target = "${launchDir}/output"

include { INPUT_CHECK } from './subworkflows/local/input_check'
include { GENEMARK_RUN } from './modules/local/genemark_run'

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

    def genemark_input = INPUT_CHECK.out.genomes.map { meta, genome_fa ->
        // mode='ES' (no RNA-seq hints), no shared_mod, force_independent=true:
        // always a fresh, independent, genome-only GeneMark-ES run.
        tuple(meta, genome_fa.toString(), 'ES', '', true, '')
    }

    GENEMARK_RUN(genemark_input)

    PUBLISH_GENEMARK(GENEMARK_RUN.out.gtf)
    PUBLISH_GENEMARK_MOD(GENEMARK_RUN.out.mod)
}
