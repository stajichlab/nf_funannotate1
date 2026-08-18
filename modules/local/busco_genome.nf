/*
 * BUSCO_GENOME — per-assembly BUSCO completeness score (genome mode).
 *
 * Containerized port of BFD/Fungi_BFD/nextflow/modules/BFD/BUSCO_GENOME/main.nf:
 *   - same storeDir layout under params.genome_stats_outdir/BUSCO_genome and the
 *     same output naming (<asmid>.BUSCO_summary.<lineage>.txt), so downstream
 *     PICK_REPRESENTATIVE_STRAIN (which reads the raw summary txt files) and any
 *     site tooling expecting BFD paths keep working.
 *   - `module load busco` replaced by the pinned busco container (assigned via
 *     withLabel: 'busco_genome' in conf/provision_singularity.config); the host
 *     $BUSCO_LINEAGES lineage tree is bound into the container there too.
 *   - runs offline against the host lineage directory via --download_path.
 *
 * Consumes the CLEANED genome (pre-masking, entry: CLEAN_GENOMES.out.genomes),
 * mirroring BFD where busco_genome reads input_clean_genomes/.
 */
process BUSCO_GENOME {
    tag      "${meta.asmid}"
    label    'busco_genome'
    storeDir { "${params.genome_stats_outdir}/BUSCO_genome" }

    input:
    tuple val(meta), path(genome)

    output:
    path("${meta.asmid}.BUSCO_summary.${meta.busco}.txt"), emit: summary

    script:
    """
    export BUSCO_LINEAGES=${params.busco_lineages}
    SCRATCH=\$(printf '%s' "\${SCRATCH:-}" | tr -d '\\n\\r')
    SCRATCH=\${SCRATCH:-/tmp}
    # inside a container the inherited \$SCRATCH host path may not exist / be
    # writable — fall back to the task workdir in that case.
    if [ ! -d "\$SCRATCH" ] || [ ! -w "\$SCRATCH" ]; then
        SCRATCH="\$PWD"
    fi
    RUNDIR="${meta.asmid}_busco_genome"
    busco -i ${genome} \\
          -l ${meta.busco} \\
          -m genome \\
          -o \$RUNDIR \\
          --out_path \$SCRATCH \\
          -c ${task.cpus} \\
          --offline --download_path \$BUSCO_LINEAGES
    SUMMARY=\$(ls \$SCRATCH/\${RUNDIR}/short_summary.specific.*.txt 2>/dev/null | head -1)
    [ -z "\$SUMMARY" ] && SUMMARY=\$(ls \$SCRATCH/\${RUNDIR}/short_summary*.txt 2>/dev/null | head -1)
    if [ -z "\$SUMMARY" ]; then
        echo "[ERROR] BUSCO genome: no short_summary file found for ${meta.asmid}" >&2
        exit 1
    fi
    cp "\$SUMMARY" "${meta.asmid}.BUSCO_summary.${meta.busco}.txt"
    rm -rf "\$SCRATCH/\$RUNDIR"
    """

    stub:
    """
    printf '# BUSCO version: 5.x\\n# The lineage dataset is: ${meta.busco}\\n\\tC:99.0%%[S:98.0%%,D:1.0%%],F:0.5%%,M:0.5%%,n:758\\n' \\
        > "${meta.asmid}.BUSCO_summary.${meta.busco}.txt"
    """
}
