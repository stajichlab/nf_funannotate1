// Write zero-byte paired FASTQ placeholder files for species with no SRA data.
// Runs in the driver process (no SLURM job) for the no-data branch from SRA branching.
// storeDir caches results so placeholders are not re-created on resume.
process WRITE_EMPTY_READS {
    label 'setup'
    label 'process_single'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads"

    input:
    val(species_tag)

    output:
    tuple val(species_tag), path("${species_tag}_norm_R1.fastq.gz"), path("${species_tag}_norm_R2.fastq.gz"),
          path("${species_tag}_norm_SE.fastq.gz"), emit: reads

    script:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    echo "[INFO] No SRA data for ${species_tag}; created empty read placeholders"
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    """
}
