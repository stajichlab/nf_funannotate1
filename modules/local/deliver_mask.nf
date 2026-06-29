// Copy a finished soft-masked genome into params.masked_dir (input_clean_genomes/),
// OVERWRITING any pre-existing tantan mask. Emits tuple(val(asmid), path) so the
// downstream MASK_GENOME subworkflow can reconstruct the final path.
// Resources: label 'deliver' in conf/profile_earlgrey.config (trivial I/O).
process DELIVER_MASK {
    label 'deliver'
    label 'process_single'
    tag  "${asmid}"

    publishDir path: { workflow.stubRun ? "${launchDir}/work/stub_masked" : params.masked_dir },
               mode: 'copy', overwrite: true

    input:
    tuple val(asmid), path(masked)

    output:
    tuple val(asmid), path("${asmid}.masked.fasta.gz"), emit: delivered

    script:
    'true'

    stub:
    'true'
}
