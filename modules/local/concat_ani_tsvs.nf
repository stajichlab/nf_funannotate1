/*
 * CONCAT_ANI_TSVS — merge per-group ANI TSVs into one table + asmid manifest.
 *
 * Port of BFD/Fungi_BFD/nextflow/modules/ani/report/CONCAT_ANI_TSVS/main.nf.
 * Each input TSV: query<TAB>ref<TAB>ANI[<TAB>AF_ref<TAB>AF_query]; the merged
 * output keeps the 3-column format with a header row, so
 * PICK_REPRESENTATIVE_STRAIN's reader sees query/ref/ANI.
 */
process CONCAT_ANI_TSVS {
    tag   "CONCAT_ANI_TSVS"
    label 'report'
    publishDir "${params.ani_outdir}/${params.ani_method}/${params.compare}", mode: 'copy', overwrite: true

    input:
    path manifest    // file listing .ani.tsv paths, one per line
    path asmid_file  // one ASMID per line, already deduped/sorted by the caller

    output:
    path("all_pairs_merged.tsv"),           emit: out
    path("all_pairs_merged.asmid_manifest.txt"), emit: asmid_manifest

    script:
    """
    printf 'query\\tref\\tANI\\n' > all_pairs_merged.tsv
    xargs -a "${manifest}" -I{} sh -c '[ -s "{}" ] && cat "{}"' >> all_pairs_merged.tsv

    cp "${asmid_file}" all_pairs_merged.asmid_manifest.txt
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > all_pairs_merged.tsv
    touch all_pairs_merged.asmid_manifest.txt
    """
}
