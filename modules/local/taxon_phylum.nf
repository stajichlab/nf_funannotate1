// Resolve the numeric NCBI taxid of a genome's phylum (kingdom as fallback)
// with taxonkit. AAFTF fcs_gx_purge -t needs that taxid.
// taxonkit runs in its own container (params.container_taxonkit) because the
// AAFTF >= 0.7 image no longer ships it.
process TAXON_PHYLUM {
    label 'taxonkit'
    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), val(genome), val(taxondb)

    output:
    tuple val(meta), val(genome), path('phylum.txt'), emit: phylum

    script:
    def taxonid = meta.taxonid
    """
    set -euo pipefail
    TAXONKIT_DB=${taxondb}
    phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1 || true)
    if [ -z "\$phylum" ]; then
        # name2taxid returns 2 lines when the input is Fungi; uniq/head keep one.
        phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1 || true)
    fi
    if ! [[ "\$phylum" =~ ^[0-9]+\$ ]]; then
        echo "ERROR: no numeric phylum taxid for ${meta.asmid} (taxonid=${taxonid}); got '\$phylum'." >&2
        echo "The phylum or kingdom name may be missing from names.dmp in TAXONKIT_DB=\$TAXONKIT_DB (outdated taxdump)." >&2
        exit 1
    fi
    echo "[INFO] Phylum for ${meta.asmid} (taxonid=${taxonid}): \$phylum" >&2
    echo "\$phylum" > phylum.txt
    ${projectDir}/bin/record_versions.sh "${params.target}/${meta.id}/versions" "${task.process}" "${meta.id}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'taxonkit=taxonkit version' >&2
    """

    stub:
    """
    echo 4890 > phylum.txt
    """
}
