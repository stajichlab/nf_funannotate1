// FCS-GX contaminant purge + length-filter for a single genome assembly.
// Uses storeDir so the task is skipped if the cleaned .fa.gz already exists.
// Runs under the 'genome_clean' provisioning label (AAFTF + taxonkit + FCS-GX);
// resources are overridden per-profile by withName: '.*:GENOME_CLEAN' in
// conf/profile_annotate.config (500 GB highmem with FCS, 8 GB short without).
process GENOME_CLEAN {
    label 'genome_clean'
    label 'process_high'
    tag "${meta.id}"

    storeDir "${launchDir}/input_clean_genomes"

    input:
    tuple val(meta), path(genome_gz), val(taxondb)

    output:
    tuple val(meta), path("${meta.asmid}.fa.gz"), emit: genome

    script:
    def out      = meta.id
    def asmid    = meta.asmid
    def taxonid  = meta.taxonid
    """
    if [ ! -f "${genome_gz}" ]; then
        echo "ERROR: genome_gz not found at path: ${genome_gz}" >&2
        exit 1
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    SCRATCH=\$(printf '%s' "\${SCRATCH:-.}" | tr -d '\\n\\r')

    echo "[INFO] Decompressing genome for ${asmid}..."
    # Accept either a gzipped (NCBI_ASM .fna.gz) or plain (local GENOME column) FASTA.
    if printf '%s' "${genome_gz}" | grep -qiE '\\.gz\$'; then
        pigz -dc ${genome_gz} > \$SCRATCH/${asmid}.raw.fa
    else
        cat ${genome_gz} > \$SCRATCH/${asmid}.raw.fa
    fi

    if [ "${params.skip_fcs}" = "true" ]; then
        # --skip_fcs: bypass AAFTF FCS-GX contaminant purge (no 470 GB gxdb needed);
        # just length-filter the assembly.
        echo "[INFO] --skip_fcs set: skipping FCS-GX purge for ${asmid}"
        ${params.clean_script} --len ${params.min_contig_len} \
            -i \$SCRATCH/${asmid}.raw.fa -o ${asmid}.fa
    else
        # Ensure /dev/shm/gxdb is present on this node; register for cleanup when done.
        source ${params.fcs_shm_script}
        TAXONKIT_DB=${taxondb}
        phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
        if [ -z "\$phylum" ]; then
            phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1)
            # weird we are getting 2 lines from name2taxid when input is Fungi add the uniq/head -n 1 to ensure only one line
        fi
        echo "[INFO] Phylum for ${asmid} (taxonid=${taxonid}): \$phylum"
        echo "[INFO] FCS-GX purge + cleaning genome for ${asmid}..."
        AAFTF fcs_gx_purge --db /dev/shm/gxdb/all \
            -i \$SCRATCH/${asmid}.raw.fa --cpus ${task.cpus} \
            -o \$SCRATCH/${asmid}.purge.fasta \
            -t "\$phylum" -w \$SCRATCH/${asmid}.fcs_report
        mkdir -p ${launchDir}/input_clean_genomes/clean
        cat \$SCRATCH/${asmid}.purge.fasta | \
            ${params.clean_script} --len ${params.min_contig_len} > ${asmid}.fa
        pigz \$SCRATCH/${asmid}.purge.fasta
        [ -f \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv ] && pigz \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv
        mv \$SCRATCH/${asmid}.purge.fasta.gz ${launchDir}/input_clean_genomes/clean/
        [ -f \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv.gz ] && \
            mv \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv.gz ${launchDir}/input_clean_genomes/clean/
    fi
    # Deliver the clean genome gzip-compressed to save space in input_clean_genomes;
    # downstream tools inflate it on the fly (they cannot read a gzipped FASTA via -i).
    pigz -f ${asmid}.fa
    echo "[INFO] Clean genome written: ${asmid}.fa.gz (\$(du -sh ${asmid}.fa.gz | cut -f1))"
    rm -f \$SCRATCH/${asmid}.raw.fa
    """

    stub:
    def asmid = meta.asmid
    """
    echo ">stub_${asmid}" | pigz -c > ${asmid}.fa.gz
    mkdir -p ${launchDir}/input_clean_genomes/clean
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fasta
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fcs_gx-taxonomy.tsv
    """
}
