// FCS-GX contaminant purge + length-filter for a single genome assembly.
// Uses storeDir so the task is skipped if the cleaned .fa.gz already exists.
// Runs under the 'genome_clean' provisioning label (AAFTF + FCS-GX). The phylum
// taxid comes from TAXON_PHYLUM (meta.phylum); taxonkit is no longer in this image;
// resources are overridden per-profile by withName: '.*:GENOME_CLEAN' in
// conf/profile_annotate.config (500 GB highmem with FCS, 8 GB short without).
process GENOME_CLEAN {
    label 'genome_clean'
    label 'process_high'
    tag "${meta.id}"

    storeDir "${launchDir}/input_clean_genomes"

    input:
    tuple val(meta), path(genome_gz), val(taxondb)   // taxondb is unused now; kept so the call site is unchanged

    output:
    tuple val(meta), path("${meta.asmid}.fa.gz"), emit: genome

    script:
    def out      = meta.id
    def asmid    = meta.asmid
    def phylum   = meta.phylum
    """
    if [ ! -f "${genome_gz}" ]; then
        echo "ERROR: genome_gz not found at path: ${genome_gz}" >&2
        exit 1
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    SCRATCH=\$(printf '%s' "\${SCRATCH:-.}" | tr -d '\\n\\r')
    # An inherited \$SCRATCH may point at a path that doesn't exist / isn't
    # writable on this node — fall back to the task workdir.
    if [ ! -d "\$SCRATCH" ] || [ ! -w "\$SCRATCH" ]; then
        SCRATCH="."
    fi

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
        echo "[INFO] Phylum for ${asmid} (from TAXON_PHYLUM): ${phylum}"
        echo "[INFO] FCS-GX purge + cleaning genome for ${asmid}..."
        # AAFTF >= 0.7: -o is an output DIRECTORY. It gets <input prefix>.fcsgx_cleaned.fasta and
        # <input prefix>.fcs_gx_report.txt; the input prefix here is ${asmid}.raw.
        FCS_OUT=\$SCRATCH/${asmid}.fcs_out
        mkdir -p \$FCS_OUT
        PURGED=\$FCS_OUT/${asmid}.raw.fcsgx_cleaned.fasta
        # AAFTF releases before 0.7 (e.g. the AAFTF/0.6.0 module) take -o <file>; 0.7 takes -o <dir>.
        # The help goes to a variable first so grep cannot close the pipe early.
        FCS_HELP=\$(AAFTF fcs_gx_purge --help 2>&1 || true)
        if [[ "\$FCS_HELP" == *"--outdir"* ]]; then FCS_O=\$FCS_OUT; else FCS_O=\$PURGED; fi
        AAFTF fcs_gx_purge -d /dev/shm/gxdb/all \
            -i \$SCRATCH/${asmid}.raw.fa -c ${task.cpus} \
            -o \$FCS_O \
            -t "${phylum}" -w \$SCRATCH/${asmid}.fcs_report
        # fcs_gx_purge's exit code is not a reliable success signal; require the cleaned file.
        if [ ! -s "\$PURGED" ]; then
            echo "ERROR: fcs_gx_purge produced no \$PURGED for ${asmid}" >&2
            exit 1
        fi
        mkdir -p ${launchDir}/input_clean_genomes/clean
        ${params.clean_script} --len ${params.min_contig_len} < \$PURGED > ${asmid}.fa
        if [ ! -s ${asmid}.fa ]; then
            echo "ERROR: no contigs left in ${asmid}.fa after the FCS-GX purge and length filter" >&2
            exit 1
        fi
        pigz -c \$PURGED > \$SCRATCH/${asmid}.purge.fasta.gz
        [ -f \$FCS_OUT/${asmid}.raw.fcs_gx_report.txt ] && pigz -c \$FCS_OUT/${asmid}.raw.fcs_gx_report.txt > \$SCRATCH/${asmid}.fcs_gx_report.txt.gz
        mv \$SCRATCH/${asmid}.purge.fasta.gz ${launchDir}/input_clean_genomes/clean/
        [ -f \$SCRATCH/${asmid}.fcs_gx_report.txt.gz ] && \
            mv \$SCRATCH/${asmid}.fcs_gx_report.txt.gz ${launchDir}/input_clean_genomes/clean/
        rm -rf \$FCS_OUT
    fi
    # Normalize FASTA deflines to the accession (first whitespace token).
    # funannotate rejects headers >24 chars and, more importantly, requires ONE
    # naming scheme across the genemark GTF, repeat-mask output and predict
    # input; NCBI deflines (long descriptions) survive clean_genome_fa.py
    # verbatim and produced inconsistent seqids in EVM. Short accessions match
    # what the classic `funannotate clean` step produced.
    awk '/^>/{print \$1; next} {print}' ${asmid}.fa > ${asmid}.fa.hdr && mv ${asmid}.fa.hdr ${asmid}.fa
    # Deliver the clean genome gzip-compressed to save space in input_clean_genomes;
    # downstream tools inflate it on the fly (they cannot read a gzipped FASTA via -i).
    pigz -f ${asmid}.fa
    echo "[INFO] Clean genome written: ${asmid}.fa.gz (\$(du -sh ${asmid}.fa.gz | cut -f1))"
    rm -f \$SCRATCH/${asmid}.raw.fa
    ${projectDir}/bin/record_versions.sh "${params.target}/${meta.id}/versions" "${task.process}" "${meta.id}" "${task.container ?: 'native'}" "${workflow.manifest.version ?: ''}@${workflow.commitId ?: workflow.revision ?: 'local'}" \\
        'AAFTF=AAFTF --version 2>&1 | grep -E '\\''^AAFTF [0-9]'\\'' | head -1' \\
        'fcs_gx=command -v gx >/dev/null && ls \$(dirname \$(command -v gx))/../conda-meta 2>/dev/null | grep ncbi-fcs-gx | head -1' \\
        'pigz=pigz --version'
    """

    stub:
    def asmid = meta.asmid
    """
    echo ">stub_${asmid}" | pigz -c > ${asmid}.fa.gz
    mkdir -p ${launchDir}/input_clean_genomes/clean
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fasta
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.fcs_gx_report.txt
    """
}
