// Batched FCS-GX purge for a list of genome assemblies.
// Stages the ~470 GB FCS-GX database into /dev/shm ONCE (~30 min), then
// cleans every genome in the batch sequentially against that in-memory DB,
// amortizing the staging cost over up to clean_batch_size genomes.
// The phylum taxid per genome comes from TAXON_PHYLUM (meta.phylum); taxonkit
// is not in this image.
// Only used when params.skip_fcs is false; GENOME_CLEAN handles the per-genome
// fallback (clean_batch_size=0 or --skip_fcs).
//
// Resources overridden per-profile by withName: '.*:GENOME_CLEAN_BATCH' in
// conf/profile_annotate.config (500 GB highmem, 7d walltime).
process GENOME_CLEAN_BATCH {
    label 'genome_clean'
    label 'process_long'
    tag "clean_batch_${task.index}"

    input:
    tuple val(items), val(taxondb)

    output:
    path "clean_batch_*.manifest.tsv", emit: manifest

    script:
    def batch_tsv = items.collect { row -> "${row[0].asmid}\t${row[1]}\t${row[0].taxonid}\t${row[0].phylum ?: ""}" }.join('\n')
    """
    set -uo pipefail
    source /etc/profile.d/modules.sh 2>/dev/null || true

    SCRATCH=\$(printf '%s' "\${SCRATCH:-.}" | tr -d '\\n\\r')
    # An inherited \$SCRATCH may point at a path that doesn't exist / isn't
    # writable on this node — fall back to the task workdir.
    if [ ! -d "\$SCRATCH" ] || [ ! -w "\$SCRATCH" ]; then
        SCRATCH="."
    fi
    DEST=${launchDir}/input_clean_genomes
    mkdir -p \$DEST/clean

    MANIFEST=clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST

    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF

    n_total=\$(grep -c . batch.tsv || true)
    echo "[INFO] batch ${task.index}: \$n_total genomes to consider"

    # Stage the FCS-GX DB into /dev/shm ONCE for the whole batch (~30 min). FCS_GX_KEEP_SHM=1
    # tells the staging script not to register its own per-shell EXIT cleanup; we remove the
    # RAM copy ourselves when the batch finishes (or aborts) via the trap below.
    export FCS_GX_KEEP_SHM=1
    source ${params.fcs_shm_script}
    trap 'rm -rf /dev/shm/gxdb 2>/dev/null || true' EXIT
    if [ ! -f /dev/shm/gxdb/all.gxi ]; then
        echo "[ERROR] FCS-GX DB not staged into /dev/shm/gxdb; aborting batch" >&2
        exit 1
    fi

    i=0
    fcs_fails=0
    while IFS=\$'\\t' read -r asmid gz taxonid phylum; do
        [ -z "\$asmid" ] && continue
        i=\$((i+1))
        target=\$DEST/\${asmid}.fa.gz
        if [ -s "\$target" ]; then
            echo "[\$i/\$n_total][SKIP] \$asmid already cleaned"
            printf '%s\\t%s\\n' "\$asmid" "\$target" >> \$MANIFEST
            continue
        elif [ -s "\$DEST/\${asmid}.fa" ]; then
            # Back-compat: a prior run may have left an uncompressed .fa.
            echo "[\$i/\$n_total][SKIP] \$asmid already cleaned (uncompressed)"
            printf '%s\\t%s\\n' "\$asmid" "\$DEST/\${asmid}.fa" >> \$MANIFEST
            continue
        fi
        if [ ! -f "\$gz" ]; then
            echo "[\$i/\$n_total][WARN] missing genome for \$asmid: \$gz" >&2
            continue
        fi

        echo "[\$i/\$n_total][INFO] \$asmid taxonid=\$taxonid phylum=\$phylum"

        # Accept gzipped (NCBI_ASM .fna.gz) or plain (local GENOME column) FASTA input.
        if printf '%s' "\$gz" | grep -qiE '\\.gz\$'; then
            pigz -dc "\$gz" > \$SCRATCH/\${asmid}.raw.fa
        else
            cat "\$gz" > \$SCRATCH/\${asmid}.raw.fa
        fi
        # AAFTF >= 0.7: fcs_gx_purge -o is an output DIRECTORY. It gets
        # <input prefix>.fcsgx_cleaned.fasta and <input prefix>.fcs_gx_report.txt; the input
        # prefix here is <asmid>.raw. Its exit code is NOT a reliable success signal (it can
        # return 0 although run_gx.py produced no report; a stale or corrupt /dev/shm/gxdb did
        # this on 2026-09-26), so success is gated on the cleaned file existing and being non-empty.
        FCS_OUT=\$SCRATCH/\${asmid}.fcs_out
        PURGED=\$FCS_OUT/\${asmid}.raw.fcsgx_cleaned.fasta
        mkdir -p \$FCS_OUT
        if AAFTF fcs_gx_purge -d /dev/shm/gxdb/all \\
            -i \$SCRATCH/\${asmid}.raw.fa -c ${task.cpus} \\
            -o \$FCS_OUT \\
            -t "\$phylum" -w \$SCRATCH/\${asmid}.fcs_report \\
            && [ -s "\$PURGED" ] ; then
            ${params.clean_script} --len ${params.min_contig_len} < \$PURGED > \$SCRATCH/\${asmid}.clean.fa \\
                && pigz -c \$SCRATCH/\${asmid}.clean.fa > \${target}.tmp \\
                && mv \${target}.tmp \$target
            rm -f \$SCRATCH/\${asmid}.clean.fa
            if [ -s \$target ] && [ "\$(pigz -dc \$target | grep -c '^>')" -gt 0 ]; then
                echo "[\$i/\$n_total][OK] \$asmid -> \$target (\$(du -sh \$target | cut -f1))"
                pigz -c \$PURGED > \$DEST/clean/\${asmid}.purge.fasta.gz 2>/dev/null || true
                [ -f \$FCS_OUT/\${asmid}.raw.fcs_gx_report.txt ] && pigz -c \$FCS_OUT/\${asmid}.raw.fcs_gx_report.txt > \$DEST/clean/\${asmid}.fcs_gx_report.txt.gz
                printf '%s\\t%s\\n' "\$asmid" "\$target" >> \$MANIFEST
            else
                echo "[\$i/\$n_total][FAIL] \$asmid: \$target was not produced or has no contigs (clean/compress/move step failed)" >&2
                rm -f \$target
                fcs_fails=\$((fcs_fails+1))
            fi
        else
            echo "[\$i/\$n_total][FAIL] fcs_gx_purge failed for \$asmid (no report / empty cleaned fasta)" >&2
            fcs_fails=\$((fcs_fails+1))
        fi
        rm -rf \$FCS_OUT
        rm -f \$SCRATCH/\${asmid}.raw.fa
    done < batch.tsv

    echoed_manifest=\$(grep -c . \$MANIFEST || echo 0)
    if [ "\${fcs_fails}" -gt 0 ] || [ "\${echoed_manifest}" -eq 0 ]; then
        echo "[ERROR] batch ${task.index}: \${fcs_fails} purge failure(s); manifest has \${echoed_manifest} entries" >&2
        exit 1
    fi
    echo "[INFO] batch ${task.index} complete: \${echoed_manifest} cleaned genomes in manifest"
    """

    stub:
    def batch_tsv = items.collect { row -> "${row[0].asmid}\t${row[1]}\t${row[0].taxonid}\t${row[0].phylum ?: ""}" }.join('\n')
    """
    DEST=${launchDir}/input_clean_genomes
    mkdir -p \$DEST/clean
    MANIFEST=clean_batch_${task.index}.manifest.tsv
    : > \$MANIFEST
    cat > batch.tsv <<'BATCH_EOF'
${batch_tsv}
BATCH_EOF
    while IFS=\$'\\t' read -r asmid gz taxonid phylum; do
        [ -z "\$asmid" ] && continue
        echo ">stub_\${asmid}" | pigz -c > \$DEST/\${asmid}.fa.gz
        touch \$DEST/clean/\${asmid}.purge.fasta \$DEST/clean/\${asmid}.fcs_gx_report.txt
        printf '%s\\t%s\\n' "\$asmid" "\$DEST/\${asmid}.fa.gz" >> \$MANIFEST
    done < batch.tsv
    """
}
