#!/usr/bin/env nextflow

/*
 * SOURCE: ../../../1KFG/common_annotate/pipeline/nextflow/funannotate.nf
 * Last synced: 2026-05-23
 * Changes vs source: removed nextflow.enable.dsl=2; params block moved to
 *                    conf/profile_annotate.config.
 *
 * Usage (from project root — a pipeline profile is REQUIRED; without it
 * params.taxondb / params.funannotate_db are null and parsing fails):
 *   sbatch nextflow/run_annotate.sh
 *   nextflow run nextflow/funannotate.nf -c nextflow/nextflow.config \
 *       -profile annotate,slurm,ucr_hpcc -resume
 */

// Data contract: every channel element is `tuple val(meta), val/path(genome)`.
// meta is a Map built by SampleUtils.makeMeta(row) — see lib/SampleUtils.groovy.
//   meta.id is the ONLY field used for tag{} and file naming.
//   meta.asmid, meta.species, meta.strain, meta.locustag, meta.busco,
//   meta.transl_table, meta.taxonid carry payload used inside process scripts.
//   header_length is NOT in meta — it comes from params.header_length (default 24).
//
// GENOME_CLEAN receives: tuple val(meta), path(genome_gz), val(taxondb)
//   → emits: tuple val(meta), path(genome_fa)   [storeDir writes input_clean_genomes/<asmid>.fa.gz]
// MASKREPEAT_TANTAN_RUN receives: tuple val(meta), val(genome_fa)
//   → emits: tuple val(meta), path(masked_fa)   [storeDir caches input_clean_genomes/<asmid>.masked.fasta.gz]
// SRA_FETCH receives: val(species_tag), val(taxonid)   [only when --run_sra_fetch]
//   → emits: val(species_tag), path(norm_R1.fastq.gz), path(norm_R2.fastq.gz), path(se)
// RNASEQ_PREPARE receives: tuple val(species_tag), val(meta), val(genome_fa), path(r1), path(r2), path(se)
//   → emits: val(species_tag), path(trinity-GG.fasta)   [storeDir caches in rnaseq_data/]
// FUNANNOTATE_TRAIN receives: tuple val(meta), val(genome_fa), path(r1), path(r2), path(se), path(trinity_fa)
//   → emits: tuple val(meta), val(genome_fa)
// FUNANNOTATE_PREDICT receives: tuple val(meta), val(genome_fa)
//   → emits metadata: tuple val(meta)

// Soft-mask each assembly using funannotate mask with tantan.
// storeDir caches the masked FASTA alongside the clean genome.
process MASKREPEAT_TANTAN_RUN {
    label 'funannotate'
    tag "${meta.id}"

    storeDir "${launchDir}/input_clean_genomes"

    cpus   8
    memory '16 GB'
    time   '2h'

    input:
    tuple val(meta), val(genome_fa)

    output:
    tuple val(meta), path("${meta.asmid}.masked.fasta.gz"), emit: masked

    script:
    def asmid = meta.asmid
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac
    funannotate mask -i "\$GENOME_IN" -o ${asmid}.masked.fasta -m tantan --cpus ${task.cpus}
    # Deliver the soft-masked genome gzip-compressed to save space; consumers inflate it.
    pigz -f ${asmid}.masked.fasta
    """

    stub:
    def asmid = meta.asmid
    """
    echo ">stub_${asmid}_masked" | pigz -c > ${asmid}.masked.fasta.gz
    """
}

// Query NCBI SRA for available paired-end RNA-seq accessions per species.
// Lightweight: runs the esearch/efetch query only — no downloading.
// Records up to 5 candidates (sorted by spot count desc) in a per-species CSV.
// storeDir caches results so re-runs skip the network query.
// To invalidate the cache for a species, delete rnaseq_reads/sra_query/<species_tag>.sra_query.csv
process SRA_QUERY {
    label 'edirect'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads/sra_query"

    cpus   1
    memory '4 GB'
    time   '30m'

    input:
    tuple val(species_tag), val(taxonid)

    output:
    tuple val(species_tag), path("${species_tag}.sra_query.csv"), emit: query_result

    script:
    """
    set -euo pipefail

    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${species_tag}.sra_query.csv

    esearch -db sra \\
        -query "txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND (BGISEQ[Platform] OR Illumina[Platform])" | \\
        efetch -format runinfo > _runinfo.tmp

    # col 1=Run, col 4=spots, col 13=LibraryStrategy, col 16=LibraryLayout, col 19=Platform
    # Prepend a platform rank (0=Illumina, 1=BGI/other) so the top 5 prefer Illumina,
    # then by spot count desc; BGI/other only fill remaining slots when Illumina runs out.
    awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ && \$4+0>=250000 {rank=(\$19~/[Ii]llumina/)?0:1; printf "%d,%s,%s,%s\\n", rank, \$1, \$4, \$19}' _runinfo.tmp | \\
        sort -t',' -k1,1n -k3,3rn | \\
        head -n 5 | \\
        while IFS=',' read -r rank acc spots platform; do
            printf '%s,%s,%s,%s,%s,PAIRED\\n' "${species_tag}" "${taxonid}" "\$acc" "\$spots" "\$platform"
        done >> ${species_tag}.sra_query.csv

    rm -f _runinfo.tmp
    NHITS=\$(awk 'END{print NR-1}' ${species_tag}.sra_query.csv)
    echo "[INFO] Found \$NHITS paired-end SRA accessions for ${species_tag} (taxonid=${taxonid})"

    # SE fallback: if no PE hits found and enable_single_end is true, query SINGLE layout
    if [ "${params.enable_single_end}" = "true" ] && [ "\$NHITS" -eq 0 ]; then
        esearch -db sra \\
            -query "txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND SINGLE[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND Illumina[Platform]" | \\
            efetch -format runinfo > _runinfo_se.tmp
        awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="SINGLE" && \$1~/^[SDE]RR/ && \$4+0>=250000 {printf "%s,%s,%s\\n", \$1, \$4, \$19}' _runinfo_se.tmp | \\
            sort -t',' -k2 -rn | \\
            head -n ${params.max_rnaseq_se_runs} | \\
            while IFS=',' read -r acc spots platform; do
                printf '%s,%s,%s,%s,%s,SINGLE\\n' "${species_tag}" "${taxonid}" "\$acc" "\$spots" "\$platform"
            done >> ${species_tag}.sra_query.csv
        rm -f _runinfo_se.tmp
        NHITS=\$(awk 'END{print NR-1}' ${species_tag}.sra_query.csv)
        echo "[INFO] SE fallback: found \$NHITS single-end accessions for ${species_tag}"
    fi
    """

    stub:
    """
    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${species_tag}.sra_query.csv
    printf '%s,%s,SRR000001,1000000,ILLUMINA,PAIRED\n' "${species_tag}" "${taxonid}" >> ${species_tag}.sra_query.csv
    echo "[STUB] SRA_QUERY for ${species_tag}"
    """
}

// Batched SRA query: handles params.sra_query_batch_size species per SLURM job.
// maxForks 4 caps concurrent jobs to avoid overwhelming NCBI.
// Per-species esearch/efetch is retried up to 3 times inline with exponential
// backoff before writing an empty CSV.  Existing per-species CSVs already
// present in rnaseq_reads/sra_query/ are reused without re-querying.
process SRA_QUERY_BATCH {
    label 'edirect'
    tag "${species_tags[0]}_+${species_tags.size() - 1}_more"

    publishDir "${launchDir}/rnaseq_reads/sra_query", mode: 'copy', overwrite: false

    maxForks 4
    cpus   1
    memory '4 GB'
    time   '4h'

    input:
    tuple val(species_tags), val(taxonids)

    output:
    path("*.sra_query.csv"), emit: query_results

    script:
    def cache_dir  = "${launchDir}/rnaseq_reads/sra_query"
    def batch_args = [species_tags, taxonids].transpose()
                         .collect { st, tid -> "${st}\\t${tid}" }
                         .join('\\n')
    """
    set -uo pipefail

    printf '${batch_args}\\n' > batch_input.tsv

    query_species() {
        local stag="\$1" tid="\$2" attempt

        for attempt in 1 2 3; do
            rm -f "_runinfo_\${stag}.tmp"
            if timeout 120 bash -c \\
                    "esearch -db sra -query 'txid\${tid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND (BGISEQ[Platform] OR Illumina[Platform])' | efetch -format runinfo" \\
                    < /dev/null > "_runinfo_\${stag}.tmp"; then
                return 0
            fi
            echo "[WARN] Attempt \${attempt}/3 failed or timed out for \${stag}"
            [ "\${attempt}" -lt 3 ] && sleep \$((attempt * 30))
        done
        rm -f "_runinfo_\${stag}.tmp"
        return 1
    }

    while IFS=\$(printf '\\t') read -r species_tag taxonid; do
        cached="${cache_dir}/\${species_tag}.sra_query.csv"
        if [ -s "\$cached" ]; then
            cp "\$cached" "\${species_tag}.sra_query.csv"
            echo "[INFO] Reusing cached result for \${species_tag}"
            continue
        fi

        if query_species "\${species_tag}" "\${taxonid}"; then
            printf 'species_tag,taxonid,sra_accession,spots,platform,layout\\n' > "\${species_tag}.sra_query.csv"
            # col 1=Run, col 4=spots, col 13=LibraryStrategy, col 16=LibraryLayout, col 19=Platform
            # Prepend a platform rank (0=Illumina, 1=BGI/other) so the top 5 prefer Illumina,
            # then by spot count desc; BGI/other only fill remaining slots when Illumina runs out.
            awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ && \$4+0>=250000 {rank=(\$19~/[Ii]llumina/)?0:1; printf "%d,%s,%s,%s\\n", rank, \$1, \$4, \$19}' "_runinfo_\${species_tag}.tmp" | \\
                sort -t',' -k1,1n -k3,3rn | \\
                head -n 5 | \\
                while IFS=',' read -r rank acc spots platform; do
                    printf '%s,%s,%s,%s,%s,PAIRED\\n' "\${species_tag}" "\${taxonid}" "\$acc" "\$spots" "\$platform"
                done >> "\${species_tag}.sra_query.csv"
            rm -f "_runinfo_\${species_tag}.tmp"
            NHITS=\$(awk 'END{print NR-1}' "\${species_tag}.sra_query.csv")
            echo "[INFO] Found \$NHITS paired-end accessions for \${species_tag} (taxonid=\${taxonid})"
            # SE fallback: if no PE hits and enable_single_end, query SINGLE layout
            if [ "${params.enable_single_end}" = "true" ] && [ "\$NHITS" -eq 0 ]; then
                rm -f "_runinfo_se_\${species_tag}.tmp"
                if timeout 120 bash -c \\
                        "esearch -db sra -query 'txid\${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND SINGLE[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND Illumina[Platform]' | efetch -format runinfo" \\
                        < /dev/null > "_runinfo_se_\${species_tag}.tmp"; then
                    awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="SINGLE" && \$1~/^[SDE]RR/ && \$4+0>=250000 {printf "%s,%s,%s\\n", \$1, \$4, \$19}' "_runinfo_se_\${species_tag}.tmp" | \\
                        sort -t',' -k2 -rn | \\
                        head -n ${params.max_rnaseq_se_runs} | \\
                        while IFS=',' read -r acc spots platform; do
                            printf '%s,%s,%s,%s,%s,SINGLE\\n' "\${species_tag}" "\${taxonid}" "\$acc" "\$spots" "\$platform"
                        done >> "\${species_tag}.sra_query.csv"
                fi
                rm -f "_runinfo_se_\${species_tag}.tmp"
                NHITS=\$(awk 'END{print NR-1}' "\${species_tag}.sra_query.csv")
                echo "[INFO] SE fallback: \$NHITS single-end accessions for \${species_tag}"
            fi
        else
            printf 'species_tag,taxonid,sra_accession,spots,platform,layout\\n' > "\${species_tag}.sra_query.csv"
            echo "[WARN] All 3 attempts failed for \${species_tag}; writing empty CSV"
        fi
    done < batch_input.tsv
    """

    stub:
    def stub_args = [species_tags, taxonids].transpose()
                        .collect { st, tid -> "${st}\\t${tid}" }
                        .join('\\n')
    """
    printf '${stub_args}\\n' > batch_input.tsv
    while IFS=\$(printf '\\t') read -r species_tag taxonid; do
        printf 'species_tag,taxonid,sra_accession,spots,platform,layout\\n' > "\${species_tag}.sra_query.csv"
        printf '%s,%s,SRR000001,1000000,ILLUMINA,PAIRED\\n' "\${species_tag}" "\${taxonid}" >> "\${species_tag}.sra_query.csv"
    done < batch_input.tsv
    echo "[STUB] SRA_QUERY_BATCH (${species_tags.size()} species)"
    """
}

// Merge all per-species SRA query CSVs into a single named manifest.
// Output: {stem}.rnaseq_sra.csv written alongside the input samples file.
// Columns: species_tag, taxonid, sra_accession, spots, platform
process COLLECT_SRA_QUERY {
    label 'setup'
    publishDir { file(params.samples).parent.toAbsolutePath().toString() }, mode: 'copy'

    cpus   1
    memory '1 GB'
    time   '10m'

    input:
    path(query_csvs)
    val(stem)

    output:
    path("${stem}.rnaseq_sra.csv"), emit: manifest

    script:
    """
    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${stem}.rnaseq_sra.csv
    for f in ${query_csvs}; do
        tail -n +2 "\$f" >> ${stem}.rnaseq_sra.csv
    done
    NSPECIES=\$(awk -F',' 'NR>1{print \$1}' ${stem}.rnaseq_sra.csv | sort -u | wc -l)
    NACCESSIONS=\$(awk 'NR>1' ${stem}.rnaseq_sra.csv | wc -l)
    echo "[INFO] ${stem}.rnaseq_sra.csv: \$NACCESSIONS accessions across \$NSPECIES species with RNA-seq data"
    """

    stub:
    """
    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${stem}.rnaseq_sra.csv
    """
}

// Write zero-byte paired FASTQ placeholder files for species with no SRA data.
// Called only for species whose SRA_QUERY CSV has no data rows, avoiding a
// SLURM job allocation for what would be an immediate empty-file write.
process WRITE_EMPTY_READS {
    label 'setup'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads"

    cpus   1
    memory '1 GB'
    time   '5m'

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

// Download and normalize up to params.max_rnaseq_runs SRA accessions for species
// that have RNA-seq data. Accessions are read from the pre-queried per-species CSV
// produced by SRA_QUERY, so no NCBI network call is made here.
// Only invoked for species with data rows in their SRA_QUERY CSV; WRITE_EMPTY_READS
// handles the no-data case at the channel level.
process SRA_FETCH {
    label 'sra'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads"

    cpus   32
    memory '96 GB'
    time   '2h'

    input:
    tuple val(species_tag), path(sra_query_csv)

    output:
    tuple val(species_tag), path("${species_tag}_norm_R1.fastq.gz"), path("${species_tag}_norm_R2.fastq.gz"),
          path("${species_tag}_norm_SE.fastq.gz"), emit: reads
    path("${species_tag}.se_candidates.csv"), optional: true, emit: se_candidates
    path("${species_tag}.blacklist_candidates.csv"), optional: true, emit: blacklist_candidates

    script:
    """

    # Output files must always exist (storeDir requirement).
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz

    # Read pre-queried accessions from SRA_QUERY CSV (up to max_rnaseq_runs).
    # Consult the project override file for per-accession actions:
    #   skip           — exclude entirely
    #   rename_headers — parallel-fastq-dump as normal, then seqkit renumbers headers
    #                    with sequential integers so R1/R2 names match exactly.
    #                    Fixes BGI/MGISEQ read-name mismatch from block-parallel splits.
    # Accessions absent from the override file use the default parallel-fastq-dump path.
    BLACKLIST="${launchDir}/rnaseq_blacklist.csv"
    RAW_ACCESSIONS=\$(awk -F',' 'NR>1 {print \$3}' ${sra_query_csv} | head -n ${params.max_rnaseq_runs})
    TAXONID=\$(awk -F',' 'NR==2 {print \$2; exit}' ${sra_query_csv})

    # Helper: record an accession that yielded a non-empty _1 but no _2 (single-end data
    # mislabeled PAIRED in SRA). Writes a row in blacklist format so it can be pasted
    # straight into rnaseq_blacklist.csv (acc,species_tag,taxonid,SE_trinity); a trailing
    # spots column (SRA-reported, col 4 of the query CSV; empty if unknown) rides along as
    # info and is ignored by the blacklist parser. Collected to rnaseq_se_candidates.csv.
    flag_se_candidate() {
        local acc="\$1" spots
        spots=\$(awk -F',' -v a="\$acc" 'NR>1 && \$3==a {print \$4; exit}' ${sra_query_csv})
        echo "[SE_CANDIDATE] \$acc ${species_tag} spots=\${spots:-NA} — add to rnaseq_blacklist.csv as SE_trinity"
        echo "\$acc,${species_tag},\$TAXONID,SE_trinity,\${spots}" >> ${species_tag}.se_candidates.csv
    }

    # Helper: record an accession whose download failed outright (parallel-fastq-dump and the
    # EBI FTP fallback both produced nothing usable). Writes a row in rnaseq_blacklist.csv
    # column order (acc,species_tag,taxonid,skip) with action=skip so it can be pasted
    # straight into the blacklist after review. Trailing spots column is info only and is
    # ignored by the blacklist parser. Collected to rnaseq_blacklist_candidates.csv.
    flag_blacklist_candidate() {
        local acc="\$1" spots
        spots=\$(awk -F',' -v a="\$acc" 'NR>1 && \$3==a {print \$4; exit}' ${sra_query_csv})
        echo "[BLACKLIST_CANDIDATE] \$acc ${species_tag} spots=\${spots:-NA} — download failed; add to rnaseq_blacklist.csv as skip"
        echo "\$acc,${species_tag},\$TAXONID,skip,\${spots}" >> ${species_tag}.blacklist_candidates.csv
    }

    # Helper: look up the explicit override action for an accession from the blacklist
    # (col 4 = action: skip | rename_headers).  Returns empty string if not listed.
    acc_action() {
        local acc="\$1"
        [ -f "\$BLACKLIST" ] && awk -F',' -v a="\$acc" 'NR>1 && \$1==a {print \$4; exit}' "\$BLACKLIST" || true
    }

    # Helper: look up sequencing platform from the per-species sra_query CSV (col 5).
    # Returns the raw SRA platform string, e.g. BGISEQ, ILLUMINA.
    acc_platform() {
        local acc="\$1"
        awk -F',' -v a="\$acc" 'NR>1 && \$3==a {print \$5; exit}' ${sra_query_csv}
    }

    # Resolve the effective download strategy for an accession:
    #   explicit blacklist action takes precedence;
    #   BGISEQ platform auto-maps to rename_headers;
    #   everything else is empty (default parallel-fastq-dump path).
    acc_strategy() {
        local acc="\$1" action platform
        action=\$(acc_action "\$acc")
        if [ -n "\$action" ]; then
            echo "\$action"
        else
            platform=\$(acc_platform "\$acc")
            case "\$platform" in
                BGISEQ|BGI*|MGI*) echo "rename_headers" ;;
                *)                 echo "" ;;
            esac
        fi
    }

    # Compute the EBI FTP directory URL for an SRA accession.
    # Path formula: vol1/fastq/<first-6-chars>/<subdir>/<acc>
    # subdir is omitted for ≤9-char accessions; for longer ones it zero-pads
    # chars 10-onward of the accession string into a 3-char field.
    ebi_ftp_dir() {
        local acc="\$1" len prefix base
        len="\${#acc}"
        prefix="\${acc:0:6}"
        base="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/\${prefix}"
        if   [ "\$len" -le 9  ]; then printf '%s/%s'      "\$base" "\$acc"
        elif [ "\$len" -eq 10 ]; then printf '%s/00%s/%s'  "\$base" "\${acc:9:1}" "\$acc"
        elif [ "\$len" -eq 11 ]; then printf '%s/0%s/%s'   "\$base" "\${acc:9:2}" "\$acc"
        else                          printf '%s/%s/%s'    "\$base" "\${acc:9:3}" "\$acc"
        fi
    }

    # Partition accessions into skip / active.
    # SE_trinity entries are treated as skip here: they route to SRA_FETCH_SE instead.
    ACCESSIONS=""
    for ACC in \$RAW_ACCESSIONS; do
        STRATEGY=\$(acc_strategy "\$ACC")
        ACTION=\$(acc_action "\$ACC")
        if [ "\$STRATEGY" = "skip" ] || [ "\$ACTION" = "SE_trinity" ]; then
            echo "[INFO] Skipping accession \$ACC for ${species_tag} (action: \${ACTION:-\$STRATEGY})"
        else
            ACCESSIONS="\$ACCESSIONS \$ACC"
        fi
    done
    ACCESSIONS=\$(echo "\$ACCESSIONS" | xargs)   # trim whitespace

    if [ -z "\$ACCESSIONS" ]; then
        echo "[INFO] No paired-end RNA-seq runs found for ${species_tag} (no accessions in query CSV or all skipped)"
    else
        echo "[INFO] SRA accessions for ${species_tag}: \$ACCESSIONS"
        TMPDIR=\${SCRATCH:-/tmp}
        mkdir -p reads

        # Download and concatenate in accession order so R1/R2 stay matched.
        for ACC in \$ACCESSIONS; do
            STRATEGY=\$(acc_strategy "\$ACC")
            echo "[INFO] Downloading \$ACC (strategy: \${STRATEGY:-default}) ..."
	    maxspot=""
	    if [ ${params.max_rnaseq_spot} -gt "0" ]; then
	   	maxspot="-X ${params.max_rnaseq_spot}"
	    fi

            if [ "\$STRATEGY" = "rename_headers" ]; then
                # BGI/MGISEQ: parallel-fastq-dump block-splitting desynchronises spot IDs
                # between R1 and R2.  Download identically to the default path, then pipe
                # each file through seqkit to replace every header with a sequential integer
                # before handing off — guarantees R1 read N and R2 read N have matching names.
                parallel-fastq-dump --sra-id \$ACC --threads ${task.cpus} \
                    --outdir reads/ --split-files --gzip --tmpdir \$TMPDIR || {
                    echo "[WARN] Download failed for \$ACC (\$maxspot), skipping"
                    flag_blacklist_candidate "\$ACC"
                    continue
                }
                if [ -f reads/\${ACC}_1.fastq.gz ] && [ -f reads/\${ACC}_2.fastq.gz ]; then
                    seqkit replace -j ${task.cpus} -p '.+' -r '{nr}' reads/\${ACC}_1.fastq.gz \
                        | ${params.fastq_hdr_script} --read 1 /dev/stdin \
                            --max-reads ${params.max_rnaseq_reads} \
                        | pigz -c >> \$TMPDIR/${species_tag}_R1.fastq.gz
                    seqkit replace -j ${task.cpus} -p '.+' -r '{nr}' reads/\${ACC}_2.fastq.gz \
                        | ${params.fastq_hdr_script} --read 2 /dev/stdin \
                            --max-reads ${params.max_rnaseq_reads} \
                        | pigz -c >> \$TMPDIR/${species_tag}_R2.fastq.gz
                    rm reads/\${ACC}_[12].fastq.gz
                elif [ -s reads/\${ACC}_1.fastq.gz ]; then
                    flag_se_candidate "\$ACC"
                    rm -f reads/\${ACC}_1.fastq.gz
                else
                    echo "[WARN] Missing pair for \$ACC after download, skipping"
                    flag_blacklist_candidate "\$ACC"
                fi
            else
                parallel-fastq-dump --sra-id \$ACC --threads ${task.cpus} \
                    --outdir reads/ --split-files \$maxspot --gzip --tmpdir \$TMPDIR || {
                    echo "[WARN] Download failed for \$ACC (\$maxspot), skipping"
                    flag_blacklist_candidate "\$ACC"
                    continue
                }
                if [ -f reads/\${ACC}_1.fastq.gz ] && [ -f reads/\${ACC}_2.fastq.gz ]; then
                    parallel -j 2 ${params.fastq_hdr_script} --read {} reads/\${ACC}_{}.fastq.gz \
		    --max-reads ${params.max_rnaseq_reads} \
		    \\| pigz -c \\>\\> \$TMPDIR/${species_tag}_R{}.fastq.gz  ::: 1 2
                    rm reads/\${ACC}_[12].fastq.gz
                elif [ -s reads/\${ACC}_1.fastq.gz ]; then
                    flag_se_candidate "\$ACC"
                    rm -f reads/\${ACC}_1.fastq.gz
                else
                    echo "[WARN] Missing pair for \$ACC after download, skipping"
                    flag_blacklist_candidate "\$ACC"
                fi
            fi
        done
        rm -rf reads
        ENFORCE="${params.readlen_script}"
        [[ -x "\$ENFORCE" ]] || { echo "[ERROR] enforce_seqpair_readlen not found or not executable at \$ENFORCE"; exit 1; }
        if ! "\$ENFORCE" in=\$TMPDIR/${species_tag}_R1.fastq.gz \
                in2=\$TMPDIR/${species_tag}_R2.fastq.gz \
                out=\$TMPDIR/${species_tag}_trunc_R1.fastq.gz \
                out2=\$TMPDIR/${species_tag}_trunc_R2.fastq.gz minlen=75; then
            echo "[WARN] enforce_seqpair_readlen failed for ${species_tag} (likely parallel-fastq-dump read-name mismatch); retrying via EBI FTP..."
            rm -f \$TMPDIR/${species_tag}_R1.fastq.gz \$TMPDIR/${species_tag}_R2.fastq.gz \
                  \$TMPDIR/${species_tag}_trunc_R1.fastq.gz \$TMPDIR/${species_tag}_trunc_R2.fastq.gz
            mkdir -p reads_ebi
            for ACC in \$ACCESSIONS; do
                EBI_DIR=\$(ebi_ftp_dir "\$ACC")
                aria2c --max-connection-per-server=4 --min-split-size=1M --max-tries=3 --retry-wait=5 \
                    "\${EBI_DIR}/\${ACC}_1.fastq.gz" -d reads_ebi/ -o "\${ACC}_1.fastq.gz" || true
                aria2c --max-connection-per-server=4 --min-split-size=1M --max-tries=3 --retry-wait=5 \
                    "\${EBI_DIR}/\${ACC}_2.fastq.gz" -d reads_ebi/ -o "\${ACC}_2.fastq.gz" || true
                if [ ! -s reads_ebi/\${ACC}_2.fastq.gz ]; then
                    if [ -s reads_ebi/\${ACC}_1.fastq.gz ]; then
                        echo "[WARN] \$ACC: no paired-end R2 at EBI but R1 present (single-end data)"
                        flag_se_candidate "\$ACC"
                    else
                        echo "[WARN] \$ACC: no paired-end R2 at EBI (single-end or accession absent); skipping"
                        flag_blacklist_candidate "\$ACC"
                    fi
                    rm -f reads_ebi/\${ACC}_1.fastq.gz reads_ebi/\${ACC}_2.fastq.gz
                    continue
                fi
                if [ ! -s reads_ebi/\${ACC}_1.fastq.gz ]; then
                    echo "[WARN] \$ACC: R1 missing at EBI; skipping"
                    flag_blacklist_candidate "\$ACC"
                    rm -f reads_ebi/\${ACC}_2.fastq.gz
                    continue
                fi
                parallel -j 2 ${params.fastq_hdr_script} --read {} reads_ebi/\${ACC}_{}.fastq.gz \
                    --max-reads ${params.max_rnaseq_reads} \
                    \\| pigz -c \\>\\> \$TMPDIR/${species_tag}_R{}.fastq.gz  ::: 1 2
                rm -f reads_ebi/\${ACC}_[12].fastq.gz
            done
            rm -rf reads_ebi
            "\$ENFORCE" in=\$TMPDIR/${species_tag}_R1.fastq.gz \
                in2=\$TMPDIR/${species_tag}_R2.fastq.gz \
                out=\$TMPDIR/${species_tag}_trunc_R1.fastq.gz \
                out2=\$TMPDIR/${species_tag}_trunc_R2.fastq.gz minlen=75 \
                || { echo "[ERROR] enforce_seqpair_readlen failed for ${species_tag} even after EBI FTP fallback"; exit 1; }
        fi
        bbnorm.sh in=\$TMPDIR/${species_tag}_trunc_R1.fastq.gz in2=\$TMPDIR/${species_tag}_trunc_R2.fastq.gz \
            out1=\$TMPDIR/${species_tag}_norm_R1.fastq.gz \
            out2=\$TMPDIR/${species_tag}_norm_R2.fastq.gz target=30 ecc=t

        fastp   --in1 \$TMPDIR/${species_tag}_norm_R1.fastq.gz \
                --in2 \$TMPDIR/${species_tag}_norm_R2.fastq.gz \
                --out1 ${species_tag}_norm_R1.fastq.gz --out2 ${species_tag}_norm_R2.fastq.gz \
                --thread ${task.cpus} --detect_adapter_for_pe \
                --cut_front --cut_front_window_size 1 --cut_front_mean_quality 5 \
                --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 5 \
                --cut_right --cut_right_window_size 4 --cut_right_mean_quality 5 \
                --length_required 25

        rm \$TMPDIR/${species_tag}_*

        NPAIRS=\$(zcat ${species_tag}_norm_R1.fastq.gz 2>/dev/null | awk 'NR%4==1' | wc -l || echo 0)
        echo "[INFO] Combined \$NPAIRS normalized read pairs for ${species_tag}"

        mkdir -p "${launchDir}/rnaseq_reads"
    fi
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    echo "[STUB] SRA_FETCH for ${species_tag}"
    """
}

// Download and normalize single-end RNA-seq for species with no paired-end data,
// or where rnaseq_blacklist.csv lists accessions with action=SE_trinity (SRA metadata
// says PAIRED but the data is actually single-end).
// Outputs zero-byte PE stubs so all read channels carry the same 4-tuple shape.
// storeDir note: if a species was previously fetched via SRA_FETCH (PE), re-routing
// it to SRA_FETCH_SE requires deleting rnaseq_reads/${species_tag}_norm_* files first.
process SRA_FETCH_SE {
    label 'sra'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads"

    cpus   32
    memory '96 GB'
    time   '2h'

    input:
    tuple val(species_tag), path(sra_query_csv)

    output:
    tuple val(species_tag),
          path("${species_tag}_norm_R1.fastq.gz"),
          path("${species_tag}_norm_R2.fastq.gz"),
          path("${species_tag}_norm_SE.fastq.gz"), emit: reads

    script:
    """

    # Zero-byte PE stubs (this process produces SE output only).
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz

    BLACKLIST="${launchDir}/rnaseq_blacklist.csv"

    acc_action() {
        local acc="\$1"
        [ -f "\$BLACKLIST" ] && awk -F',' -v a="\$acc" 'NR>1 && \$1==a {print \$4; exit}' "\$BLACKLIST" || true
    }

    ebi_ftp_dir() {
        local acc="\$1" len prefix base
        len="\${#acc}"
        prefix="\${acc:0:6}"
        base="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/\${prefix}"
        if   [ "\$len" -le 9  ]; then printf '%s/%s'      "\$base" "\$acc"
        elif [ "\$len" -eq 10 ]; then printf '%s/00%s/%s'  "\$base" "\${acc:9:1}" "\$acc"
        elif [ "\$len" -eq 11 ]; then printf '%s/0%s/%s'   "\$base" "\${acc:9:2}" "\$acc"
        else                          printf '%s/%s/%s'    "\$base" "\${acc:9:3}" "\$acc"
        fi
    }

    # Collect SE accessions from the query CSV:
    #   SE_trinity entries: col 6 (layout) may be PAIRED in SRA but blacklist says SE_trinity.
    #                       Download with --split-files and take _1 only (real SE data).
    #   SINGLE layout entries: col 6 == SINGLE; pfd gives ACC.fastq.gz or ACC_1.fastq.gz.
    SE_ACCESSIONS=""
    while IFS=',' read -r stag tid acc spots platform layout rest; do
        [ "\$stag" = "species_tag" ] && continue
        ACTION=\$(acc_action "\$acc")
        if [ "\$ACTION" = "SE_trinity" ] || [ "\${layout:-PAIRED}" = "SINGLE" ]; then
            SE_ACCESSIONS="\$SE_ACCESSIONS \$acc"
        fi
    done < ${sra_query_csv}
    SE_ACCESSIONS=\$(echo "\$SE_ACCESSIONS" | tr ' ' '\\n' | grep -v '^\$' | head -n ${params.max_rnaseq_se_runs} | tr '\\n' ' ' | xargs)

    if [ -z "\$SE_ACCESSIONS" ]; then
        echo "[INFO] No SE accessions found for ${species_tag}"
        exit 0
    fi

    echo "[INFO] SE accessions for ${species_tag}: \$SE_ACCESSIONS"
    TMPDIR=\${SCRATCH:-/tmp}
    mkdir -p reads

    for ACC in \$SE_ACCESSIONS; do
        ACTION=\$(acc_action "\$ACC")
        echo "[INFO] Downloading \$ACC (SE mode, action: \${ACTION:-default}) ..."

        # Download with --split-files. For SE_trinity (mislabeled PAIRED) this yields _1/_2;
        # we take only _1 (the actual SE reads) and discard _2. For genuine SINGLE layout,
        # pfd typically produces ACC_1.fastq.gz or ACC.fastq.gz.
        parallel-fastq-dump --sra-id "\$ACC" --threads ${task.cpus} \\
            --outdir reads/ --split-files --gzip --tmpdir "\$TMPDIR" || {
            echo "[WARN] pfd failed for \$ACC; trying EBI FTP..."
            EBI_DIR=\$(ebi_ftp_dir "\$ACC")
            # SE_trinity (mislabeled PAIRED): EBI has ACC_1.fastq.gz
            # Genuine SINGLE: EBI has ACC.fastq.gz
            aria2c --max-connection-per-server=4 --min-split-size=1M --max-tries=3 --retry-wait=5 \\
                "\${EBI_DIR}/\${ACC}_1.fastq.gz" -d reads/ -o "\${ACC}_1.fastq.gz" 2>/dev/null || true
            if [ ! -s "reads/\${ACC}_1.fastq.gz" ]; then
                aria2c --max-connection-per-server=4 --min-split-size=1M --max-tries=3 --retry-wait=5 \\
                    "\${EBI_DIR}/\${ACC}.fastq.gz" -d reads/ -o "\${ACC}.fastq.gz" 2>/dev/null || true
            fi
        }

        SE_FILE=""
        if   [ -s "reads/\${ACC}_1.fastq.gz" ]; then SE_FILE="reads/\${ACC}_1.fastq.gz"
        elif [ -s "reads/\${ACC}.fastq.gz"   ]; then SE_FILE="reads/\${ACC}.fastq.gz"
        fi

        if [ -n "\$SE_FILE" ]; then
            ${params.fastq_hdr_script} --read 1 "\$SE_FILE" \\
                --max-reads ${params.max_rnaseq_reads} \\
                | pigz -c >> "\$TMPDIR/${species_tag}_SE.fastq.gz"
            rm -f "reads/\${ACC}_1.fastq.gz" "reads/\${ACC}_2.fastq.gz" "reads/\${ACC}.fastq.gz"
        else
            echo "[WARN] No SE file found for \$ACC after download; skipping"
            rm -f "reads/\${ACC}"*.fastq.gz
        fi
    done
    rm -rf reads

    if [ ! -s "\$TMPDIR/${species_tag}_SE.fastq.gz" ]; then
        echo "[WARN] No SE reads downloaded for ${species_tag}; leaving empty SE placeholder"
        exit 0
    fi

    bbnorm.sh in="\$TMPDIR/${species_tag}_SE.fastq.gz" \\
        out="\$TMPDIR/${species_tag}_norm_SE_bbn.fastq.gz" target=30 ecc=f

    fastp --in1 "\$TMPDIR/${species_tag}_norm_SE_bbn.fastq.gz" \\
          --out1 ${species_tag}_norm_SE.fastq.gz \\
          --thread ${task.cpus} \\
          --cut_front --cut_front_window_size 1 --cut_front_mean_quality 5 \\
          --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 5 \\
          --cut_right --cut_right_window_size 4 --cut_right_mean_quality 5 \\
          --length_required 25

    rm -f "\$TMPDIR/${species_tag}_SE.fastq.gz" "\$TMPDIR/${species_tag}_norm_SE_bbn.fastq.gz"

    NREADS=\$(zcat ${species_tag}_norm_SE.fastq.gz 2>/dev/null | awk 'NR%4==1' | wc -l || echo 0)
    echo "[INFO] \$NREADS normalized SE reads for ${species_tag}"

    mkdir -p "${launchDir}/rnaseq_reads"
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    echo "[STUB] SRA_FETCH_SE for ${species_tag}"
    """
}

// Run funannotate train on the representative (first) assembly of each species, then
// archive the Trinity-GG transcripts (normalized reads are in rnaseq_reads)
// reads into rnaseq_data/ so all other strains can skip those expensive steps.
// storeDir skips this process entirely if all five output files already exist.
process RNASEQ_PREPARE {
    label 'funannotate'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_data"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(species_tag), val(meta), val(genome_fa), path(r1), path(r2), path(se)

    output:
    tuple val(species_tag),
            path("${species_tag}.trinity-GG.fasta"), emit: shared

    script:
    def out           = meta.id
    def species       = meta.species
    def strain        = meta.strain
    def header_length = params.header_length
    """
    # ── Empty-reads sentinel: no RNA-seq found by SRA_FETCH / SRA_FETCH_SE ──
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ]; then
        echo "[INFO] No RNAseq reads for ${species_tag}; writing empty shared markers"
        touch ${species_tag}.trinity-GG.fasta
        exit 0
    fi

    # ── If representative was already trained, just extract shared files ──────
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    if [ -f "\$TRAIN_GFF3" ]; then
        echo "[INFO] Training already complete for ${out}; extracting shared files to rnaseq_data"
        TRAINDIR="${params.training_target}/${out}/training"
        TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
        if [ -n "\$TRINITY_FA" ]; then
            cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
        else
            touch ${species_tag}.trinity-GG.fasta
        fi
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    # ── Run full funannotate train on the representative genome ───────────────
    # Use SCRATCH for the funannotate output dir so Trinity/HISAT2/normalize
    # intermediates land on fast local storage and don't consume project quota.
    echo "[INFO] RNASEQ_PREPARE: running funannotate train for representative ${out} (species: ${species_tag})"

    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    if [ -s "${r1}" ]; then
        funannotate train -i "\$GENOME_IN" -o \$SCRATCH/${out} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            --stop_after_trinity --no_trimmomatic
    else
        echo "[INFO] RNASEQ_PREPARE: using single-end reads for ${out}"
        funannotate train -i "\$GENOME_IN" -o \$SCRATCH/${out} \\
            --single_norm ${se} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            --stop_after_trinity --no_trimmomatic
    fi

    # ── Copy shared outputs to rnaseq_data/ ──────────────────────────────────
    TRAINDIR="\$SCRATCH/${out}/training"
    TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
    if [ -n "\$TRINITY_FA" ]; then
        cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
    else
        echo "[WARN] No trinity.fasta found under \$TRAINDIR for ${out}"
        touch ${species_tag}.trinity-GG.fasta
    fi

    # ── Clean up scratch output dir (all intermediates were temporary) ────────
    rm -rf "\$SCRATCH/${out}"
    echo "[INFO] RNASEQ_PREPARE complete for ${species_tag}"
    """

    stub:
    def out = meta.id
    """
    echo ">stub_trinity_${species_tag}" > ${species_tag}.trinity-GG.fasta
    mkdir -p ${params.training_target}/${out}/training
    touch ${params.training_target}/${out}/training/funannotate_train.pasa.gff3
    """
}

// For non-representative strains: funannotate train --trinity <shared_fasta> runs only
// PASA (skips Trimmomatic, normalization, HISAT2, and Trinity-GG assembly).
// Falls back to a full train when no shared Trinity is available (e.g. species with
// a single strain or when run_sra_fetch is false).
process FUNANNOTATE_TRAIN {
    label 'funannotate'
    tag "${meta.id}"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(meta), val(genome_fa), path(r1), path(r2), path(se), path(trinity_fa)

    output:
    tuple val(meta), val(genome_fa)

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def header_length = params.header_length
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no RNA-seq data at all ────────────────────────────────────────
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ] && [ ! -s "${trinity_fa}" ]; then
        echo "[INFO] No RNAseq data for ${out}, skipping funannotate train"
        exit 0
    fi

    # ── Skip if training output already present and rnaseq is not newer than GBK ──
    # Accept a compressed prediction (.gbk.gz) as "done" so folders can be space-saved.
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk"
    [ -f "\$PREDICT_GBK" ] || PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk.gz"
    if [ -f "\$TRAIN_GFF3" ]; then
        RETRAIN=0
        if [ -f "\$PREDICT_GBK" ]; then
            # Re-train if the rnaseq reads are newer than the existing prediction GBK.
            if [ -s "${r1}" ] && [ "${r1}" -nt "\$PREDICT_GBK" ]; then
                echo "[INFO] RNAseq R1 reads newer than predict GBK for ${out}; retraining"
#                rm -rf "${params.training_target}/${out}/training"
                RETRAIN=1
            elif [ -s "${se}" ] && [ "${se}" -nt "\$PREDICT_GBK" ]; then
                echo "[INFO] RNAseq SE reads newer than predict GBK for ${out}; retraining"
#                rm -rf "${params.training_target}/${out}/training"
                RETRAIN=1
            fi
        fi
        if [ \$RETRAIN -eq 0 ]; then
            echo "[INFO] Training already complete for ${out}; skipping"
            exit 0
        fi
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load funannotate/dev-1.8.18 
    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    pasa_db_arg="--pasa_db sqlite"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        MYSQL_SCRATCH=${params.training_target}/${out}/training/mysql_db
        if [ ! -f \$MYSQL_SCRATCH/mysql/conf/my.cnf ]; then
            echo "[INFO] Setting up temporary MariaDB for PASA at \$MYSQL_SCRATCH"
            mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
            rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
                { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
            cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
                { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        fi
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        # ──  may be unnecessary if overridden by -B option later? ──
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH/mysql_db
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        module load singularity
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb_${asmid} /usr/bin/mysqld_safe
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # ── Use shared Trinity transcripts (PASA only) or run full train ──────────
    if [ -s "${trinity_fa}" ]; then
        if [ -s "${r1}" ]; then
            echo "[INFO] Running funannotate train (PASA+PE) for ${out} using shared Trinity"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --left_norm ${r1} --right_norm ${r2} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        elif [ -s "${se}" ]; then
            echo "[INFO] Running funannotate train (PASA+SE) for ${out} using shared Trinity"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --single_norm ${se} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        else
            echo "[INFO] Running funannotate train (PASA only, no reads) for ${out} using shared Trinity"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --left_norm ${r1} --right_norm ${r2} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        fi
    elif [ -s "${r1}" ]; then
        echo "[INFO] Running funannotate train (full PE, no shared Trinity) for ${out}"
        funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    else
        echo "[INFO] Running funannotate train (full SE, no shared Trinity) for ${out}"
        funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
            --single_norm ${se} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    fi

    # ── Remove large intermediates not needed for predict or update ─────────────
    # Keeps: *.bam, *.bai, *.pasa.gff3, *.stringtie.gtf, *.transcripts.gff3
    TRAINDIR="${params.training_target}/${out}/training"
    echo "[INFO] Removing large training intermediates in \$TRAINDIR"
    rm -rf "\$TRAINDIR/hisat2"
    rm -rf "\$TRAINDIR/trinity_gg"
    echo "[INFO] Training cleanup complete for ${out}"
    echo "mysql is ${params.pasa_mysql}"
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] FUNANNOTATE_TRAIN stub for ${out}"
    mkdir -p ${params.training_target}/${out}/training
    touch ${params.training_target}/${out}/training/funannotate_train.pasa.gff3
    """
}

// Option B persistence model: funannotate predict computes DIRECTLY into the persistent
// per-genome dir (${params.target}/${out}), symmetric with FUNANNOTATE_TRAIN writing to
// training_target. funannotate checkpoints into predict_misc/, so a restart after an
// OOM/timeout/orchestrator death resumes completed steps in place rather than starting
// over. There is no publishDir copy and no work-dir<->target rsync: the durable output is
// written where downstream steps already read it. Large intermediates still go to the
// node-local --tmpdir. The Nextflow output is a small marker file (nothing consumes the
// predict dir as a channel; downstream rebuilds metadata from the CSV and gates on the
// on-disk GBK), so emitting a marker keeps the DAG edge without copying the result tree.
process FUNANNOTATE_PREDICT {
    label 'funannotate'
    tag "${meta.id}"

    cpus   16
    memory '32 GB'
    time   '32h'

    input:
    tuple val(meta), val(genome_fa)

    output:
    val meta, emit: metadata
    path("${meta.id}.predict.done"), emit: done

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    def header_length = params.header_length
    def transl_table  = meta.transl_table
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load funannotate/dev-1.8.18
    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    PREDICTDIR="${params.target}/${out}"
    PREDICT_GBK="\$PREDICTDIR/predict_results/${out}.gbk"

    if [ "${params.debug.toBoolean()}" = "true" ]; then
        echo "[DEBUG] out=${out} asmid=${asmid} species=${species} strain=${strain}"
        echo "[DEBUG] locustag=${locustag} busco=${busco_lineage} transl_table=${transl_table}"
        echo "[DEBUG] proteins=${params.proteins} genome_fa=${genome_fa}"
        echo "[DEBUG] PREDICTDIR=\$PREDICTDIR TMPDIR=\$TMPDIR pwd=\$(pwd)"
    fi

    # ── Skip vs. refresh decision ─────────────────────────────────────────────
    # The workflow schedules this process when the GBK is missing OR stale (rnaseq/trinity
    # newer than the GBK, per staleRnaseq()). Re-derive staleness here from the same on-disk
    # timestamps so a current GBK short-circuits, but a stale one forces a clean re-predict.
    if [ -s "\$PREDICT_GBK" ]; then
        SPECIES_TAG=\$(printf '%s' "${species}" | sed -E 's/[[:space:]]+/_/g')
        STALE=0
        for f in "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_R1.fastq.gz" \\
                 "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_SE.fastq.gz" \\
                 "${launchDir}/rnaseq_data/\${SPECIES_TAG}.trinity-GG.fasta"; do
            if [ -s "\$f" ] && [ "\$f" -nt "\$PREDICT_GBK" ]; then STALE=1; fi
        done
        if [ "\$STALE" -eq 0 ]; then
            echo "[INFO] Prediction already complete and current for ${out}; nothing to do"
            touch ${out}.predict.done
            exit 0
        fi
        echo "[INFO] Stale prediction for ${out}: rnaseq/trinity newer than GBK — clearing predict outputs for a fresh run"
        rm -rf "\$PREDICTDIR/predict_results" "\$PREDICTDIR/predict_misc"
    fi

    mkdir -p "\$PREDICTDIR"

    # ── Guard against a corrupt partial from a previous attempt ───────────────
    # funannotate resumes from predict_misc/. If predict_results/ exists without a
    # predict_misc/ (a half-written tree with no checkpoints and no GBK), clear it so
    # predict starts the consensus/output step from a clean state instead of choking on it.
    if [ ! -d "\$PREDICTDIR/predict_misc" ] && [ -d "\$PREDICTDIR/predict_results" ]; then
        echo "[WARN] predict_results/ present without predict_misc/ for ${out}; clearing stale partial"
        rm -rf "\$PREDICTDIR/predict_results"
    fi

    # funannotate predict expects training data at <outdir>/training; point it at the
    # persistent training dir. The symlink lives in the persistent project tree (no
    # publishDir to recursively copy the target), so it is left in place.
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "\$PREDICTDIR/training"
    fi

    TBL2ASN_PARAMS="-l paired-ends"

    # Inflate a gzipped clean/masked genome to a local uncompressed copy; funannotate
    # cannot read a gzipped FASTA via -i, and the pre-flight awk below also needs plain
    # text. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # ── Too-small-genome pre-flight guard ────────────────────────────────────
    # Assemblies that are both small AND fragmented cannot yield funannotate's
    # required 30 training models; predict would run for hours then abort with
    # "Not enough gene models N to train Augustus (30 required), exiting". Detect
    # that up front from cheap contig stats and skip cleanly (flag, no crash).
    # Requires BOTH gates so complete small genomes (e.g. Malassezia) are unaffected.
    # Disabled when predict_min_asm_bp=0.
    SKIP_REPORT="${params.target}/predict_skipped_too_small.tsv"
    if [ "${params.predict_min_asm_bp}" -gt 0 ]; then
        # Per-contig lengths -> sort descending -> N50 (portable; no gawk asort).
        read ASM_BP ASM_CTG ASM_N50 < <(
            awk '/^>/{if(len)print len;len=0;next}{len+=length(\$0)}END{if(len)print len}' "\$GENOME_IN" \\
            | sort -rn \\
            | awk '{L[NR]=\$1;tot+=\$1}END{half=tot/2;run=0;n50=0;for(i=1;i<=NR;i++){run+=L[i];if(run>=half){n50=L[i];break}}print tot, NR, n50}')
        echo "[INFO] Pre-flight assembly stats for ${out}: \${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}"
        SMALL=0; FRAG=0
        [ "\$ASM_BP" -lt "${params.predict_min_asm_bp}" ] && SMALL=1
        { [ "\$ASM_N50" -lt "${params.predict_frag_max_n50}" ] || [ "\$ASM_CTG" -gt "${params.predict_frag_max_contigs}" ]; } && FRAG=1
        if [ "\$SMALL" -eq 1 ] && [ "\$FRAG" -eq 1 ]; then
            echo "[WARN] ${out} is too small/fragmented for funannotate training (\${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}); skipping predict" >&2
            mkdir -p "${params.target}"
            [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "preflight_small_fragmented" "\$ASM_BP" "\$ASM_CTG" "\$ASM_N50" >> "\$SKIP_REPORT"
            touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
            touch ${out}.predict.done
            exit 0
        fi
    fi

    funannotate predict --name ${locustag} -i "\$GENOME_IN" --strain "${strain}" \\
        -o "\$PREDICTDIR" -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w codingquarry:0 glimmerhmm:0 \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence ${params.proteins} \\
        --max_intronlen ${params.max_intronlen} --min_intronlen ${params.min_intronlen} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table} --auto-skip-genemark || true

    # ── Post-predict catch ────────────────────────────────────────────────────
    # If predict produced no GBK, distinguish the known "too few training models"
    # outcome (an unfixable property of the assembly) from a genuine error. The
    # former is flagged and skipped so it does not abort the batch; anything else
    # still hard-fails so real problems surface.
    if [ ! -s "\$PREDICT_GBK" ]; then
        PLOG="\$PREDICTDIR/logfiles/funannotate-predict.log"
        if [ -f "\$PLOG" ] && grep -q "Not enough gene models .* to train Augustus" "\$PLOG"; then
            NMODELS=\$(grep -oE "Not enough gene models [0-9]+" "\$PLOG" | grep -oE "[0-9]+" | tail -1)
            echo "[WARN] ${out}: funannotate found only \${NMODELS:-<min} training models (needs 30); too small/fragmented to annotate — skipping" >&2
            mkdir -p "${params.target}"
            [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "funannotate_too_few_models:\${NMODELS:-NA}" "" "" "" >> "\$SKIP_REPORT"
            touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
            touch ${out}.predict.done
            exit 0
        fi
        echo "ERROR: funannotate predict did not produce expected GBK: \$PREDICT_GBK" >&2
        exit 1
    fi
    if [ -d "\$PREDICTDIR/predict_misc/ab_initio_parameters" ]; then
        mv "\$PREDICTDIR/predict_misc/ab_initio_parameters" "\$PREDICTDIR"
        mv "\$PREDICTDIR/predict_misc/trnascan.no-overlaps.gff3" "\$PREDICTDIR"
        rm -rf "\$PREDICTDIR/predict_misc"
        mkdir -p "\$PREDICTDIR/predict_misc"
        mv "\$PREDICTDIR/ab_initio_parameters" "\$PREDICTDIR/trnascan.no-overlaps.gff3" "\$PREDICTDIR/predict_misc"
    fi
    find "\$PREDICTDIR/predict_results/" -maxdepth 1 \\( -name "*.txt" -o -name "*.mrna-transcripts.fa" \\) -print0 \
        | xargs -0 --no-run-if-empty pigz
    sync
    touch ${out}.predict.done
    echo "[INFO] Prediction complete for ${out} at \$PREDICTDIR"
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] Would run funannotate predict for ${out} using ${genome_fa}"
    [ -f "${genome_fa}" ] || [ -f "${genome_fa}.gz" ] || { echo "ERROR: genome not found at ${genome_fa}[.gz]" >&2; exit 1; }
    mkdir -p ${params.target}/${out}/predict_results ${params.target}/${out}/predict_misc
    # non-empty so downstream size>0 gating (predict_ch / postpredict) is exercised
    echo "LOCUS stub_${out}" > ${params.target}/${out}/predict_results/${out}.gbk
    echo ">stub_${out}_p1" > ${params.target}/${out}/predict_results/${out}.proteins.fa
    touch ${out}.predict.done
    """
}

process FUNANNOTATE_ANNOTATE {
    label 'funannotate'
    tag "${meta.id}"

    cpus   16
    memory '32 GB'
    time   '48h'

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}.annotate.done"), emit: marker

    script:
    def out           = meta.id
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    def header_length = params.header_length
    def antiSm    = file("${params.target}/${meta.id}/antismash_local/${meta.id}.gbk")
    def antiSmArg = antiSm.exists() ? "--antismash ${antiSm}" : ""
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    funannotate annotate -i ${params.target}/${out} -o ${params.target}/${out} \\
        --species "${species}" --strain "${strain}" \\
        --busco_db ${busco_lineage} --rename ${locustag} \\
        --sbt ${params.sbt_template} \\
        --header_length ${header_length} \\
        ${antiSmArg} \\
        --cpu ${task.cpus} --tmpdir \$TMPDIR

    EXPECTED_GBK="${params.target}/${out}/annotate_results/${out}.gbk"
    if [ ! -f "\$EXPECTED_GBK" ]; then
        echo "ERROR: funannotate annotate did not produce expected GBK: \$EXPECTED_GBK" >&2
        exit 1
    fi
    touch ${out}.annotate.done
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] Would run funannotate annotate for ${out}"
    mkdir -p ${params.target}/${out}/annotate_results ${params.target}/${out}/annotate_misc
    touch ${params.target}/${out}/annotate_results/${out}.gbk
    touch ${out}.annotate.done
    """
}

process FUNANNOTATE_UPDATE {
    label 'funannotate'
    tag "${meta.id}"

    cpus   16
    memory '96 GB'
    time   '48h'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    val meta

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    def header_length = params.header_length
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no reads (empty marker file from SRA_FETCH) ──────────────────
    if [ ! -s "${r1}" ]; then
        echo "[INFO] No RNAseq reads for ${out}, skipping funannotate update"
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    pasa_db_arg="--pasa_db sqlite"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        MYSQL_SCRATCH=${params.training_target}/${out}/training/mysql_db
        if [ ! -f \$MYSQL_SCRATCH/mysql/conf/my.cnf ]; then
            echo "[INFO] Setting up temporary MariaDB for PASA at \$MYSQL_SCRATCH"
            mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
            rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
                { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
            cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
                { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        fi
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH/mysql_db
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        module load singularity
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb_${asmid} /usr/bin/mysqld_safe
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # Link training data into work dir so funannotate update finds it at the relative path it expects.
    mkdir -p ${out}
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "${out}/training"
    fi

    # r1/r2 are pre-normalized reads from SRA_FETCH (fastp-trimmed + bbnorm-normalized).
    # funannotate update will still run its internal alignment step against these.
    echo "[INFO] Running funannotate update for ${out}"
    funannotate update -i ${params.target}/${out} \\
        --left ${r1} --right ${r2} \\
        --cpus ${task.cpus} \\
        \$pasa_db_arg
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    EXPECTED="${params.target}/${out}/update_results/${out}.gbk"
    if [ ! -f "\$EXPECTED" ]; then
        echo "ERROR: funannotate update did not produce expected GBK: \$EXPECTED" >&2
        exit 1
    fi
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] FUNANNOTATE_UPDATE stub for ${out} (r1=${r1}, r2=${r2})"
    mkdir -p ${params.target}/${out}/update_results
    touch ${params.target}/${out}/update_results/${out}.tbl
    touch ${params.target}/${out}/update_results/${out}.gbk
    touch ${params.target}/${out}/update_results/${out}.gff3
    """
}

// A funannotate step's GenBank output may be stored uncompressed (.gbk) or
// gzip-compressed (.gbk.gz) so completed folders can be compressed to save space.
// Returns the existing non-empty file (preferring .gbk), or null if neither exists.
// Use this for completion/skip gating so a compressed result still counts as "done".
def gbkResult(String dir, String out) {
    def plain = file("${dir}/${out}.gbk")
    if (plain.exists() && plain.size() > 0) return plain
    def gz = file("${dir}/${out}.gbk.gz")
    if (gz.exists() && gz.size() > 0) return gz
    return null
}

// Clean/masked genomes in input_clean_genomes may be stored gzip-compressed (.gz) to
// save space. Given the uncompressed base path (e.g. .../<asmid>.fa or
// .../<asmid>.masked.fasta), returns the existing non-empty file, preferring the
// compressed form. Falls back to the plain path object when neither exists, so callers'
// .exists() checks still report missing.
def genomeFile(String base) {
    def gz = file("${base}.gz")
    if (gz.exists() && gz.size() > 0) return gz
    return file(base)
}

def staleRnaseq(String out, String species) {
    def species_tag = species.replaceAll(/\s+/, '_')
    def gbk = gbkResult("${params.target}/${out}/predict_results", out)
    if (gbk == null) return false  // predict hasn't run yet; normal path handles it
    def r1      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_R1.fastq.gz")
    def se      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_SE.fastq.gz")
    def trinity = file("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
    def r1_newer      = r1.exists()      && r1.size() > 0      && r1.lastModified()      > gbk.lastModified()
    def se_newer      = se.exists()      && se.size() > 0      && se.lastModified()      > gbk.lastModified()
    def trinity_newer = trinity.exists() && trinity.size() > 0 && trinity.lastModified() > gbk.lastModified()
    if (r1_newer || se_newer || trinity_newer) {
        log.info "stale prediction for ${out}: rnaseq/trinity newer than GBK — scheduling retrain+repredict"
        return true
    }
    return false
}

include { validateParameters; paramsSummaryLog; paramsHelp } from 'plugin/nf-schema'
include { ASM_STATS }        from './modules/local/asm_stats'
include { INPUT_CHECK }      from './subworkflows/local/input_check'
include { SETUP_DBS }        from './subworkflows/local/setup_dbs'
include { CLEAN_GENOMES }    from './subworkflows/local/clean_genomes'
include { ANTISMASH_RUN }    from './modules/local/antismash_run'
include { INTERPROSCAN_RUN } from './modules/local/interproscan_run'
include { SIGNALP_RUN }      from './modules/local/signalp_run'

workflow {
    // `--help` prints schema-driven parameter help (grouped, with types/defaults) and exits.
    if (params.help) {
        log.info paramsHelp()
        exit 0
    }
    // Type-check params against nextflow_schema.json and log the resolved set.
    // (Unrecognised params warn rather than fail — see nextflow.config.)
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Fail fast with an actionable message when a pipeline profile was not selected
    // (these params come from conf/profile_annotate.config). Without it, downstream
    // file(params.funannotate_db) calls throw a cryptic "file() ... cannot be null".
    if( !params.taxondb || !params.funannotate_db )
        error "Missing params.taxondb / params.funannotate_db — add a pipeline profile, e.g. -profile annotate,slurm,module (or use: sbatch nextflow/run_annotate.sh)"

    // ── Samplesheet ingestion (INPUT_CHECK) ──────────────────────────────────
    // Parses samples CSV, applies taxon/asmid/suppress/n_test filters, builds
    // meta maps, and resolves genome paths. Two outputs:
    //   jobs        — tuple(meta, gz)  with genome existence filter (cleaning path)
    //   postpredict — meta only        no genome filter (annotate/update paths)
    INPUT_CHECK()
    def jobs = INPUT_CHECK.out.genomes

    def ch_versions = Channel.empty()
    if (params.debug.toBoolean()) {
        jobs.view { meta, gz -> "[CHANNEL] Submitting: out=${meta.id}, asmid=${meta.asmid}, transl_table=${meta.transl_table}, gz=${gz}" }
    }

    // Build/seed the three run-once databases. All use storeDir so they are no-ops
    // on any run where their target directories already exist.
    SETUP_DBS()
    def taxondb_ch = SETUP_DBS.out.taxondb

    CLEAN_GENOMES(jobs, taxondb_ch)

    if (!params.only_clean.toBoolean()) {
        def clean_genome_ch = CLEAN_GENOMES.out.genomes

        // ── Generate assembly statistics (for earlgrey_mask.nf SELECT_REPS) ────────
        // Generate asm_stats.tsv if --gen_asm_stats is true and the file doesn't exist.
        // This is used by earlgrey_mask.nf to select representative genomes per species.
        if (params.gen_asm_stats.toBoolean()) {
            def asm_stats_path = file(params.tables_dir).toAbsolutePath()
            def asm_stats_gz = file("${asm_stats_path}/asm_stats.tsv.gz")
            if (!asm_stats_gz.exists()) {
                log.info "Generating assembly statistics: ${asm_stats_gz}"
                ASM_STATS(
                    file(params.samples),
                    file(params.genome_dir)
                )
                ch_versions = ch_versions.mix(ASM_STATS.out.versions)
            } else {
                log.info "Assembly statistics already exist: ${asm_stats_gz}"
            }
        }

        // ── Repeat masking ────────────────────────────────────────────────────────
        // predict_genome_ch carries the genome path to use for prediction — either
        // the tantan soft-masked genome (default) or the clean unmasked genome
        // (--run_repeatmasker false).
        def predict_genome_ch
        if (params.run_repeatmasker.toBoolean()) {
            MASKREPEAT_TANTAN_RUN(clean_genome_ch)
            predict_genome_ch = MASKREPEAT_TANTAN_RUN.out.masked
                .map { meta, masked_fa ->
                    tuple(meta, masked_fa.toAbsolutePath().toString())
                }
        } else {
            // --run_repeatmasker false: use masked genome if a prior run produced it, else unmasked.
            predict_genome_ch = clean_genome_ch
                .map { meta, genome_fa ->
                    def masked = genomeFile("${launchDir}/input_clean_genomes/${meta.asmid}.masked.fasta")
                    def use_fa = masked.exists() ? masked.toString() : genome_fa
                    if (params.debug.toBoolean()) {
                        log.info "[DEBUG] ${meta.asmid}: genome_fa=${use_fa} (masked=${masked.exists()})"
                    }
                    tuple(meta, use_fa)
                }
        }

        // Gate the predict chain on funannotate DB + augustus config being ready.
        // SETUP_DBS was already called above; its storeDir-cached outputs are free
        // on resumed runs. Gating here threads the dependency through the entire
        // downstream funannotate subgraph (train, predict, update, annotate).
        // (MASKREPEAT uses `funannotate mask`, which needs neither, so it is intentionally
        // left ungated and can run in parallel with these setup steps.)
        predict_genome_ch = predict_genome_ch
            .combine(SETUP_DBS.out.db)
            .combine(SETUP_DBS.out.config)
            .map { row -> row[0..-3] }

        // FUNANNOTATE_PREDICT input tuple drops taxonid (not needed after masking/clean).
        // When SRA is enabled: SRA_FETCH fetches reads once per species; RNASEQ_PREPARE runs
        // funannotate train on the representative assembly and archives Trinity-GG, trimmed, and
        // normalized reads to rnaseq_data/; all other strains run FUNANNOTATE_TRAIN --trinity.
        def predict_input_ch
        def reads_ch = Channel.empty()
        if (params.run_sra_fetch.toBoolean()) {
            // Build per-species input: group assemblies, keep first taxonid per species.
            def sra_input = predict_genome_ch
                .map { meta, genome_fa ->
                    def species_tag = meta.species.replaceAll(/\s+/, '_')
                    tuple(species_tag, meta.taxonid)
                }
                .groupTuple(by: 0)
                .map { species_tag, taxonids -> tuple(species_tag, taxonids[0]) }

            // Step 1: query or reuse cached per-species SRA query results.
            // skip_sra_query=true reads existing CSVs from rnaseq_reads/sra_query/ directly,
            // bypassing SRA_QUERY_BATCH entirely (no SLURM jobs submitted).
            def sra_query_results
            if (params.skip_sra_query.toBoolean()) {
                sra_query_results = sra_input
                    .map { species_tag, _taxonid ->
                        def csv = file("${launchDir}/rnaseq_reads/sra_query/${species_tag}.sra_query.csv")
                        if (!csv.exists()) {
                            log.warn "skip_sra_query: no cached CSV for ${species_tag} — skipping this species"
                            return null
                        }
                        tuple(species_tag, csv)
                    }
                    .filter { it != null }
            } else {
                def sra_batched = sra_input
                    .collate(params.sra_query_batch_size)
                    .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }) }
                SRA_QUERY_BATCH(sra_batched)
                sra_query_results = SRA_QUERY_BATCH.out.query_results
                    .flatten()
                    .map { csv -> tuple(csv.baseName.replaceAll(/\.sra_query$/, ''), csv) }
            }

            // Step 2: Collect all per-species results into {stem}.rnaseq_sra.csv
            def stem = file(params.samples).baseName
            COLLECT_SRA_QUERY(
                sra_query_results.map { _stag, csv -> csv }.collect(),
                stem
            )

            if (!params.stop_after_sra_query.toBoolean()) {
            // Step 3: Classify each species CSV for routing.
            // Read blacklist once so closures below can check SE_trinity accessions.
            // Uses a Map<accession, action> for O(1) lookup.
            def blPath = file("${launchDir}/rnaseq_blacklist.csv")
            def blMap = blPath.exists()
                ? blPath.readLines().drop(1)
                      .findAll { it.trim() && !it.startsWith('#') }
                      .collectEntries { line ->
                          def cols = line.split(',')
                          cols.size() >= 4 ? [(cols[0].trim()): cols[3].trim()] : [:]
                      }
                : [:]

            // csvHasPE: CSV has at least one PAIRED accession not blocked or overridden to SE.
            def csvHasPE = { csv ->
                csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                    def cols = line.split(',')
                    if (cols.size() < 3) return false
                    def layout = cols.size() > 5 ? cols[5].trim() : 'PAIRED'
                    def action = blMap.get(cols[2].trim(), '')
                    layout == 'PAIRED' && action != 'skip' && action != 'SE_trinity'
                }
            }

            // csvHasSEtrinity: CSV has PAIRED accessions overridden to SE via SE_trinity blacklist.
            // These bypass the enable_single_end gate — they are a manual per-accession override.
            def csvHasSEtrinity = { csv ->
                csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                    def cols = line.split(',')
                    cols.size() >= 3 && blMap.get(cols[2].trim(), '') == 'SE_trinity'
                }
            }

            // csvHasSingleLayout: CSV has at least one genuine SINGLE-layout accession.
            // Only active when enable_single_end=true.
            def csvHasSingleLayout = { csv ->
                csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                    def cols = line.split(',')
                    cols.size() > 5 && cols[5].trim() == 'SINGLE' && blMap.get(cols[2].trim(), '') != 'skip'
                }
            }

            // Three-way branch:
            //   has_pe  → SRA_FETCH  (PE wins; SE_trinity entries ignored here, handled by SRA_FETCH)
            //   has_se  → SRA_FETCH_SE (SE_trinity always; SINGLE layout only if enable_single_end)
            //   no_data → WRITE_EMPTY_READS
            def branched_sra = sra_query_results
                .branch {
                    has_pe: csvHasPE.call(it[1])
                    has_se: csvHasSEtrinity.call(it[1]) ||
                            (params.enable_single_end.toBoolean() && csvHasSingleLayout.call(it[1]))
                    no_data: true
                }

            SRA_FETCH(branched_sra.has_pe)
            SRA_FETCH_SE(branched_sra.has_se)
            WRITE_EMPTY_READS(branched_sra.no_data.map { stag, _csv -> stag })

            // Accessions found to be single-end (non-empty _1, no _2) during the PE fetch are
            // recorded as blacklist-ready rows; merge all per-task notes into one reviewable
            // file at the project root. Add these to rnaseq_blacklist.csv as SE_trinity and
            // rerun to route them through SRA_FETCH_SE.
            SRA_FETCH.out.se_candidates
                .collectFile(name: 'rnaseq_se_candidates.csv', storeDir: launchDir, newLine: false)

            // Accessions whose download failed outright (parallel-fastq-dump + EBI FTP both
            // produced nothing) are recorded in rnaseq_blacklist.csv column order so they can be
            // reviewed and pasted straight into the blacklist as skip entries; merged at root.
            SRA_FETCH.out.blacklist_candidates
                .collectFile(name: 'rnaseq_blacklist_candidates.csv', storeDir: launchDir, newLine: false)
            reads_ch = SRA_FETCH.out.reads
                .mix(SRA_FETCH_SE.out.reads)
                .mix(WRITE_EMPTY_READS.out.reads)

            if (!params.stop_after_sra_fetch.toBoolean()) {
            // Build per-assembly channel keyed by species_tag with SRA reads joined.
            // reads_ch is now a 4-tuple: (species_tag, r1, r2, se)
            def assembly_with_reads = predict_genome_ch
                .map { meta, genome_fa ->
                    def species_tag = meta.species.replaceAll(/\s+/, '_')
                    tuple(species_tag, meta, genome_fa)
                }
                .combine(reads_ch, by: 0)
            // assembly_with_reads: (species_tag, meta, genome_fa, r1, r2, se)

            // RNASEQ_PREPARE: run funannotate train --stop_after_trinity once per species on
            // the representative (first) assembly, then cache the Trinity-GG FASTA in rnaseq_data/
            // so all other strains share it. Normalized reads stay in rnaseq_reads/ (SRA_FETCH storeDir).
            // pasa.gff3 is NOT produced here (--stop_after_trinity stops before PASA);
            // it is produced by FUNANNOTATE_TRAIN for every strain including the representative.
            // Species whose representative r1 and se are both zero-length skip RNASEQ_PREPARE
            // entirely; an empty trinity FASTA is written locally without submitting a SLURM job.
            def repr_ch = assembly_with_reads
                .groupTuple(by: 0)
                .map { species_tag, metas, genomes, r1s, r2s, ses ->
                    tuple(species_tag, metas[0], genomes[0], r1s[0], r2s[0], ses[0])
                }

            def repr_branched = repr_ch.branch {
                has_reads: it[3].size() > 0 || it[5].size() > 0   // r1=[3] or se=[5]
                no_reads:  true
            }

            RNASEQ_PREPARE(repr_branched.has_reads)

            // For species with no RNA-seq reads, write an empty trinity FASTA to rnaseq_data/
            // in the driver process (no SLURM job) and emit it directly as a shared channel item.
            def empty_shared_ch = repr_branched.no_reads
                .map { species_tag, _meta, _gfa, _r1, _r2, _se ->
                    def empty_fa = file("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
                    if (!empty_fa.exists()) {
                        empty_fa.parent.mkdirs()
                        empty_fa.text = ''
                    }
                    tuple(species_tag, empty_fa)
                }

            def shared_ch = RNASEQ_PREPARE.out.shared.mix(empty_shared_ch)

            // Join shared Trinity from rnaseq_data back to every assembly for FUNANNOTATE_TRAIN.
            // Normalized reads (r1/r2/se) come from SRA_FETCH/SRA_FETCH_SE via assembly_with_reads.
            def train_input = assembly_with_reads
                .combine(shared_ch, by: 0)
                .map { species_tag, meta, genome_fa, r1, r2, se, trinity_fa ->
                    tuple(meta, genome_fa, r1, r2, se, trinity_fa)
                }
            // train_input: meta=0, genome_fa=1, r1=2, r2=3, se=4, trinity_fa=5

            // Branch on r1 (idx 2), se (idx 4), or trinity_fa (idx 5) sizes.
            // Assemblies with no RNA-seq bypass FUNANNOTATE_TRAIN entirely.
            def branched = train_input.branch {
                has_rnaseq: it[2].size() > 0 || it[4].size() > 0 || it[5].size() > 0
                no_rnaseq:  true
            }
            def predict_no_rnaseq = branched.no_rnaseq
                .map { meta, genome_fa, _r1, _r2, _se, _tf ->
                    tuple(meta, genome_fa)
                }

            // Skip TRAIN at the channel level when pasa.gff3 already exists and is non-empty,
            // UNLESS the rnaseq reads or trinity FASTA is newer than the existing prediction GBK
            // (staleRnaseq), in which case we re-run training so predict can be refreshed too.
            def train_todo = branched.has_rnaseq.filter { meta, _gfa, _r1, _r2, _se, _tf ->
                def gff3 = file("${params.training_target}/${meta.id}/training/funannotate_train.pasa.gff3")
                !gff3.exists() || gff3.size() == 0 || staleRnaseq(meta.id as String, meta.species as String)
            }
            def train_done = branched.has_rnaseq
                .filter { meta, _gfa, _r1, _r2, _se, _tf ->
                    def gff3 = file("${params.training_target}/${meta.id}/training/funannotate_train.pasa.gff3")
                    gff3.exists() && gff3.size() > 0 && !staleRnaseq(meta.id as String, meta.species as String)
                }
                .map { meta, genome_fa, _r1, _r2, _se, _tf ->
                    tuple(meta, genome_fa)
                }
            FUNANNOTATE_TRAIN(train_todo)
            predict_input_ch = FUNANNOTATE_TRAIN.out.mix(train_done).mix(predict_no_rnaseq)
            } // end if (!params.stop_after_sra_fetch)
            } // end if (!params.stop_after_sra_query)
        } else {
            predict_input_ch = predict_genome_ch
        }

        if ((!params.stop_after_sra_fetch.toBoolean() && !params.stop_after_sra_query.toBoolean()) || !params.run_sra_fetch.toBoolean()) {
        def predict_ch = predict_input_ch
            .filter { meta, _gfa ->
                gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) == null || staleRnaseq(meta.id as String, meta.species as String)
            }
        FUNANNOTATE_PREDICT(predict_ch)

        // ── Post-predict steps and annotation ────────────────────────────────────
        // postpredict: all samples with a completed predict_results/*.gbk, whether
        // produced in this run or a prior one. This is the source for all optional
        // pre-annotate steps and for FUNANNOTATE_ANNOTATE itself.
        // INPUT_CHECK.out.samples already has taxon/asmid/suppress/n_test filters applied;
        // we just add the predict-results existence check on top.
        def postpredict = INPUT_CHECK.out.samples
            // Only genomes whose prediction was already complete AND current in a PRIOR run.
            // This is the exact logical complement of the predict_ch filter, so this set is
            // disjoint from the genomes (re)predicted in THIS run (which arrive via
            // FUNANNOTATE_PREDICT.out.metadata below). Keeping them disjoint means no genome
            // is fed downstream twice and stale genomes correctly wait for the fresh predict.
            .filter { meta ->
                gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) != null && !staleRnaseq(meta.id as String, meta.species as String)
            }

        // annotate_ready_ch threads through optional pre-annotate steps. Each optional
        // step splits the channel into "needs to run" vs "already done", processes the
        // former, then mixes the freshly-completed items back. FUNANNOTATE_ANNOTATE only
        // fires once all requested optional steps are complete for a given sample.
        // Joining ANTISMASH/INTERPRO/SIGNALP output back through predict_meta reconstructs
        // the metadata tuple while encoding the dependency edge in the channel DAG.
        //
        // Same-run completion gate: genomes predicted in THIS run flow in via
        // FUNANNOTATE_PREDICT.out.metadata (a real channel edge, so downstream waits for
        // predict to finish), while prior-run genomes flow in via postpredict (available
        // immediately). The two sets are disjoint by the filters above, so a plain mix
        // needs no dedup. (The optional steps below are still each gated behind their
        // params — run_antismash/interpro/signalp/update/annotate, all default false.)
        def predict_meta = postpredict.mix(FUNANNOTATE_PREDICT.out.metadata)
        def annotate_ready_ch = predict_meta

        if (params.run_antismash.toBoolean()) {
            def as_todo = annotate_ready_ch.filter { meta ->
                def asDir = file("${params.target}/${meta.id}/antismash_local")
                !(asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') })
            }
            def as_done = annotate_ready_ch.filter { meta ->
                def asDir = file("${params.target}/${meta.id}/antismash_local")
                asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') }
            }
            ANTISMASH_RUN(as_todo)
            ch_versions = ch_versions.mix(ANTISMASH_RUN.out.versions)
            def as_completed = ANTISMASH_RUN.out.results
                .map { meta, _files -> meta }
            annotate_ready_ch = as_completed.mix(as_done)
        }

        if (params.run_interpro.toBoolean()) {
            def ipr_todo = annotate_ready_ch.filter { meta ->
                !file("${params.target}/${meta.id}/annotate_misc/iprscan.xml").exists()
            }
            def ipr_done = annotate_ready_ch.filter { meta ->
                file("${params.target}/${meta.id}/annotate_misc/iprscan.xml").exists()
            }
            INTERPROSCAN_RUN(ipr_todo)
            ch_versions = ch_versions.mix(INTERPROSCAN_RUN.out.versions)
            def ipr_completed = INTERPROSCAN_RUN.out.results
                .map { meta, _xml -> meta }
            annotate_ready_ch = ipr_completed.mix(ipr_done)
        }

        if (params.run_signalp.toBoolean()) {
            def sp_todo = annotate_ready_ch.filter { meta ->
                !file("${params.target}/${meta.id}/annotate_misc/signalp.results.txt").exists()
            }
            def sp_done = annotate_ready_ch.filter { meta ->
                file("${params.target}/${meta.id}/annotate_misc/signalp.results.txt").exists()
            }
            SIGNALP_RUN(sp_todo)
            ch_versions = ch_versions.mix(SIGNALP_RUN.out.versions)
            def sp_completed = SIGNALP_RUN.out.results
                .map { meta, _txt -> meta }
            annotate_ready_ch = sp_completed.mix(sp_done)
        }

        if (params.run_update.toBoolean()) {
            if (params.run_sra_fetch.toBoolean()) {
                // UPDATE runs from predict results in parallel with antismash/interpro/signalp.
                // Reads are joined from SRA_FETCH (storeDir-cached, so prior-run reads are reused).
                // The join on upd_signal gates annotate_ready_ch so ANNOTATE waits for UPDATE.
                def upd_input = predict_meta
                    .map { meta ->
                        def species_tag = meta.species.replaceAll(/\s+/, '_')
                        tuple(species_tag, meta)
                    }
                    .combine(reads_ch, by: 0)
                    .map { _st, meta, r1, r2 ->
                        tuple(meta, r1, r2)
                    }
                def upd_todo = upd_input.filter { meta, _r1, _r2 ->
                    gbkResult("${params.target}/${meta.id}/update_results", meta.id as String) == null
                }
                def upd_done_signal = upd_input
                    .filter { meta, _r1, _r2 ->
                        gbkResult("${params.target}/${meta.id}/update_results", meta.id as String) != null
                    }
                    .map { meta, _r1, _r2 -> tuple(meta.id, 'upd') }
                FUNANNOTATE_UPDATE(upd_todo)
                def upd_signal = FUNANNOTATE_UPDATE.out
                    .map { meta -> tuple(meta.id, 'upd') }
                    .mix(upd_done_signal)
                annotate_ready_ch = annotate_ready_ch
                    .map { meta -> tuple(meta.id, meta) }
                    .join(upd_signal)
                    .map { _id, meta, _flag -> meta }
            } else {
                log.warn "run_update=true but run_sra_fetch=false; funannotate update skipped (no reads available)"
            }
        }

        if (params.run_annotate.toBoolean()) {
            FUNANNOTATE_ANNOTATE(annotate_ready_ch.filter { meta ->
                gbkResult("${params.target}/${meta.id}/annotate_results", meta.id as String) == null
            })
        }
        } // end if (!params.stop_after_sra_fetch || !params.run_sra_fetch)
    }

    // Collect software versions from all processes that emit versions.yml.
    // Written to logs/software_versions.yml alongside the trace file.
    ch_versions
        .unique()
        .collectFile(
            name:     'software_versions.yml',
            storeDir: "${launchDir}/logs/nextflow",
            newLine:  true
        )
}

