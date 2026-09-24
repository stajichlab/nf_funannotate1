// Single-species NCBI SRA query. storeDir caches the per-species CSV so re-runs
// skip the network query. To invalidate, delete rnaseq_reads/sra_query/<tag>.sra_query.csv.
// The batched variant (SRA_QUERY_BATCH) is preferred for many species at once.
//
// The runinfo fetch is capped to the first 250 UIDs (efetch -start 1 -stop
// 250 -- esearch itself has no -retmax flag) and given a 300s timeout: an
// uncapped query for a heavily-sequenced species (900+ matching SRA records)
// can take longer than a short timeout to fetch runinfo for, and would
// otherwise fail outright. esearch/sra also has no date -sort option, so
// "most recent first" is applied after the fetch, over whatever subset of
// UIDs esearch's default (not necessarily chronological) ordering returned
// in that first-250 window. See SRA_QUERY_BATCH/main.nf for the fuller
// rationale (both modules were fixed together, 2026-08-30).
//
// A best-effort keyword filter on LibraryName/SampleName excludes likely
// host-associated/co-infection or predation-model samples (mouse, macrophage,
// in vivo, blood, tissue, Acanthamoeba, Galleria, zebrafish, C. elegans,
// etc.) so Trinity isn't assembled against a mix of fungal + host/predator
// reads -- e.g. Acanthamoeba castellanii co-culture is a standard
// Cryptococcus virulence assay and would otherwise slip through untagged.
// This is a heuristic, not a guarantee -- SRA runinfo has no dedicated "pure
// culture vs host-associated" field; treat it as a first pass, not a
// substitute for a curated accession list on heavily-studied species.
//
// Community samples are excluded too: LibrarySource (runinfo col 15)
// METATRANSCRIPTOMIC / METAGENOMIC, rumen/microbiome/fecal keywords, and an
// explicit study denylist. Confirmed 2026-09-24: two Sichuan Agricultural
// University goat-rumen metatranscriptome studies -- SRP500439/PRJNA1097788 and
// SRP606484/PRJNA1301529 -- file one run per microbe taxon ("Microbe sample
// from <species>", SampleName like H.RUEMN.5.1) under 12 taxa each, 6 of them
// fungi (Aspergillus calidoustus, Batrachochytrium dendrobatidis, Rhizophagus
// irregularis, Rhizopus microsporus, Rozella allomycis, Spizellomyces
// punctatus). Their recent ReleaseDate ranked them first: R. microsporus got
// 100% rumen reads (26 of 49.2M hisat2-mapped, Trinity-GG 0 transcripts) and
// B. dendrobatidis 86%.
process SRA_QUERY {
    label 'edirect'
    label 'process_single'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads/sra_query"

    input:
    tuple val(species_tag), val(taxonid)

    output:
    tuple val(species_tag), path("${species_tag}.sra_query.csv"), emit: query_result

    script:
    """
    set -euo pipefail

    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${species_tag}.sra_query.csv

    timeout 300 bash -c "esearch -db sra \\
        -query 'txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND (BGISEQ[Platform] OR Illumina[Platform])' | \\
        efetch -format runinfo -start 1 -stop 250" > _runinfo.tmp

    # col 1=Run, col 2=ReleaseDate, col 4=spots, col 12=LibraryName, col 13=LibraryStrategy,
    # col 16=LibraryLayout, col 19=Platform, col 30=SampleName.
    # Prepend a platform rank (0=Illumina, 1=BGI/other) so the top 5 prefer Illumina,
    # then most-recent ReleaseDate first, then spot count desc as a final tiebreaker.
    # Best-effort keyword filter on LibraryName/SampleName excludes likely
    # host-associated/co-infection samples -- see module header comment.
    awk -F',' -v BL="${launchDir}/rnaseq_blacklist.csv" '
        BEGIN {
            # rnaseq_blacklist.csv 'skip' rows (action col 4; col 5 tolerated for rows that
            # carry spots in col 4). Dropped here so the top-N backfills with a clean run.
            while ((getline l < BL) > 0) { if (l ~ /^#/) continue; split(l, f, ","); if (f[4] == "skip" || f[5] == "skip") skip[f[1]] = 1 }
            close(BL)
        }
        NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ && \$4+0>=250000 {
            meta = " " tolower(\$12 " " \$30) " "
            if (meta ~ /(mouse|murine| mice | rat |rabbit|macrophage|phagocyt|in.?vivo|infect|co.?infec|co.?cultur|amoeba|acanthamoeba|galleria|zebrafish|elegans| host |blood|serum|plasma|csf|cerebrospinal|lung|brain|spleen|kidney| liver |tissue|biopsy|patient|clinical|autopsy|necropsy|bronch|rumen|ruemn|rumin|microbiom|metagenom|metatranscript|faec|fecal|feces|sludge|compost)/) next
            # LibrarySource (col 15) METATRANSCRIPTOMIC/METAGENOMIC = community sample, not this species.
            if (\$15 ~ /^META/) next
            # Known mislabeled studies (SRAStudy col 21 / BioProject col 22) -- see header comment.
            if (\$21 ~ /^(SRP500439|SRP606484)\$/ || \$22 ~ /^(PRJNA1097788|PRJNA1301529)\$/) next
            if (\$1 in skip) next
            # Same submitter as those studies (CenterName; runinfo has no author field): keep, but warn.
            if (toupper(\$0) ~ /SICHUAN AGRICULTURAL UNIVERSITY/) print "[WARN] SRA_QUERY: " \$1 " (" \$21 ") is from SICHUAN AGRICULTURAL UNIVERSITY -- check it is not a rumen/community sample" > "/dev/stderr"
            rank = (\$19 ~ /[Ii]llumina/) ? 0 : 1
            printf "%d,%s,%s,%s,%s\\n", rank, \$2, \$1, \$4, \$19
        }' _runinfo.tmp | \\
        sort -t',' -k1,1n -k2,2r -k4,4rn | \\
        head -n 5 | \\
        while IFS=',' read -r rank reldate acc spots platform; do
            printf '%s,%s,%s,%s,%s,PAIRED\\n' "${species_tag}" "${taxonid}" "\$acc" "\$spots" "\$platform"
        done >> ${species_tag}.sra_query.csv

    rm -f _runinfo.tmp
    NHITS=\$(awk 'END{print NR-1}' ${species_tag}.sra_query.csv)
    echo "[INFO] Found \$NHITS paired-end SRA accessions for ${species_tag} (taxonid=${taxonid})"

    # SE fallback: if no PE hits found and enable_single_end is true, query SINGLE layout
    if [ "${params.enable_single_end}" = "true" ] && [ "\$NHITS" -eq 0 ]; then
        timeout 300 bash -c "esearch -db sra \\
            -query 'txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND SINGLE[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND Illumina[Platform]' | \\
            efetch -format runinfo -start 1 -stop 250" > _runinfo_se.tmp
        awk -F',' -v BL="${launchDir}/rnaseq_blacklist.csv" '
            BEGIN {
                # rnaseq_blacklist.csv 'skip' rows (action col 4; col 5 tolerated for rows that
                # carry spots in col 4). Dropped here so the top-N backfills with a clean run.
                while ((getline l < BL) > 0) { if (l ~ /^#/) continue; split(l, f, ","); if (f[4] == "skip" || f[5] == "skip") skip[f[1]] = 1 }
                close(BL)
            }
            NR>1 && \$13=="RNA-Seq" && \$16=="SINGLE" && \$1~/^[SDE]RR/ && \$4+0>=250000 {
                meta = " " tolower(\$12 " " \$30) " "
                if (meta ~ /(mouse|murine| mice | rat |rabbit|macrophage|phagocyt|in.?vivo|infect|co.?infec|co.?cultur|amoeba|acanthamoeba|galleria|zebrafish|elegans| host |blood|serum|plasma|csf|cerebrospinal|lung|brain|spleen|kidney| liver |tissue|biopsy|patient|clinical|autopsy|necropsy|bronch|rumen|ruemn|rumin|microbiom|metagenom|metatranscript|faec|fecal|feces|sludge|compost)/) next
                # LibrarySource (col 15) METATRANSCRIPTOMIC/METAGENOMIC = community sample, not this species.
                if (\$15 ~ /^META/) next
                # Known mislabeled studies (SRAStudy col 21 / BioProject col 22) -- see header comment.
                if (\$21 ~ /^(SRP500439|SRP606484)\$/ || \$22 ~ /^(PRJNA1097788|PRJNA1301529)\$/) next
                if (\$1 in skip) next
                # Same submitter as those studies (CenterName; runinfo has no author field): keep, but warn.
                if (toupper(\$0) ~ /SICHUAN AGRICULTURAL UNIVERSITY/) print "[WARN] SRA_QUERY: " \$1 " (" \$21 ") is from SICHUAN AGRICULTURAL UNIVERSITY -- check it is not a rumen/community sample" > "/dev/stderr"
                printf "%s,%s,%s,%s\\n", \$2, \$1, \$4, \$19
            }' _runinfo_se.tmp | \\
            sort -t',' -k1,1r -k3,3rn | \\
            head -n ${params.max_rnaseq_se_runs} | \\
            while IFS=',' read -r reldate acc spots platform; do
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
