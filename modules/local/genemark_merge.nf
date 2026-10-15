// Gather step of the chunked GeneMark route: merge the chunk GTFs of one genome in chunk order, renumber gene and
// transcript IDs, validate against the whole genome, and emit gtf + mod with the same names and keyed-tuple shapes
// as the monolithic route (genemark_sidecar.nf publishes them; this process has no publishDir, like the old
// GENEMARK_RUN). The model is emitted here, after validation, so a half-merged result never produces one.
// The process shell has no `set -e` (nextflow.config process.shell), so strict mode is set below.
process GENEMARK_MERGE {
    label 'genemark'
    tag "${meta.id}"

    input:
    tuple val(meta), val(n_chunks), val(genome_fa), path(chunk_gtfs, stageAs: 'gtfs/*'), path(model)

    output:
    tuple val(meta), path("${meta.id}.genemark.gtf"), emit: gtf
    tuple val(meta), path("${meta.id}.genemark.mod"), emit: mod

    script:
    def out = meta.id
    """
    set -euo pipefail
    # 1. exactly the expected chunks, numbered 0..n-1 with no gap and no duplicate
    EXPECTED=${n_chunks}
    GOT=\$(ls gtfs/*.gtf | wc -l)
    [ "\$GOT" -eq "\$EXPECTED" ] || { echo "ERROR: GENEMARK_MERGE ${out}: got \$GOT chunk GTFs, expected \$EXPECTED" >&2; exit 1; }
    for i in \$(seq 0 \$((EXPECTED - 1))); do
        [ -e "\$(printf 'gtfs/${out}.chunk_%04d.gtf' \$i)" ] || { echo "ERROR: GENEMARK_MERGE ${out}: chunk \$i is missing" >&2; exit 1; }
    done
    GTFS=\$(for i in \$(seq 0 \$((EXPECTED - 1))); do printf 'gtfs/${out}.chunk_%04d.gtf ' \$i; done)

    python "${workflow.projectDir}/bin/genemark_merge_gtf.py" --out "${out}.genemark.gtf" \$GTFS
    [ -s "${out}.genemark.gtf" ] || { echo "ERROR: GENEMARK_MERGE ${out}: merged GTF is empty" >&2; exit 1; }

    # 2. every GTF contig name must be a contig of the whole genome; an empty contig list is a failure
    case "${genome_fa}" in
        *.gz) gzip -dc "${genome_fa}" > genome_for_names.fa ;;
        *)    cp "${genome_fa}" genome_for_names.fa ;;
    esac
    awk '/^>/ { print substr(\$1, 2) }' genome_for_names.fa > contig_names.txt
    [ -s contig_names.txt ] || { echo "ERROR: GENEMARK_MERGE ${out}: no contig names read from ${genome_fa}" >&2; exit 1; }
    if ! awk 'NR == FNR { ids[\$1] = 1; next }
              !/^#/ && NF >= 9 && !(\$1 in ids) { print "unknown sequence " \$1 > "/dev/stderr"; bad = 1; exit 1 }
              END { exit bad }' contig_names.txt "${out}.genemark.gtf"; then
        echo "ERROR: GENEMARK_MERGE ${out}: merged GTF names sequences that are not in the genome" >&2
        exit 1
    fi
    rm -f genome_for_names.fa contig_names.txt

    [ -s "${model}" ] || { echo "ERROR: GENEMARK_MERGE ${out}: model is empty" >&2; exit 1; }
    cp "${model}" "${out}.genemark.mod"
    # same predict_misc persistence as the monolithic route (see genemark_train.nf)
    LOWER_OUT=\$(printf '%s' "${out}" | tr '[:upper:]' '[:lower:]')
    AB_INITIO_DIR="${params.target}/${out}/predict_misc/ab_initio_parameters"
    mkdir -p "\$AB_INITIO_DIR"
    cp "${out}.genemark.mod" "\$AB_INITIO_DIR/\${LOWER_OUT}.genemark.mod"
    echo "[INFO] GENEMARK_MERGE complete for ${out}"
    """

    stub:
    """
    touch "${meta.id}.genemark.gtf" "${meta.id}.genemark.mod"
    """
}
