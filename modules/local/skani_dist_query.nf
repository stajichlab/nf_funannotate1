/*
 * SKANI_DIST_QUERY — asymmetric query-vs-reference skani ANI (orphan placement).
 *
 * Port of BFD/Fungi_BFD/nextflow/modules/ani/compare/SKANI_DIST_QUERY/main.nf.
 * skani 0.3.x: sketch query and reference genomes into separate directories,
 * then `skani dist -q <query_sketch> -rl <ref_list>`; keeps the top
 * `params.query_top_n` ref hits per query. Output 3-column query/ref/ANI TSV.
 */
process SKANI_DIST_QUERY {
    tag   "${group_name} q=${query_genomes.size()} r=${ref_genomes.size()}"
    label 'skani'

    cpus   { (ref_genomes.size() + query_genomes.size()) > 2000 ? 32 : 8 }
    memory {
        def baseGB  = (ref_genomes.size() + query_genomes.size()) > 2000 ? 64 : 16
        def floorGB = task.attempt >= 3 ? 100 : task.attempt == 2 ? 24 : 6
        Math.max(baseGB, floorGB) + ' GB'
    }
    time {
        def n = ref_genomes.size() + query_genomes.size()
        def baseH  = n > 500 ? 6 : n > 200 ? 3 : 1
        def floorH = task.attempt >= 2 ? 6 : 1
        Math.min(6, Math.max(baseH, floorH)) + ' h'
    }
    stageInMode('symlink')

    publishDir { "${params.ani_outdir}/skani_query/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
    tuple val(group_name), path(query_genomes), path(ref_genomes)

    output:
    tuple val(group_name), path("${group_name}.query.ani.tsv")

    script:
    def q_list = query_genomes.collect { it.toString() }.join('\n')
    def r_list = ref_genomes.collect   { it.toString() }.join('\n')
    """
    printf '%s\\n' "${q_list}" > query_list.txt
    printf '%s\\n' "${r_list}" > ref_list.txt
    skani sketch -t ${task.cpus} -l query_list.txt -o query_sketch
    skani sketch -t ${task.cpus} -l ref_list.txt   -o ref_sketch

    skani dist -q query_sketch -rl ref_list.txt \\
        --min-af ${params.skani_min_af} -n ${params.query_top_n} -t ${task.cpus} \\
        -o skani_raw.tsv

    awk -F'\\t' 'NR>1 { print \$1"\\t"\$2"\\t"\$3 }' skani_raw.tsv > ${group_name}.query.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.query.ani.tsv
    """
}
