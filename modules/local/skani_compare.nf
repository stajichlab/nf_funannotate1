/*
 * SKANI_COMPARE — symmetric all-vs-all skani ANI within one group.
 *
 * Port of BFD/Fungi_BFD/nextflow/modules/ani/compare/SKANI_COMPARE/main.nf.
 * skani 0.3.x: sketch everything into one directory, then `skani triangle -E`
 * on it; sparse mode outputs name1<TAB>name2<TAB>ANI<TAB>AF_ref<TAB>AF_query.
 * Self-comparisons (diagonal) are stripped, leaving a 3-column query/ref/ANI TSV.
 *
 * Resource closures inline BFD's skaniCpusFor/skaniMemoryFor/aniTimeFor (the
 * k8s caps in capCpus/capMemGB are no-ops here). `params.skani_preset` /
 * `params.skani_compression` / `params.skani_min_af` control the run.
 */
process SKANI_COMPARE {
    tag   "${group_name} [${n_genomes} genomes]"
    label 'skani'

    cpus   { n_genomes > 500 ? 32 : n_genomes > 200 ? 16 : 8 }
    memory {
        def baseGB  = n_genomes > 500 ? 48 : n_genomes > 200 ? 32 : 16
        def floorGB = task.attempt >= 3 ? 100 : task.attempt == 2 ? 24 : 6
        Math.max(baseGB, floorGB) + ' GB'
    }
    time {
        def baseH  = n_genomes > 500 ? 6 : n_genomes > 200 ? 3 : 1
        def floorH = task.attempt >= 2 ? 6 : 1
        Math.min(6, Math.max(baseH, floorH)) + ' h'
    }
    stageInMode('symlink')

    publishDir { "${params.ani_outdir}/${params.ani_method}/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
    tuple val(group_name), val(n_genomes), path(genomes)

    output:
    tuple val(group_name), path("${group_name}.full.ani.tsv")

    script:
    def presetFlags = ['fast': '--fast', 'medium': '', 'slow': '--slow']
    def presetKey   = (params.skani_preset ?: 'medium').toString().toLowerCase()
    def preset      = presetFlags.containsKey(presetKey) ? presetFlags[presetKey] : error("--skani_preset must be one of: fast, medium, slow")
    def cflag       = (params.skani_compression as int) > 0 ? "-c ${params.skani_compression}" : ''
    def genome_list = genomes.collect { it.toString() }.join('\n')
    """
    printf '%s\\n' "${genome_list}" > genome_list.txt
    skani sketch ${preset} ${cflag} -t ${task.cpus} -l genome_list.txt -o sketch_db

    skani triangle -l genome_list.txt -E \\
        --min-af ${params.skani_min_af} -t ${task.cpus} \\
        -o skani_raw.tsv

    awk -F'\\t' 'NR>1 && \$1 != \$2 {
        print \$1"\\t"\$2"\\t"\$3
    }' skani_raw.tsv > ${group_name}.full.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.full.ani.tsv
    """
}
