#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samples      = "${launchDir}/samples.csv"
params.target       = "${launchDir}/annotate"
params.inputfolder  = "predict_results"
params.taxon        = "fungi"

process ANTISMASH_RUN {
    tag "$out"

    cpus   8
    memory '16 GB'
    time   '60h'

    publishDir "${params.target}/${out}/antismash_local", mode: 'copy', overwrite: true

    input:
    tuple val(out), path(gbk)

    output:
    tuple val(out), path("antismash_local/**")

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load antismash
    antismash --taxon ${params.taxon} \\
        --output-dir antismash_local \\
        --genefinding-tool none \\
        --fullhmmer --clusterhmmer --cb-general --pfam2go \\
        -c ${task.cpus} \\
        ${gbk}
    pigz antismash_local/*.json
    """

    stub:
    """
    mkdir -p antismash_local
    touch antismash_local/${out}.json
    touch antismash_local/index.html
    """
}

workflow {
    def target = file(params.target)

    channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def species = row.SPECIES?.trim()?.replaceAll(/['"]/, '')
            def strain  = row.STRAIN?.trim()?.replaceAll(/['"]/, '')
            strain = strain?.replaceAll(/;.*$/, '')?.trim()
            def out = [species, strain].findAll { it }.join('_').replaceAll(/\s+/, '_')
            [out, row]
        }
        .filter { out, _row -> out }
        .map { out, _row ->
            def gbks   = file("${target}/${out}/${params.inputfolder}/*.gbk")
            def jsons  = file("${target}/${out}/antismash_local/*.json")
            def gzs    = file("${target}/${out}/antismash_local/*.json.gz")
            [out, gbks, jsons + gzs]
        }
        .filter { out, gbks, existing ->
            if (!gbks) {
                log.warn "No gbk for ${out} in ${target}/${out}/${params.inputfolder}"
                return false
            }
            if (existing) {
                log.info "Skipping ${out}: antismash json already present"
                return false
            }
            return true
        }
        .map { out, gbks, _existing -> tuple(out, gbks[0]) }
        .set { jobs }

    ANTISMASH_RUN(jobs)
}
