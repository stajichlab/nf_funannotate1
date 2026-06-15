#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samples         = "${launchDir}/samples.csv"
params.target          = "${launchDir}/annotate"
params.proteins	       = "${launchDir}/lib/swissprot_fungi.faa"
params.source          = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM"
params.seqcenter       = "NCBI"
params.augustus_config = "${launchDir}/lib/augustus/3.5/config"
params.funannotate_db  = "/bigdata/stajichlab/shared/lib/funannotate_db"
params.min_contig_len  = 2000
params.clean_script    = "${launchDir}/scripts/clean_genome_fa.py"
params.sbt_template    = "${launchDir}/lib/template.sbt"  // fill in correct path
params.debug           = false   // --debug: verbose logging in script + channel views
params.n_test          = 0       // --n_test N: limit to first N samples (0 = all)

process FUNANNOTATE_PREDICT {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '32h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag), val(busco_lineage), val(header_length), val(transl_table), path(genome_gz)

    output:
    tuple val(out), path("${out}/**")

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    # ---- debug block -------------------------------------------------------
    echo "[DEBUG] out          = ${out}"
    echo "[DEBUG] asmid        = ${asmid}"
    echo "[DEBUG] species      = ${species}"
    echo "[DEBUG] strain       = ${strain}"
    echo "[DEBUG] locustag     = ${locustag}"
    echo "[DEBUG] busco        = ${busco_lineage}"
    echo "[DEBUG] transl_table = ${transl_table}"
    echo "[DEBUG] proteins     = ${params.proteins}"
    echo "[DEBUG] genome_gz NF = ${genome_gz}"
    echo "[DEBUG] TMPDIR       = \$TMPDIR"
    echo "[DEBUG] pwd          = \$(pwd)"
    ls -lah .
    stat ${genome_gz} 2>&1 || echo "[DEBUG] stat failed for ${genome_gz}"
    echo "[DEBUG] is symlink: \$([ -L ${genome_gz} ] && readlink -f ${genome_gz} || echo 'not a symlink')"
    echo "[DEBUG] pigz version: \$(pigz --version 2>&1)"
    # ---- end debug block ---------------------------------------------------

    GENOME=\$TMPDIR/${asmid}.fa

    if [ ! -f "${genome_gz}" ]; then
        echo "ERROR: genome_gz not found at path: ${genome_gz}" >&2
        exit 1
    fi

    echo "[INFO] Decompressing and cleaning genome..."
    pigz -dc ${genome_gz} | ${params.clean_script} --len ${params.min_contig_len} > \$GENOME
    echo "[INFO] Genome written to \$GENOME (size: \$(du -sh \$GENOME | cut -f1))"

    TBL2ASN_PARAMS="-l paired-ends"

    funannotate predict --name ${locustag} -i \$GENOME --strain "${strain}" \\
        -o ${out} -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w codingquarry:0 \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length 24 --protein_evidence ${params.proteins} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table}

    F=\$(ls ${out}/predict_results/*.gbk 2>/dev/null | head -n 1)
    if [ -z "\$F" ]; then
        echo "ERROR: funannotate predict did not produce a .gbk in ${out}/predict_results — treating as failure" >&2
        exit 1
    fi
    mv ${out}/predict_misc/ab_initio_parameters ${out}
    rm -rf ${out}/predict_misc
    mkdir -p ${out}/predict_misc
    mv ${out}/ab_initio_parameters ${out}/predict_misc
    pigz ${out}/predict_results/*.txt ${out}/predict_results/*.mrna-transcripts.fa

    rm -f \$GENOME
    """

    stub:
    """
    echo "[STUB] Would run funannotate predict for ${out} using ${genome_gz}"
    mkdir -p ${out}/predict_results
    touch ${out}/predict_results/${out}.gbk
    """
}

// TODO: interpro scan with nextflow sub?
// TODO: signalp on gpu?

process FUNANNOTATE_ANNOTATE {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '48h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag), val(busco_lineage), val(header_length)

    output:
    tuple val(out), path("${out}/**")

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    # ---- debug block -------------------------------------------------------
    echo "[DEBUG] out      = ${out}"
    echo "[DEBUG] asmid    = ${asmid}"
    echo "[DEBUG] locustag = ${locustag}"
    echo "[DEBUG] busco    = ${busco_lineage}"
    echo "[DEBUG] pwd      = \$(pwd)"
    ls -lah .
    # ---- end debug block ---------------------------------------------------

    funannotate annotate -i ${params.target}/${out} \\
        --species "${species}" --strain "${strain}" \\
        --busco_db ${busco_lineage} --rename ${locustag} \\
        --sbt ${params.sbt_template} \\
        --cpu ${task.cpus} --tmpdir \$TMPDIR
    """

    stub:
    """
    echo "[STUB] Would run funannotate annotate for ${out}"
    mkdir -p ${out}/annotate_results
    touch ${out}/annotate_results/${out}.gbk
    """
}

// Check whether predict_results already has a gbk for this sample.
// Uses explicit directory listing rather than a glob path object, which
// is always truthy in Groovy regardless of whether files exist.
def hasExistingGbk(targetDir, out) {
    def dir = new File("${targetDir}/${out}/predict_results")
    if (!dir.exists()) return false
    return dir.list()?.any { f -> f.endsWith('.gbk') } ?: false
}

workflow {
    def target = file(params.target)

    def jobs = channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def species      = row.SPECIES?.trim()?.replaceAll(/['"]/, '')
            def strain       = row.STRAIN?.trim()?.replaceAll(/['"]/, '')
            strain = strain.replaceAll(/;.*$/,'').trim()
            def out          = [species,strain].findAll { it }.join('_').replaceAll(/\s+/, '_')
            def asmid        = row.ASMID?.trim()
            def locustag     = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def busco        = row.BUSCO_LINEAGE?.trim()
            def header_length = 24
            def transl_table = row.TRANSL_TABLE?.trim() ?: '1'
            [out, asmid, species, strain, locustag, busco, header_length, transl_table]
        }
        .filter { out, asmid, _species, _strain, _locustag, _busco, _header_length, _transl_table ->
            out && asmid
        }
        // n_test > 0 limits to first N samples; -1 means take all
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        .map { out, asmid, species, strain, locustag, busco, header_length, transl_table ->
            def gz = file("${params.source}/${asmid}/${asmid}_genomic.fna.gz")
            [out, asmid, species, strain, locustag, busco, header_length, transl_table, gz]
        }
        .filter { out, asmid, _species, _strain, _locustag, _busco, _header_length, _transl_table, gz ->
            if (hasExistingGbk(target, out)) {
                log.info "Skipping ${out}: predict_results gbk already present"
                return false
            }
            if (!gz.exists()) {
                log.warn "Missing genome for ${out} (asmid=${asmid}): ${gz}"
                return false
            }
            if (params.debug) {
                log.info "Queuing ${out}: genome=${gz} (${gz.size()} bytes)"
            }
            return true
        }
        .map { out, asmid, species, strain, locustag, busco, header_length, transl_table, gz ->
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, gz)
        }

    if (params.debug) {
        jobs.view { t -> "[CHANNEL] Submitting: out=${t[0]}, asmid=${t[1]}, transl_table=${t[7]}, gz=${t[8]}" }
    }

    FUNANNOTATE_PREDICT(jobs)
}
