// Run InterProScan 6 on one assembly's proteome, producing
// annotate_misc/iprscan.xml.gz for funannotate annotate (which reads the .gz
// directly when no plain iprscan.xml is present) plus iprscan.tsv.gz.
//
// InterProScan 6 is not a single program: it is its own Nextflow pipeline
// (ebi-pf-team/interproscan6) that runs each member-database search in a
// separate container. So this process has NO container of its own. It
// launches IPS6 as a nested `nextflow run` on the host, with IPS6's local
// executor, inside this one task's allocation. One task per genome, not one
// scheduler job per IPS6 sub-task.
//
// Shared, read-only resources (set up once by scripts/setup_interproscan6.sh):
//   params.iprscan6_pipeline  pinned local checkout of the IPS6 workflow
//   params.sif_dir            IPS6's container images, pre-pulled
//   params.iprscan6_datadir   InterPro + member-DB data (<datadir>/<db>/<ver>/)
// The nested work dir lives on node-local $SCRATCH and is removed on success.
//
// Licensed IPS6 apps (DeepTMHMM, Phobius, SignalP) stay off unless
// params.iprscan6_licensed_config is set: funannotate reads only InterPro
// entries + GO terms from this XML, and DEEPTMHMM_ANNOTATION / SIGNALP_RUN
// already feed funannotate its TM / signal-peptide calls.
process INTERPROSCAN_RUN {
    label 'interproscan6'
    label 'process_medium'
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}.iprscan.done"), emit: results
    path 'versions.yml',                              emit: versions

    script:
    def out      = meta.id
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    // Option B persistence (like DEEPTMHMM_ANNOTATION): write straight to the
    // persistent target dir that annotate_genome.nf's done-check and
    // FUNANNOTATE_ANNOTATE both read.
    def miscdir  = "${params.target}/${out}/annotate_misc"
    def pipeline = params.iprscan6_pipeline as String
    // A local checkout runs as-is; a remote name (ebi-pf-team/interproscan6)
    // needs the pinned revision.
    def revArg   = pipeline.startsWith('/') ? '' : "-r ${params.iprscan6_revision}"
    def licArg   = params.iprscan6_licensed_config ? "-c ${params.iprscan6_licensed_config}" : ''
    def applArg  = params.iprscan6_applications ? "--applications ${params.iprscan6_applications}" : ''
    def apiArg   = params.iprscan6_matches_api.toString().toBoolean() ? '' : '--no-matches-api'
    // IPS6's local executor would otherwise assume the whole node's memory.
    // Split this task's allocation: 1/4 to the nested head JVM (it runs
    // IPS6's parse/output steps in-process), the rest to its sub-tasks.
    def memGb    = task.memory ? (task.memory.toGiga() as int) : 16
    def heapGb   = Math.max(2, memGb.intdiv(4))
    def execGb   = Math.max(2, memGb - heapGb)
    def subCpus  = Math.min(params.iprscan6_task_cpus as int, task.cpus as int)
    def ctOpts   = params.iprscan6_container_options ?: ''
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    if ! command -v ${params.iprscan6_profile} >/dev/null; then
        echo "ERROR: '${params.iprscan6_profile}' (params.iprscan6_profile) is not on PATH;" \\
             "the nested InterProScan 6 run needs it to start its containers" >&2
        exit 1
    fi
    IPR_DIR="${params.iprscan6_datadir}/interpro/${params.iprscan6_interpro}"
    if [ ! -d "\$IPR_DIR" ]; then
        echo "WARNING: \$IPR_DIR not found;" \\
             "IPS6 will download the data itself (see scripts/setup_interproscan6.sh)" >&2
    else
        # A shared datadir with unreadable files makes IPS6 fail with only a
        # bare path in FIND_DATABASES; name the real problem up front.
        UNREADABLE=\$(find "${params.iprscan6_datadir}" -maxdepth 3 ! -readable 2>/dev/null | head -5)
        UNREADABLE="\$UNREADABLE\$(find "\$IPR_DIR" ! -readable 2>/dev/null | head -5)"
        if [ -n "\$UNREADABLE" ]; then
            echo "ERROR: files under ${params.iprscan6_datadir} are not readable by \$(id -un):" >&2
            echo "\$UNREADABLE" >&2
            exit 1
        fi
    fi
    # Node-local scratch may not exist / be writable if an inherited \$SCRATCH
    # points at another node's path — fall back to the task workdir.
    TMPDIR=\$(printf '%s' "\${SCRATCH:-}" | tr -d '\\n\\r')
    TMPDIR=\${TMPDIR:-/tmp}
    if [ ! -d "\$TMPDIR" ] || [ ! -w "\$TMPDIR" ]; then
        TMPDIR="\$PWD"
    fi
    export TMPDIR
    IPS6_WORK=\$(mktemp -d "\$TMPDIR/ips6_${out}.XXXXXX")

    # The containers inherit TMPDIR (node-local, e.g. /scratch/<user>/<job>)
    # but do not mount it, so tools that write temp files there fail with
    # "Read-only file system" (seen: PROSITE ps_scan.pl). Bind it.
    cat > ips6_local.config <<-END_CONFIG
    executor { memory = '${execGb} GB' }
    apptainer {
        cacheDir   = '${params.sif_dir}'
        autoMounts = true
        runOptions = '${ctOpts} -B \$TMPDIR:\$TMPDIR'
    }
    singularity {
        cacheDir   = '${params.sif_dir}'
        autoMounts = true
        runOptions = '${ctOpts} -B \$TMPDIR:\$TMPDIR'
    }
    END_CONFIG

    # Isolate the nested run from the parent run's Nextflow environment.
    unset NXF_WORK NXF_PARAMS_FILE NXF_CLI NXF_CLI_OPTS
    export NXF_APPTAINER_CACHEDIR=${params.sif_dir}
    export NXF_SINGULARITY_CACHEDIR=${params.sif_dir}
    export NXF_ANSI_LOG=false
    export NXF_TEMP=\$TMPDIR
    # IPS6 copies its SQLite sequence DB to java.io.tmpdir (node-local).
    export NXF_OPTS="-Xms1g -Xmx${heapGb}g -Djava.io.tmpdir=\$TMPDIR"

    # IPS6 requires --outdir to exist already. It can also exit 0 after a
    # failed parameter check, so success is judged by its output files below.
    mkdir -p ips6_out
    ${params.iprscan6_nextflow} run ${pipeline} ${revArg} \\
        -profile ${params.iprscan6_profile} \\
        -c ips6_local.config ${licArg} \\
        -w \$IPS6_WORK \\
        --input ${proteins} \\
        --datadir ${params.iprscan6_datadir} \\
        --interpro ${params.iprscan6_interpro} \\
        --formats xml,tsv --goterms --pathways \\
        --cpus ${subCpus} --max-workers ${task.cpus} \\
        ${applArg} ${apiArg} ${params.iprscan6_extra_args} \\
        --outdir \$PWD/ips6_out --outprefix ${out}.iprscan

    for f in xml tsv; do
        if [ ! -s ips6_out/${out}.iprscan.\$f ]; then
            echo "ERROR: InterProScan 6 wrote no ${out}.iprscan.\$f" >&2
            exit 1
        fi
    done
    IPS6_VERSION=\$(grep -m1 -o 'interproscan-version="[^"]*"' ips6_out/${out}.iprscan.xml | cut -d'"' -f2)
    IPR_VERSION=\$(grep -m1 -o 'interpro-version="[^"]*"' ips6_out/${out}.iprscan.xml | cut -d'"' -f2)

    # .part + mv: a killed task never leaves a truncated file that the
    # done-check in annotate_genome.nf would accept.
    mkdir -p ${miscdir}
    gzip -c ips6_out/${out}.iprscan.xml > ${miscdir}/iprscan.xml.gz.part
    gzip -c ips6_out/${out}.iprscan.tsv > ${miscdir}/iprscan.tsv.gz.part
    mv ${miscdir}/iprscan.tsv.gz.part ${miscdir}/iprscan.tsv.gz
    mv ${miscdir}/iprscan.xml.gz.part ${miscdir}/iprscan.xml.gz
    rm -rf \$IPS6_WORK ips6_out
    touch ${out}.iprscan.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \${IPS6_VERSION}
        interpro_data: \${IPR_VERSION}
    END_VERSIONS
    """

    stub:
    def out     = meta.id
    def miscdir = "${params.target}/${out}/annotate_misc"
    """
    mkdir -p ${miscdir}
    echo '<results/>' | gzip -c > ${miscdir}/iprscan.xml.gz
    touch ${out}.iprscan.done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: stub
        interpro_data: stub
    END_VERSIONS
    """
}
