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
include { MASK_GENOME }      from './subworkflows/local/mask_genome'
include { FETCH_RNASEQ }     from './subworkflows/local/fetch_rnaseq'
include { TRAIN_PREDICT }    from './subworkflows/local/train_predict'
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
        // predict_genome_ch: tantan soft-masked (default) or clean/prior-masked genome.
        // MASK_GENOME handles the run_repeatmasker if/else and storeDir-cached masking.
        MASK_GENOME(clean_genome_ch)
        def predict_genome_ch = MASK_GENOME.out.genomes

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

        // SRA read fetching (FETCH_RNASEQ) + training (RNASEQ_PREPARE + FUNANNOTATE_TRAIN)
        // + prediction (FUNANNOTATE_PREDICT). All three are composed in TRAIN_PREDICT.
        // When run_sra_fetch=false, reads channel is empty and genomes go straight to predict.
        def reads_ch = Channel.empty()
        if (params.run_sra_fetch.toBoolean()) {
            def sra_input = predict_genome_ch
                .map { meta, _genome_fa ->
                    def species_tag = meta.species.replaceAll(/\s+/, '_')
                    tuple(species_tag, meta.taxonid)
                }
                .groupTuple(by: 0)
                .map { species_tag, taxonids -> tuple(species_tag, taxonids[0]) }

            FETCH_RNASEQ(sra_input)
            reads_ch = FETCH_RNASEQ.out.reads
        }

        if (!params.run_sra_fetch.toBoolean() || (!params.stop_after_sra_fetch.toBoolean() && !params.stop_after_sra_query.toBoolean())) {
        TRAIN_PREDICT(predict_genome_ch, reads_ch)

        // ── Post-predict steps and annotation ────────────────────────────────────
        // postpredict: all samples with a completed predict_results/*.gbk, whether
        // produced in this run or a prior one. This is the source for all optional
        // pre-annotate steps and for FUNANNOTATE_ANNOTATE itself.
        def postpredict = INPUT_CHECK.out.samples
            // Only genomes predicted in a PRIOR run — complement of what TRAIN_PREDICT runs now.
            .filter { meta ->
                gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) != null && !staleRnaseq(meta.id as String, meta.species as String)
            }

        // Same-run completion gate: new predictions arrive via TRAIN_PREDICT.out.metadata;
        // prior-run completions arrive via postpredict (available immediately).
        def predict_meta = postpredict.mix(TRAIN_PREDICT.out.metadata)
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

