// Shared filesystem utilities for funannotate.nf, earlgrey_mask.nf, and their
// subworkflows. Uses java.io.File so this class can be loaded from lib/ without
// access to Nextflow's script-scope globals (params, launchDir, file()).
class FunannotateUtils {

    // Returns the first non-empty GBK file (.gbk preferred over .gbk.gz), or null.
    // Use for skip/completion gating so compressed result folders still count as done.
    static File gbkResult(String dir, String id) {
        def plain = new File("${dir}/${id}.gbk")
        if (plain.exists() && plain.size() > 0) return plain
        def gz = new File("${dir}/${id}.gbk.gz")
        if (gz.exists() && gz.size() > 0) return gz
        return null
    }

    // Returns the first non-empty genome file (.gz preferred), else the plain path.
    // Falls back gracefully so callers' .exists() checks still report missing files.
    static File genomeFile(String base) {
        def gz = new File("${base}.gz")
        if (gz.exists() && gz.size() > 0) return gz
        return new File(base)
    }

    // Returns true when any rnaseq/trinity input is newer than the existing predict GBK,
    // indicating that training and prediction need to be refreshed.
    // target and launchDir must be passed explicitly (not available in lib/ class scope).
    static boolean staleRnaseq(String id, String species, String target, String launchDir) {
        def species_tag = species.replaceAll(/\s+/, '_')
        def gbk = gbkResult("${target}/${id}/predict_results", id)
        if (gbk == null) return false
        def gbkMod = gbk.lastModified()
        def r1      = new File("${launchDir}/rnaseq_reads/${species_tag}_norm_R1.fastq.gz")
        def se      = new File("${launchDir}/rnaseq_reads/${species_tag}_norm_SE.fastq.gz")
        def trinity = new File("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
        return (r1.exists()      && r1.size() > 0      && r1.lastModified()      > gbkMod) ||
               (se.exists()      && se.size() > 0      && se.lastModified()      > gbkMod) ||
               (trinity.exists() && trinity.size() > 0 && trinity.lastModified() > gbkMod)
    }

    // Returns true when the genome's source assembly archive (as resolved by
    // input_check.nf's default params.source/<asmid>/<asmid>_genomic.fna.gz
    // convention) is newer than the existing predict GBK, indicating a
    // re-assembled/updated genome that needs retraining and repredicting --
    // otherwise a swapped assembly would silently keep a GBK annotated against
    // stale coordinates. Mirrors staleRnaseq's shape/purpose, keyed off the
    // source assembly file instead of the rnaseq inputs. Port of BFD
    // funannotate/utils.nf's staleGenome(). Caveat: only covers the
    // asmid-derived default path -- a per-row `GENOME=` override (see
    // input_check.nf) isn't visible here since this class has no row-level
    // context, so a locally-overridden genome path's own staleness isn't
    // detected by this check.
    static boolean staleGenome(String id, String asmid, String source, String target) {
        def gbk = gbkResult("${target}/${id}/predict_results", id)
        if (gbk == null) return false
        def gfa = new File("${source}/${asmid}/${asmid}_genomic.fna.gz")
        return gfa.exists() && gfa.size() > 0 && gfa.lastModified() > gbk.lastModified()
    }

    // ── Species-level ab-initio parameter reuse (port of BFD funannotate/utils.nf) ──
    // The shared store lives at <sharedRoot>/<species_tag>/ (writes:
    // parameters.json + <species_tag>.genemark.mod), a top-level sibling of
    // params.target shared across every project tree annotating the same species.
    // All paths are passed explicitly (not available in lib/ class scope).

    // Absolute path to a species' shared ab-initio parameters.json, or null if the
    // species has no backfilled store yet. No store on disk is the correct, safe
    // signal to fall back to independent training, not an error.
    static File sharedParamsJsonFor(String species, String sharedRoot) {
        def species_tag = species.replaceAll(/\s+/, '_')
        def json = new File("${sharedRoot}/${species_tag}/parameters.json")
        return (json.exists() && json.size() > 0) ? json : null
    }

    // Same existence/non-empty check as sharedParamsJsonFor(), pointed at the
    // species' shared GeneMark .mod instead of parameters.json. Used by GENEMARK_RUN
    // to decide fast-reuse (--predict_with) vs fresh --ES training.
    static File sharedGenemarkModFor(String species, String sharedRoot) {
        def species_tag = species.replaceAll(/\s+/, '_')
        def mod = new File("${sharedRoot}/${species_tag}/${species_tag}.genemark.mod")
        return (mod.exists() && mod.size() > 0) ? mod : null
    }

    // A strain's FUNANNOTATE_TRAIN-produced transcript-to-genome alignment BAM --
    // GENEMARK_RUN's ET mode derives RNA-seq-informed intron hints from this.
    // Returns '' (not null) when absent/empty: GENEMARK_RUN's script checks for a
    // non-empty string, and an absent file is the safe signal to fall back to ES.
    // Returns true when the genome's TRAINING evidence is newer than its existing
    // predict GBK -- i.e. it has been re-trained since it was last annotated, so
    // the published annotation no longer reflects the evidence and must be
    // regenerated.
    //
    // WHY: the predict completion filter in train_predict.nf skips any genome
    // that already has a GBK unless staleRnaseq() or staleGenome() fires. Those
    // key off the READS and the SOURCE ASSEMBLY, not the training output, so
    // re-training a genome (new PASA models from the same reads, a repaired
    // pipeline, a fixed tool) leaves the filter satisfied and predict is never
    // instantiated. That is a SEPARATE blocker from the Nextflow cache issue
    // that trainingFingerprint() addresses, and it fires FIRST: without this
    // check the predict task is never created, so the cache never even gets
    // consulted. Confirmed 2026-09-20 on Chaetomium globosum in
    // v1.9.0-beta12_container_rust -- training was modified and the rerun
    // reported only "Cached process > FETCH_RNASEQ:COLLECT_SRA_QUERY", with no
    // FUNANNOTATE_PREDICT task in the DAG at all.
    static boolean staleTraining(String id, String trainingTarget, String target) {
        def gbk = gbkResult("${target}/${id}/predict_results", id)
        if (gbk == null) return false
        def gbkMod = gbk.lastModified()
        return ['funannotate_train.pasa.gff3', 'funannotate_train.transcripts.gff3'].any { name ->
            def f = new File("${trainingTarget}/${id}/training/${name}")
            f.exists() && f.size() > 0 && f.lastModified() > gbkMod
        }
    }

    // Fingerprint of the training evidence FUNANNOTATE_PREDICT will consume.
    //
    // WHY THIS EXISTS: FUNANNOTATE_PREDICT's inputs are all `val`
    // (meta, genome_fa, genemark_gtf, other_gff) and it reaches the training
    // directory at runtime through params.training_target, symlinking it into
    // the task dir. The training evidence is therefore INVISIBLE to Nextflow's
    // task hash -- so improving a genome's training set and re-running produces
    // a CACHE HIT on the old prediction and the new evidence is silently
    // ignored. Observed three times on 2026-09-19/20 in
    // BFD/Funannotate_benchmarking: Malassezia globosa was re-trained from 335
    // to 1,893 PASA models (after re-sourcing its RNA-seq) and predict was
    // skipped as cached against the old 227k-read result; Chaetomium globosum
    // hit the same thing in both beta.12 arms after an EVM-isolation test had
    // populated the cache. Note `-w <newdir>` does NOT avoid this: the resume
    // database lives in the launch dir, not the work dir.
    //
    // Passing this string into the predict input tuple puts the evidence into
    // the hash, so a changed training set correctly invalidates the cache.
    // Cheap by design (length + mtime, not a content digest) because it is
    // evaluated at channel-build time for every genome; funannotate rewrites
    // these files wholesale on each train, so length+mtime moves whenever the
    // evidence does. Returns "no-training" when the genome has none (no-RNAseq
    // assemblies), which is itself a stable, correct key.
    static String trainingFingerprint(String out, String trainingTarget) {
        def parts = []
        ['funannotate_train.pasa.gff3', 'funannotate_train.transcripts.gff3'].each { name ->
            def f = new File("${trainingTarget}/${out}/training/${name}")
            if (f.exists() && f.size() > 0) {
                parts << "${name}:${f.size()}:${f.lastModified()}"
            }
        }
        return parts ? parts.join('|') : 'no-training'
    }

    static String trainingTranscriptBamFor(String out, String trainingTarget) {
        def bam = new File("${trainingTarget}/${out}/training/transcript.alignments.bam")
        return (bam.exists() && bam.size() > 0) ? bam.toString() : ''
    }

    // A reuse_eligible strain's GBK must be considered stale if the shared
    // parameters.json it was predicted with has since been refreshed (representative
    // re-annotated, or the ANI/reuse assignment changed) -- mirrors staleRnaseq's
    // shape/purpose. sharedJsonPath as passed (nullable string) gating on emptiness.
    static boolean staleSharedParams(String out, String sharedJsonPath, String target) {
        if (!sharedJsonPath) return false
        def sharedJson = new File(sharedJsonPath)
        if (!sharedJson.exists() || sharedJson.size() == 0) return false
        def gbk = gbkResult("${target}/${out}/predict_results", out)
        if (gbk == null) return false  // predict hasn't run yet; normal path handles it
        return sharedJson.lastModified() > gbk.lastModified()
    }

    // Eagerly loads abinitio_reuse_csv into out -> [species, reuse_eligible,
    // is_representative, ani_to_representative]. Returns an empty map when the file
    // doesn't exist yet (PICK_REPRESENTATIVE_STRAIN hasn't run for this dataset) --
    // every row then falls through to independent-training behavior unchanged.
    // The caller gates on the share_abinitio_params param before invoking this.
    static Map loadAbinitioReuseMap(String csvPath) {
        def m = [:]
        def csv = new File(csvPath)
        if (!csv.exists() || !csv.canRead()) return m
        def lines = csv.readLines()
        if (lines.size() < 2) return m
        def header = lines[0].split(',', -1)*.trim()
        def iSpecies  = header.indexOf('species')
        def iOut      = header.indexOf('out')
        def iEligible = header.indexOf('reuse_eligible')
        def iIsRep    = header.indexOf('is_representative')
        def iAni      = header.indexOf('ani_to_representative')
        lines.drop(1).each { line ->
            def f = line.split(',', -1)
            if ([iSpecies, iOut, iEligible].any { it < 0 }) return
            if (f.size() <= [iSpecies, iOut, iEligible].max()) return
            // ani_to_representative is blank for the representative's own row (implicitly
            // 100%) and for any strain ANI couldn't confidently compare -- parsed as null
            // so callers can distinguish "no signal" from 0.0.
            def aniStr = iAni >= 0 && f.size() > iAni ? f[iAni].trim() : ''
            def isRep = iIsRep >= 0 ? f[iIsRep].trim().toLowerCase() == 'true' : false
            m[f[iOut].trim()] = [
                species              : f[iSpecies].trim(),
                reuse_eligible       : f[iEligible].trim().toLowerCase() == 'true',
                is_representative    : isRep,
                ani_to_representative: isRep ? 100.0 : (aniStr ? aniStr as Double : null),
            ]
        }
        return m
    }

    // Image for the PASA MariaDB steps (SETUP_MARIADB_DATADIR and the sidecar
    // branch of FUNANNOTATE_TRAIN / FUNANNOTATE_UPDATE). An explicit
    // params.container_mariadb wins; otherwise use params.container_funannotate,
    // which bundles MariaDB (mariadbd, mariadb-install-db, /usr/bin/mysqld_safe).
    static String mariadbImage(Map params) {
        return (params.container_mariadb ?: params.container_funannotate) as String
    }

    // Local file for an image. A plain path is returned unchanged. A URI
    // (docker://, oras://, library://, ...) maps to the file Nextflow's
    // apptainer cache uses for it under cacheDir (apptainer.cacheDir =
    // params.sif_dir), with the same naming as Nextflow's
    // SingularityCache.simpleName: strip the scheme, ':' and '/' -> '-', add
    // '.img'. e.g. docker://ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1 ->
    // <cacheDir>/ghcr.io-nextgenusfs-funannotate-1.9.0-rc.1.img
    // Scripts that call `apptainer exec` directly (outside Nextflow's
    // container handling) use this so they reuse the image Nextflow already
    // pulled instead of converting the URI on every call.
    static String localImageFile(String image, String cacheDir) {
        int p = image.indexOf('://')
        if (p == -1) return image
        String name = image.substring(p + 3)
        String ext = '.img'
        if (name.contains('.sif:')) {
            ext = '.sif'
            name = name.replace('.sif:', '-')
        } else if (name.endsWith('.sif')) {
            ext = '.sif'
            name = name.substring(0, name.length() - 4)
        }
        return "${cacheDir}/${name.replace(':', '-').replace('/', '-')}${ext}" as String
    }

    // Bash that makes sure localImageFile(image, cacheDir) exists before a
    // direct `apptainer exec`. A plain path must already exist. For a URI, the
    // image is pulled ONCE into the Nextflow cache file under a lock (so
    // parallel tasks do not pull it twice); later tasks and Nextflow itself
    // reuse that file.
    static String ensureLocalImageScript(String image, String cacheDir) {
        String f = localImageFile(image, cacheDir)
        if (f == image) {
            return """[ -e '${f}' ] || { echo "ERROR: container image not found: ${f}" >&2; exit 1; }"""
        }
        return """if [ ! -s '${f}' ]; then
    mkdir -p '${cacheDir}'
    (
        flock -w 7200 9 || { echo "ERROR: timed out waiting for the pull lock on ${f}" >&2; exit 1; }
        if [ ! -s '${f}' ]; then
            echo "[INFO] Pulling ${image} once into ${f}"
            _APPT=\$(command -v apptainer || command -v singularity)
            "\$_APPT" pull '${f}.pulling.'\$\$ '${image}' && mv -f '${f}.pulling.'\$\$ '${f}' \\
                || { rm -f '${f}.pulling.'\$\$; echo "ERROR: failed to pull ${image}" >&2; exit 1; }
        fi
    ) 9>'${f}.lock' || exit 1
fi"""
    }
}
