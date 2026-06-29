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
}
