// End-to-end test of the GENEMARK_RUN subworkflow (train once, scatter predict, gather) with the real profiles.
// Two rows in one batch, both fresh ES training (force_independent = true, as genemark_sidecar.nf does):
//   chunked_k4  1,759-contig genome; with --genemark_chunk_min_contigs 1000 --genemark_chunk_contigs 500 it takes the chunked route (4 chunks)
//   mono_small  30-contig genome; below the threshold, takes the single-task route
// Equivalence (done after the run, see README.md): run `gmes_petap.pl --predict_with <published chunked_k4 .mod>` on the
// whole genome and compare with the chunked GTF. ES training is not deterministic, so the same published model is the reference.
nextflow.enable.dsl = 2
include { GENEMARK_RUN } from './subworkflows/local/genemark_run'
workflow {
    def d = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_annotation_wave2/input_clean_genomes'
    rows = channel.of(
        [[id: 'chunked_k4',  species: 'Test chunked', transl_table: '1'], "${d}/GCA_000308955.1_ETyphina_1.0.masked.fasta.gz", 'ES', '', 'true', '', ''],
        [[id: 'mono_small',  species: 'Test mono',    transl_table: '1'], "${d}/GCA_001599415.1_JCM_22320_assembly_v001.masked.fasta.gz", 'ES', '', 'true', '', ''],
    )
    GENEMARK_RUN(rows)
    GENEMARK_RUN.out.gtf.view { meta, f -> "GTF ${meta.id}" }
    GENEMARK_RUN.out.mod.view { meta, f -> "MOD ${meta.id}" }
}
