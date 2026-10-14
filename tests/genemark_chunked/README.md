# Chunked GeneMark test (HPCC)

    cd <repo>; module load nextflow
    nextflow run test_genemark_chunked.nf \
        -c nextflow.config -profile genemark_sidecar,slurm,singularity \
        -c conf/site_ucr_hpcc.config -c conf/site_ucr_hpcc_singularity.config \
        -c tests/genemark_chunked/test_overlay.config --target <outdir>

Expected: `chunked_k4` takes the chunked route (4 chunks; one row per chunk in `<outdir>/genemark_chunked.tsv`), `mono_small`
takes the single-task route. Equivalence: ES training is not deterministic, so take the published `chunked_k4.genemark.mod`, run
`gmes_petap.pl --predict_with <that .mod> --sequence <genome> --max_intron 3000 --soft_mask 2000 --fungus` on the whole
genome, and compare with the chunked GTF (identical records, order and IDs expected).

## Same-model monolithic reference (inside a SLURM job, 8 CPUs, node-local `$SCRATCH`)

    module load apptainer
    WORK=$(mktemp -d "$SCRATCH/gm_mono.XXXXXX"); cd "$WORK"
    gzip -dc <masked genome .fasta.gz> > genome.fa; cp <published chunked_k4.genemark.mod> genemark-shared.mod
    apptainer exec --bind "$WORK:$WORK,$(dirname "$(readlink -f "$HOME/.gm_key")"),$SCRATCH" <braker3 sif> \
        gmes_petap.pl --predict_with genemark-shared.mod --sequence genome.fa --max_intron 3000 --soft_mask 2000 \
        --cores 8 --fungus
    # compare genemark.gtf with <outdir>/published/chunked_k4.genemark.gtf (identical records, order and IDs expected)

Add `--genemark_local_scratch true` to the nextflow command to test the node-local-scratch path.
