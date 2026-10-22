# Running one batch on more than one cluster

A large batch can be split between two clusters, for example a SLURM cluster and
NRP Nautilus (Kubernetes). Neither needs to be UCR HPCC: any SLURM site set up
as in [adding_a_site.md](adding_a_site.md) and any NRP namespace set up as in
[k8s/README.md](../k8s/README.md) will do.

The pipeline has no scheduler that spans clusters. Each cluster runs its **own,
independent** Nextflow run over a **different set of species**. You make the
split, in the sample sheets.

## Why split

| Work | Fits best | Why |
|---|---|---|
| Per-genome fan-out: masking, BUSCO, ANI, GeneMark, predict | the cluster with the most free cores; NRP is often best | CPU-bound and independent per genome. On NRP tasks start in seconds, but concurrency is set by `executor.queueSize` in your site config (profile default 100; be considerate, see [fair use](../k8s/README.md#nrp-fair-use-and-limits)). SLURM can run more at once but may queue |
| Genomes **without** RNA-seq (`--run_sra_fetch false`) | either; usually NRP | no SRA fetch, Trinity or PASA, so they finish much faster. Note this also turns **training** off (below) |
| Genomes **with** RNA-seq (SRA fetch or your own reads, Trinity, PASA) | the cluster that allows long, large-memory tasks | PASA and Trinity are memory-heavy. NRP bare pods are capped at 16 cores / 32 GB and 6 h unless you use Jobs |
| FCS-GX purge (470 GB database) | SLURM with a highmem node | not available on k8s |
| InterProScan, SignalP, DeepTMHMM, antiSMASH | SLURM | `--run_interpro true` stops the run with an error on `k8s`/`nrp`; the others are not available there. prodigal is unverified on k8s |

These are rules of thumb, not requirements. Queue wait usually decides it:
check how many tasks sit Pending on each side before assigning a big slice.

## How to split

1. **Slice the sample sheet by species.** Make one `samples.csv` per cluster with
   the same header. Put **every strain of a species in the same slice**, and no
   species in both:
   - RNA-seq reads, Trinity and PASA training are shared per species.
   - ANI reuse picks one representative strain per species, from the samples of
     *that run*. If a species is split, each side trains its own representative
     and they can overwrite each other's parameters when you merge.

   ```bash
   # genomes of the species listed in nrp_species.txt (one name per line) -> NRP; the rest -> SLURM
   awk -F, 'NR==FNR{s[$0];next} FNR==1||($1 in s)' nrp_species.txt samples.csv > samples_nrp.csv
   awk -F, 'NR==FNR{s[$0];next} FNR==1||!($1 in s)' nrp_species.txt samples.csv > samples_slurm.csv
   ```
   (Plain `awk -F,` is enough only if no field is quoted or contains a comma.)

   Running the same genome on both sides does not corrupt anything, because the
   clusters do not share a filesystem. It wastes compute and leaves two
   conflicting results to merge.
2. **One launch directory per cluster**, each with its own `samples.csv`,
   `rnaseq_reads/`, `work/`, params and `.nextflow` cache.
3. **Set the RNA-seq mode for each slice.** `--run_sra_fetch` applies to the
   **whole run**, not per genome, and your own reads are only used when it is `true`
   (see [RNA-seq strategies](../README_Extended.md#4-choose-where-the-rna-seq-reads-come-from)):
   - the RNA-seq slice: keep `run_sra_fetch` at its default `true`. For a species
     with your own reads, put the `_norm_` files and the header-only `sra_query` CSV
     in `rnaseq_reads/`. Leave a species empty to fetch from SRA;
   - the no-RNA-seq slice: `--run_sra_fetch false`. No training runs and any reads
     in `rnaseq_reads/` are ignored.
4. **Make cleaning consistent.** SLURM may run the FCS-GX purge while NRP cannot,
   so contamination screening would differ. Either use `--skip_fcs true` on both
   sides, or clean everything on SLURM first with `--only_clean true` and copy
   `input_clean_genomes/` into the NRP run directory. Cleaning is cached there and
   NRP then skips it.
5. **Use the same versions.** Check out the same pipeline tag or commit on both
   sides: `-r <tag>` for a run by project name on SLURM, and
   `git -C /data/src/nf_funannotate1 checkout <tag>` on the NRP PVC (the run Job
   uses that checkout, so `-r` has no effect). Set `container_funannotate`,
   `container_genome_clean` and the other `container_*` to the same tags on both.
   The k8s default for `container_genome_clean` is `ghcr.io/stajichlab/aaftf:latest`,
   whereas SLURM uses `<sif_dir>/AAFTF.sif`.
6. **Launch each side.**

   ```bash
   # SLURM cluster
   nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,singularity \
       -c mysite.config --samples samples_slurm.csv -resume
   ```

   NRP side ([k8s/README.md](../k8s/README.md) has the full setup):
   - Put the slice on the PVC as `<RUN_DIR>/samples.csv` (the default), or set
     `samples:` to an absolute path in the run's params YAML. YAML does not expand
     `${launchDir}`.
   - Copy the genomes onto the PVC. `GENOME` paths are relative to the launch
     directory. Use the `s3-sync` Job or `kubectl exec -i ... -- sh -c 'cat > f' < f`.
   - Copy your own reads into `<RUN_DIR>/rnaseq_reads/`.
   - Start from [`k8s/params_nrp_test.yaml`](../k8s/params_nrp_test.yaml). The `nrp`
     profile does not turn these off for you, so set `skip_fcs: true` (unless
     `input_clean_genomes/` is copied in) and `run_interpro`, `run_signalp`,
     `run_deeptmhmm`, `run_antismash` to `false`.
   - Copy `k8s/overlays/example/runs/my-run/`, set its params file and `RUN_DIR`,
     and `kubectl apply -k` it.

## What each side needs locally

Reference data is not shared between clusters. Each launch needs its own
`funannotate_db`, BUSCO lineages, `swissprot_fungi.faa`, `template.sbt` and
`taxondb`, plus its slice's genomes, `rnaseq_reads/` and (if used)
`input_clean_genomes/`. On NRP use the
[Reference DB archive](../k8s/README.md#reference-db-archive). The ANI-reuse store
is wherever `--gene_prediction_shared_abinitio` points; a path on one cluster
is not visible on the other. A species keeps a representative only within its
own run.

## Bring the results together

- Per-genome results (`genome_annotation/<Species_Strain>/` and the training
  directories) do not overlap between slices. Top-level files of the same name
  do: `genome_annotation/predict_skipped_too_small.tsv`,
  `predict_blocked_awaiting_representative.tsv`, `pasa_train_failed.tsv` and
  `rnaseq_*_candidates.csv`. Concatenate or rename them, do not overwrite.
- NRP output is on the PVC, which a SLURM node cannot see and `rsync` cannot
  reach. Copy it out with the `s3-sync` Job (PVC to a bucket), then pull it with
  `aws s3 sync --endpoint-url=https://s3-west.nrp-nautilus.io` where you need it.
- **Steps that only run on SLURM.** To add InterProScan, SignalP, antiSMASH or
  `annotate` to the genomes NRP predicted, copy those genomes' outputs to the
  SLURM launch directory and run again with the NRP slice's sample sheet.
  Genomes with a current `.gbk` skip prediction and go to annotation. Keep the
  file times (`rsync -a`, `cp -p`) and copy `input_clean_genomes/` too, because
  the skip checks compare file times and a newer input makes prediction run again.

## Moving a genome between clusters

If a genome stalls on one side, remove its species from that side's sheet and add it
to the other. Copy with times preserved (`rsync -a` or `cp -p`):
`genome_annotation/<Species_Strain>/`, `genome_annotation_training/<Species_Strain>/`,
`input_clean_genomes/<ASMID>.*`, and the species' `rnaseq_reads/*_norm_*` and
`rnaseq_data/<Species_name>.trinity-GG.fasta`. `-resume` will not help on the new
side, because the task cache belongs to the old launch directory and cluster.
Only the output checks skip finished work, and they re-run a step when an input
is newer than its output.
