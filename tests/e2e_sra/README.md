# End-to-end SRA test (Yarrowia lipolytica + Cryptococcus deneoformans)

Runs the real chain on two small genomes: SRA query, SRA fetch, Trinity/PASA training and
predict. *Cryptococcus deneoformans* also covers the large-result and mutant/host cases that the
SRA query was fixed for (see `tests/live_sra_query.sh`, which tests the query alone in minutes).
Plan a few hours on a cluster; the genomes are about 20 Mb each.

Written for UCR HPCC (SLURM, containers) and runnable on any site that has the pipeline's
[site config](../../docs/adding_a_site.md); the NRP variant is at the end.

## One job for the whole sequence

`run_e2e.sbatch` does steps 1 to 5 below in one job (a failed check stops it):

```bash
mkdir -p /path/to/sra_e2e/logs && cd /path/to/sra_e2e      # an empty launch directory
REPO=/path/to/nf_funannotate1 PROTEINS=/path/to/swissprot_fungi.faa SIF_DIR=/path/to/singularity_cache \
    sbatch $REPO/tests/e2e_sra/run_e2e.sbatch
# UCR HPCC: SIF_DIR=/bigdata/stajichlab/shared/lib/singularity_cache PROTEINS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/swissprot_fungi.faa
# (or put the file at lib/swissprot_fungi.faa in the launch directory, the default)
```

Add `SITE_CONFIG=/path/mysite.config` away from UCR HPCC (UCR HPCC is detected and its
profile added), `STOP_AFTER=2` to stop after the first run, or `sbatch -p <partition>`.
The log is `logs/sra_e2e.<jobid>.log`. The manual steps follow.

## 1. Set up a launch directory

```bash
mkdir -p /path/to/sra_e2e && cd /path/to/sra_e2e
cp /path/to/nf_funannotate1/tests/e2e_sra/{samples.csv,params.yaml,fetch_genomes.sh,check.sh} .
bash fetch_genomes.sh                    # two RefSeq genomes into genomes/
```

Your site settings go in `mysite.config` (partitions, account, `sif_dir`; see
`conf/site_template.config`) and, on UCR HPCC, `-profile annotate,slurm,ucr_hpcc`
or `singularity` (see [docs/site_ucr_hpcc.md](../../docs/site_ucr_hpcc.md)). The pipeline also needs
`lib/swissprot_fungi.faa` or `--proteins`.

## 2. First run

```bash
sbatch /path/to/nf_funannotate1/run_annotate.sh -params-file params.yaml        # UCR HPCC launcher
# or, any site:
nextflow run /path/to/nf_funannotate1 -profile annotate,slurm,singularity -c mysite.config \
    -params-file params.yaml -resume
bash check.sh                            # query CSV + key, reads, Trinity, one GBK per genome
```

## 3. Rerun with nothing cached in Nextflow

```bash
rm -rf work .nextflow                    # keep rnaseq_reads/, rnaseq_data/, genome_annotation/
nextflow run /path/to/nf_funannotate1 ... -params-file params.yaml          # same command
bash check.sh rerun                      # queries and reads came from the caches, nothing downloaded
```

## 4. Own-reads path (the recipe in README_Extended step 4)

Treat Yarrowia's SRA reads as your own: the query CSV is made header-only and its key is kept.

```bash
f=rnaseq_reads/sra_query/Yarrowia_lipolytica.sra_query
head -1 $f.csv > $f.tmp && mv $f.tmp $f.csv      # header-only; $f.key stays
rm -rf work .nextflow
nextflow run /path/to/nf_funannotate1 ... -params-file params.yaml
bash check.sh own-reads                  # "Reusing cached result", no NCBI query for Yarrowia
```

Then delete `$f.key` and run once more: the species must now be **searched again** (the key
rule), and the CSV is rewritten with runs.

## NRP variant

Run on NRP 2026-10-22 and 2026-10-23: all five steps pass (11, 4, 2, 3 checks). Layout and commands:

- Run directory `/data/runs/sra_e2e` on the PVC: copy `samples.csv`, `params.yaml`, `check.sh` and
  `fetch_genomes.sh` there (`fetch_genomes.sh` works from a shell pod, which has `curl`), set a BUSCO lineage that is on
  the PVC (the sample sheet's `tremellomycetes_odb10` is not; use `fungi_odb10`), and start the run Job from
  [k8s/README.md](../../k8s/README.md) (`kubectl apply -k <overlay>/runs/sra_e2e`).
- The PVC checkout of this repo must have the funannotate image fix (a `docker://` default for k8s/NRP); without it
  task pods sit in `InvalidImageName`.
- Steps 3 to 5 are the same as above, done by hand: `kubectl delete -k` the finished Job, change the run directory
  (steps 3 to 5 start with `rm -rf work .nextflow .nextflow-snapshots .nf-run.heartbeat`; step 4 makes Yarrowia's
  query CSV header-only and keeps its key; step 5 also deletes the key), `kubectl apply -k` again, then run
  `bash check.sh <mode>` in the run directory from a shell pod. `check.sh` needs bash >= 4 and GNU `stat` only.
- Run 1 took about 2 hours of wall time in the first attempt and an hour after the restart; later steps about
  15 to 30 minutes each, mostly pod start-up.
- A node that hangs pulling the 4 GB funannotate image can stall a run for over an hour: see "A node hangs pulling
  the image" in k8s/README.md and `k8s_avoid_nodes`.
