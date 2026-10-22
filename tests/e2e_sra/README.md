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
mkdir -p lib && cp /path/to/swissprot_fungi.faa lib/       # protein evidence (or set PROTEINS=)
REPO=/path/to/nf_funannotate1 sbatch $REPO/tests/e2e_sra/run_e2e.sbatch
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

Same sample sheet and params on NRP: copy this directory to the PVC run directory
(`/data/runs/sra_e2e`), put the genomes in `genomes/`, and use the run Job from
[k8s/README.md](../../k8s/README.md) with `run_sra_fetch: true` (the "Turning features back on"
step 1). Needs `sra_tools` 1.4.0 (the default), and no queue limit of more than your
`executor.queueSize`. Run `check.sh` in the run directory from a pod
(`kubectl exec ... -- bash check.sh`). Trinity needs the larger task sizes, so set
`k8s.computeResourceType = 'Job'` and raise `k8s_max_cpus` / `k8s_max_memory` as k8s/README.md
describes, or expect the PASA step to be capped at 16 cores / 32 GB.
