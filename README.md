# nf_funannotate1

A Nextflow DSL2 framework for **eukaryotic genome annotation** (fungal defaults,
generic-capable) that runs on a workstation, any SLURM cluster, Kubernetes, or the
UCR HPCC. It runs the full [funannotate](https://github.com/nextgenusfs/funannotate)
workflow — genome clean → repeat mask → RNA-seq fetch/train → predict → optional
antiSMASH/InterProScan/SignalP → annotate/update — plus a standalone
[EarlGrey](https://github.com/TobyBaril/EarlGrey) curated repeat-masking pipeline.
The pipeline runs directly from GitHub; no clone is required.

## Install

- **Java 17+**
- **Nextflow 25.10.0+** (`curl -s https://get.nextflow.io | bash`)
- **Apptainer or Singularity** on `PATH` (on every compute node, if you use a cluster)

## Quick start

Run from a launch directory that holds `samples.csv`. Outputs, caches and `lib/`
assets live in the launch directory, not in Nextflow's cached copy of the repo
(`~/.nextflow/assets/`). Add `-r <branch|tag>` to pin a revision.

```bash
# Check the workflow graph first (synthetic data, no real tools):
nextflow run stajichlab/nf_funannotate1 -profile test -stub-run

# SLURM cluster (mysite.config holds your partitions and paths):
nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,singularity -c mysite.config -resume

# Workstation or single server:
nextflow run stajichlab/nf_funannotate1 -profile annotate,local,singularity -resume
```

On a cluster, run the Nextflow head process in a small, long job (about 2 CPUs and
32 GB of RAM, because image pulls run in the head process):

```bash
sbatch -c 2 --mem 32G --time 7-00:00:00 --wrap \
  "nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,singularity -c mysite.config -resume"
```

Minimal `samples.csv` (one row per genome; `GENOME` is a FASTA path, plain or gzipped):

```csv
SPECIES,STRAIN,ASMID,LOCUSTAG,BUSCO_LINEAGE,TRANSL_TABLE,NCBI_TAXONID,GENOME
Aspergillus nidulans,FGSC-A4,ANID_A4,ANIDA4,eurotiomycetes_odb10,1,227321,genomes/Aspergillus_nidulans_FGSC-A4.fasta
```

Minimal `mylab_params.yaml` (pass with `-params-file`; use absolute paths):

```yaml
sif_dir: /shared/containers          # container images and Apptainer pull cache
# proteins: /data/db/uniprot_sprot_fungi.faa   # default: <launch dir>/lib/swissprot_fungi.faa
```

RNA-seq comes from **one of two sources, chosen per species** (normally not mixed):

- **Your own reads:** place `<Species_name>_norm_R1/R2/SE.fastq.gz` and a header-only
  `sra_query/<Species_name>.sra_query.csv` plus its `.key` in `rnaseq_reads/`. Nothing is downloaded.
- **Public reads from NCBI SRA:** leave `rnaseq_reads/` empty for that species. The
  pipeline queries SRA, downloads and normalizes the reads into `rnaseq_reads/`.

Exact file names and an example are in
[README_Extended.md](README_Extended.md#4-choose-where-the-rna-seq-reads-come-from).
To test one genome first add
`--n_test 1`; to pick one use `--asmid ANID_A4`. After a stop or failure, run the
same command again and `-resume` continues. The full walkthrough (read layout,
`samples.csv` columns, params, container pulls) is in
[README_Extended.md](README_Extended.md#quick-start-local-genomes--local-rna-seq--singularity).

## Profiles

Compose one option from each axis: `-profile <pipeline>,<executor>,<provisioning>`.

| Axis | Options |
|---|---|
| **pipeline** | `annotate` · `earlgrey` · `test` / `stub` |
| **executor** | `slurm` · `local` |
| **provisioning** | `singularity` (containers; portable) · `pixi` · `conda` (shared frozen envs) · `ucr_hpcc` (UCR Lmod modules; institutional example) |

- **Your own SLURM cluster:** put partitions, account and the container cache in a
  small `mysite.config` (start from `conf/site_template.config`) and pass it with
  `-c`. See [docs/adding_a_site.md](docs/adding_a_site.md) and
  [Running on your own SLURM cluster](README_Extended.md#running-on-your-own-slurm-cluster).
- **Kubernetes (NRP Nautilus):** `-profile annotate,nrp` replaces the executor and
  provisioning picks. Any NRP user can use it, not only UCR. See
  [k8s/README.md](k8s/README.md).
- **Two clusters at once** (for example SLURM plus NRP): give each its own
  `samples.csv` slice and launch directory. See
  [docs/multi_cluster.md](docs/multi_cluster.md).

### UCR HPCC

UCR HPCC is the institutional example, not a requirement. Its Lmod profile,
launcher auto-detection and shared paths are in
[docs/site_ucr_hpcc.md](docs/site_ucr_hpcc.md).

## Outputs

In the launch directory: `input_clean_genomes/` (cleaned and masked genomes),
`rnaseq_data/` (Trinity transcripts), `genome_annotation_training/`, and
`genome_annotation/<Species_Strain>/predict_results/` (GFF3, GBK, protein and
transcript FASTA), plus `funannotate_db/` and `work/` caches (keep them for `-resume`).

## More documentation

- [README_Extended.md](README_Extended.md) — full walkthrough and reference:
  - [Quick start: local genomes + RNA-seq](README_Extended.md#quick-start-local-genomes--local-rna-seq--singularity)
  - [Orthogonal profiles](README_Extended.md#orthogonal-profiles) and
    [funannotate version / EVM selection](README_Extended.md#selecting-a-funannotate-version--evm-backend-conda-vs-singularity)
  - [Singularity images](README_Extended.md#singularity-images-to-build) and
    [PASA MySQL backend](README_Extended.md#pasa-mysql-backend-mariadb)
  - [Throughput, resumability & storage](README_Extended.md#throughput-resumability--storage)
  - [EarlGrey repeat masking](README_Extended.md#earlgrey-repeat-masking)
  - [Layout](README_Extended.md#layout), [Testing](README_Extended.md#testing)
- [docs/usage.md](docs/usage.md) — usage and parameters
- [docs/output.md](docs/output.md) — full list of outputs
- [docs/adding_a_site.md](docs/adding_a_site.md) — SLURM site configuration
- [docs/multi_cluster.md](docs/multi_cluster.md) — splitting one batch across clusters (e.g. SLURM + NRP)
- [docs/site_ucr_hpcc.md](docs/site_ucr_hpcc.md) — UCR HPCC as a worked example site
- [CHANGELOG.md](CHANGELOG.md) · [CITATIONS.md](CITATIONS.md)
