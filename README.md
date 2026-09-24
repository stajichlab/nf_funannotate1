# nf_funannotate1

A Nextflow DSL2 framework for **eukaryotic genome annotation** (fungal defaults,
generic-capable) on the UCR HPCC. It runs the full
[funannotate](https://github.com/nextgenusfs/funannotate) workflow — genome
clean → repeat mask → RNA-seq fetch/train → predict → optional
antiSMASH/InterProScan/SignalP → annotate/update — plus a standalone
[EarlGrey](https://github.com/TobyBaril/EarlGrey) curated repeat-masking pipeline
([see below](#earlgrey-repeat-masking)).

The pipeline lives at the repo root (`funannotate.nf` + `nextflow.config`), so it
runs directly from GitHub — no clone required:

```bash
nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,ucr_hpcc -resume
```

Nextflow caches the repo under `~/.nextflow/assets/`; add `-r <branch|tag>` to pin
a revision and `-latest` to pull updates. Outputs and the `samples.csv` /
`lib/` assets are read from your **launch directory**, not the cached checkout.

## Quick start: local genomes + local RNA-seq + Singularity

This example is for a new lab on its own Linux machine or cluster. It uses:

- a folder of genome FASTA files,
- a folder of RNA-seq FASTQ files,
- known species, NCBI taxon IDs and BUSCO lineages,
- Singularity/Apptainer containers for all tools.

It does not use NCBI downloads, environment modules or any site-specific paths.

### 1. Install the requirements

You need:

- Linux with **Apptainer** or **Singularity** on `PATH` (on every compute node, if you use a cluster).
- **Java 17 or newer**.
- **Nextflow 25.10.0 or newer**.
- Internet access on the first run. The run pulls container images and downloads the funannotate, taxonomy and BUSCO databases.

```bash
curl -s https://get.nextflow.io | bash      # installs ./nextflow
mv nextflow ~/bin/                          # or any directory on your PATH
nextflow -version

# Pull the pipeline. Nextflow caches it under ~/.nextflow/assets/.
nextflow pull stajichlab/nf_funannotate1

# Optional: test the workflow graph with synthetic data and no real tools.
nextflow run stajichlab/nf_funannotate1 -profile test -stub-run
```

You do not need a clone of the repository. For a pinned version, add `-r <tag>`.

### 2. Make a launch directory

Run the pipeline from a project directory. All outputs and caches go in this
directory. Example layout before the first run:

```
my_annotation/                          # launch directory: run nextflow from here
├── samples.csv                         # sample sheet (step 3)
├── mylab.config                        # site config (step 5)
├── lib/
│   └── swissprot_fungi.faa             # protein evidence FASTA (required, step 5)
├── genomes/                            # your genome FASTA folder
│   ├── Aspergillus_nidulans_FGSC-A4.fasta
│   ├── Aspergillus_nidulans_WT2.fa.gz
│   └── Fusarium_graminearum_PH-1.fasta
└── rnaseq_raw/                         # your RNA-seq FASTQ folder
    ├── Anid_mycelium_R1.fastq.gz
    ├── Anid_mycelium_R2.fastq.gz
    ├── Anid_conidia_R1.fastq.gz
    ├── Anid_conidia_R2.fastq.gz
    ├── Fgra_R1.fastq.gz
    └── Fgra_R2.fastq.gz
```

### 3. Write `samples.csv`

Use one row per genome. The `GENOME` column holds the path to the FASTA file.
A relative path resolves against the launch directory. The FASTA can be plain
or gzipped (`.fa`, `.fna`, `.fasta`, `.gz`).

```csv
SPECIES,STRAIN,ASMID,LOCUSTAG,BUSCO_LINEAGE,TRANSL_TABLE,NCBI_TAXONID,GENOME
Aspergillus nidulans,FGSC-A4,ANID_A4,ANIDA4,eurotiomycetes_odb10,1,227321,genomes/Aspergillus_nidulans_FGSC-A4.fasta
Aspergillus nidulans,WT2,ANID_WT2,ANIDW2,eurotiomycetes_odb10,1,162425,genomes/Aspergillus_nidulans_WT2.fa.gz
Fusarium graminearum,PH-1,FGRA_PH1,FGRAPH,hypocreales_odb10,1,229533,genomes/Fusarium_graminearum_PH-1.fasta
```

| Column | Meaning |
|---|---|
| `SPECIES` | Species name. It also sets the RNA-seq file name (step 4). |
| `STRAIN` | Strain or isolate name. The output directory is `<SPECIES>_<STRAIN>` with spaces changed to `_`. Give each genome of one species a different strain. |
| `ASMID` | Unique ID for this genome. It names the cleaned genome files. |
| `LOCUSTAG` | Locus-tag prefix for the gene IDs. |
| `BUSCO_LINEAGE` | BUSCO dataset, for example `fungi_odb10` or `eurotiomycetes_odb10`. |
| `TRANSL_TABLE` | NCBI genetic code. Use `1` for most fungi, `12` for CUG-Ser yeasts. |
| `NCBI_TAXONID` | NCBI taxonomy ID. Use the closest valid taxon for an unpublished genome. |
| `GENOME` | Path to the genome FASTA. |

### 4. Put the RNA-seq reads where the pipeline finds them

The pipeline looks for RNA-seq reads **per species** in `rnaseq_reads/`, by
file name. The name is `SPECIES` with spaces changed to `_`. For
`Aspergillus nidulans`, the name is `Aspergillus_nidulans`. All strains of one
species share the same reads.

Each species needs three read files and one query file:

| File | Content |
|---|---|
| `rnaseq_reads/<Species_name>_norm_R1.fastq.gz` | paired-end read 1 (all libraries joined) |
| `rnaseq_reads/<Species_name>_norm_R2.fastq.gz` | paired-end read 2, in the same order as R1 |
| `rnaseq_reads/<Species_name>_norm_SE.fastq.gz` | single-end reads, or an empty file |
| `rnaseq_reads/sra_query/<Species_name>.sra_query.csv` | a header-only CSV. It stops the NCBI SRA search and download for this species. |

```bash
cd my_annotation
mkdir -p rnaseq_reads/sra_query

# gzip files can be joined with cat. Keep R1 and R2 in the same library order.
cat rnaseq_raw/Anid_mycelium_R1.fastq.gz rnaseq_raw/Anid_conidia_R1.fastq.gz \
    > rnaseq_reads/Aspergillus_nidulans_norm_R1.fastq.gz
cat rnaseq_raw/Anid_mycelium_R2.fastq.gz rnaseq_raw/Anid_conidia_R2.fastq.gz \
    > rnaseq_reads/Aspergillus_nidulans_norm_R2.fastq.gz
: > rnaseq_reads/Aspergillus_nidulans_norm_SE.fastq.gz     # no single-end reads

cp rnaseq_raw/Fgra_R1.fastq.gz rnaseq_reads/Fusarium_graminearum_norm_R1.fastq.gz
cp rnaseq_raw/Fgra_R2.fastq.gz rnaseq_reads/Fusarium_graminearum_norm_R2.fastq.gz
: > rnaseq_reads/Fusarium_graminearum_norm_SE.fastq.gz

for sp in Aspergillus_nidulans Fusarium_graminearum; do
    echo 'species_tag,taxonid,sra_accession,spots,platform,layout' \
        > rnaseq_reads/sra_query/${sp}.sra_query.csv
done
```

For single-end data only, put the reads in `_norm_SE.fastq.gz` and make
`_norm_R1` and `_norm_R2` empty files.

Notes:

- Use real files, not symbolic links. A link to a path outside the launch
  directory may not be visible inside the containers.
- The pipeline does **not** trim or normalize these reads. It gives them to
  `funannotate train` as normalized reads. Trim adapters first. For very
  deep libraries, subsample or normalize them first.
- Training is skipped when a species has less than 50 MB of compressed reads
  in total (`--train_min_rnaseq_bytes`, default `50000000`).
- A species with no files in `rnaseq_reads/` makes the pipeline search NCBI
  SRA for public RNA-seq. To stop this for all species, use
  `--run_sra_fetch false`. That setting also turns off RNA-seq training for
  all species.

### 5. Write a site config (`mylab.config`)

The defaults contain some paths for the developers' cluster. Put your own
settings in one config file in the launch directory. Pass it with `-c`.

This example works on one workstation or a SLURM cluster. Change the values
that are marked.

```groovy
// mylab.config
params {
    // Directory for container images. Use shared storage on a cluster.
    sif_dir        = "/data/containers/singularity_cache"      // CHANGE

    // Protein evidence for funannotate predict (required).
    // Any protein FASTA works, e.g. UniProt/Swiss-Prot for your group.
    proteins       = "${launchDir}/lib/swissprot_fungi.faa"

    // Skip the NCBI FCS-GX contamination screen. It needs a ~470 GB database
    // and 500 GB of RAM. With this setting, cleaning only removes short contigs.
    skip_fcs       = true
    min_contig_len = 2000

    // Use the SQLite PASA backend. The MySQL backend (pasa_mysql = true)
    // needs no extra image: MariaDB is inside the funannotate image. See
    // "PASA MySQL backend (MariaDB)" below.
    pasa_mysql     = false

    // The default genome-clean image (AAFTF.sif) is not public. With
    // skip_fcs = true, this step needs only python3 and pigz. The public
    // funannotate image has both. Write the image name in full here:
    // other params are not yet defined when this file is read.
    container_genome_clean = 'docker://ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1'
}

process {
    // ASM_STATS runs in a minimal image by default. That image does not
    // have the tools it needs. Use the funannotate image instead.
    withName: 'ASM_STATS|.*:ASM_STATS' { container = 'docker://ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1' }

    // Maximum resources for one task. Set these to your largest node or to
    // your workstation. Tasks that ask for more are reduced to these values.
    resourceLimits = [ cpus: 32, memory: '128.GB', time: '72.h' ]      // CHANGE

    // ---- SLURM only (remove on a workstation) ----
    // queue          = 'batch'                     // your default partition
    // clusterOptions = '-N 1 -n 1 --account=mylab' // keep "-N 1 -n 1"
    // withName: '.*:FUNANNOTATE_PREDICT' { queue = 'long' }
}

// Optional: bind shared directories into every container. Use this when
// databases or genomes live outside the launch directory and $HOME.
// Keep the default options at the start of the string.
// apptainer.runOptions = '--nv ${SCRATCH:+-B $SCRATCH:$SCRATCH} -B /data:/data'

// Optional, SLURM only: limit how many jobs Nextflow submits.
// executor { queueSize = 100; submitRateLimit = '20/1min' }
```

On a workstation, the `local` profile runs at most 4 tasks at the same time.
To change this, add `executor { queueSize = 2 }`.

The pipeline makes no other site-specific assumptions:

- GeneMark runs from the public `teambraker/braker3` image. It does not need a license key.
- The pipeline downloads the funannotate, taxonomy and BUSCO databases on the first run. It caches them in the launch directory. To use copies that you already have, set `funannotate_db`, `taxondb` or `busco_lineages`.
- The functional steps are off by default: `--run_annotate`, `--run_antismash`, `--run_interpro` and `--run_signalp`. `--run_annotate` also needs an eggNOG database (`--eggnog_db`). SignalP 6 and DeepTMHMM need licensed images that you must build. See [One-time site artifacts](#one-time-site-artifacts-checklist).

For more settings, copy `conf/site_template.config` from the pipeline
(`~/.nextflow/assets/stajichlab/nf_funannotate1/conf/`).

### 6. Run

The profile has three parts: `annotate`, then an executor (`local` or `slurm`),
then `singularity`.

```bash
cd my_annotation

# Workstation or single server:
nextflow run stajichlab/nf_funannotate1 \
    -profile annotate,local,singularity \
    -c mylab.config -resume

# SLURM cluster: run the Nextflow head process in a long, small job.
# It needs about 2 CPUs and 32 GB of RAM, because image pulls run in the head process.
sbatch -c 2 --mem 32G --time 7-00:00:00 --wrap \
  "nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,singularity -c mylab.config -resume"
```

To test one genome first, add `--n_test 1`. To select one genome, use
`--asmid ANID_A4`. After a stop or a failure, run the same command again.
`-resume` continues from the last completed step.

### 7. Outputs

After the run, the launch directory contains:

```
my_annotation/
├── input_clean_genomes/        # cleaned and repeat-masked genomes (<ASMID>.fa.gz, <ASMID>.masked.fasta.gz)
├── rnaseq_data/                # Trinity transcripts per species (<Species_name>.trinity-GG.fasta)
├── genome_annotation_training/ # funannotate train output per genome
├── genome_annotation/          # funannotate predict output per genome:
│   └── Aspergillus_nidulans_FGSC-A4/predict_results/  # GFF3, GBK, protein and transcript FASTA
├── funannotate_db/  work/      # database and Nextflow caches (keep them for -resume)
└── .nextflow.log
```

See [`docs/output.md`](docs/output.md) for the full list of outputs.

### Other launch modes

```bash
# UCR HPCC (Lmod modules on SLURM), from a launch dir with samples.csv:
nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,ucr_hpcc -resume --n_test 1

# or, from a local checkout, use the UCR sbatch launcher
sbatch /path/to/nf_funannotate1/run_annotate.sh --n_test 1
```

## Orthogonal profiles

Compose one option from each of three axes: `-profile <pipeline>,<executor>,<provisioning>`

| Axis | Options |
|---|---|
| **pipeline** | `annotate` · `earlgrey` · `test` / `stub` |
| **executor** | `slurm` · `local` |
| **provisioning** | `ucr_hpcc` (default; institutional Lmod modules) · `conda` (shared frozen envs) · `pixi` · `singularity` (containers) |

```bash
nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,ucr_hpcc -resume
nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,conda -resume
nextflow run stajichlab/nf_funannotate1 -profile annotate,local,singularity -resume
```

The `run_annotate.sh` launcher honours `EXECUTOR=` / `PROVISION=` (default
`slurm` / `ucr_hpcc`) and `PIPELINE=` / `REVISION=` env vars:

```bash
sbatch run_annotate.sh                          # default: slurm + ucr_hpcc (Lmod modules)
PROVISION=conda       sbatch run_annotate.sh     # shared conda envs
PROVISION=singularity sbatch run_annotate.sh     # portable containers
EXECUTOR=local         sbatch run_annotate.sh     # head + tasks local
```

`PROVISION=conda`/`singularity` **layers on top of** `ucr_hpcc` rather than
replacing it (`-profile annotate,slurm,ucr_hpcc,conda`) — `ucr_hpcc` still
carries the SLURM SCRATCH/TMPDIR safety net regardless of which provisioning
axis supplies the tools; only the last-loaded axis's `beforeScript`/`container`
setting wins per process.

It runs the pipeline **by project name** (`nextflow run stajichlab/nf_funannotate1`)
rather than by file path, so it is safe under `sbatch` (which copies the script
to a spool dir). For development, point it at a local checkout:
`PIPELINE=$PWD sbatch run_annotate.sh`.

Process scripts carry **no `module load`** — provisioning is supplied per process
`label` by the provisioning profile (`conf/provision_*.config`): a `beforeScript`
(ucr_hpcc/pixi/conda) or a `container` (singularity).

### Selecting a funannotate version / EVM backend (conda vs. singularity)

Neither axis auto-selects a version — pick one explicitly per run:

| | conda (`--conda_env`) | singularity (`--container_funannotate`) |
|---|---|---|
| 1.8.17 | `funannotate-1.8.17` | local `.sif` pulled from `docker://nextgenusfs/funannotate:v1.8.17` (Docker Hub only — ghcr has no 1.8.17 tag; build with `scripts/pull_funannotate_image.sh` or see `conf/release_1_8.config`) |
| 1.9.0-rc.1, perl EVM | `funannotate-1.9.0-rc.1` (`params.conda_env` default) | `docker://ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1-norust` |
| 1.9.0-rc.1, rust EVM | `funannotate-1.9.0-rc.1-rust` | `params.container_funannotate` default: `docker://ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1` (rust-enabled) |

```bash
# 1.8.17 via conda
sbatch run_annotate.sh --conda_env funannotate-1.8.17 -c conf/release_1_8.config
PROVISION=conda sbatch run_annotate.sh --conda_env funannotate-1.9.0-rc.1-rust

# 1.9.0-rc.1, perl EVM, via singularity (no-rust image)
PROVISION=singularity sbatch run_annotate.sh \
    --container_funannotate docker://ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1-norust
```

`conda_env` resolves under `--conda_envs_root` (`$CONDA_ENVS_ROOT`, default
`/bigdata/stajichlab/shared/condaenv` when set by `run_annotate.sh`); see
`conf/provision_conda.config`. `container_funannotate` resolves relative to
`--sif_dir` when not given as an absolute path/URI; see
`conf/provision_singularity.config`.

### Singularity images to build

Nextflow pulls these public images automatically into the `sif_dir` cache:

- `funannotate` — `ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1` (includes MariaDB, see below)
- `sra` — `ghcr.io/hyphaltip/sra_tools_container/sra_tools:1.3.1`
- edirect, prodigal, skani, busco, interproscan, setup, braker3 (biocontainers / Docker Hub)

These images are not public. Build them and point at them with
`--container_*` (defaults under `/bigdata/stajichlab/shared/lib/singularity_cache`):
`AAFTF` (genome_clean), `signalp6-fast.sif` (fast mode, licensed),
`DeepTMHMM-1.0.sif` (licensed), and the antismash-procps image. Build
one-liners live in `conf/provision_singularity.config` next to each
`container_*` param.

### PASA MySQL backend (MariaDB)

MariaDB is part of the funannotate image. `ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1`
(and `-norust`) contains `mariadbd` 11.8.6, `mariadb-install-db`, `mariadb`
and `/usr/bin/mysqld_safe`. With `--pasa_mysql true` under `-profile singularity`,
`FUNANNOTATE_TRAIN` and `FUNANNOTATE_UPDATE` run inside that image, find these
binaries on `PATH`, and start MariaDB in the image. They do not use a separate
MariaDB container.

The other MariaDB steps use `params.container_mariadb`: `SETUP_MARIADB_DATADIR`
(runs `mariadb-install-db`), and the sidecar branch of `FUNANNOTATE_TRAIN` /
`FUNANNOTATE_UPDATE` (used when the task itself has no MariaDB, e.g. under
`-profile ucr_hpcc` or `conda`). These steps call `apptainer exec` directly.

- **Default:** `container_mariadb` is unset, so these steps use
  `params.container_funannotate`. No separate MariaDB image is needed at any site.
- **Override:** `--container_mariadb <path.sif | docker://uri>`. The old UCR
  image (MariaDB 10.3.9) still works:
  `--container_mariadb /bigdata/stajichlab/shared/lib/singularity_cache/mariadb.sif`.
- **`docker://` URIs are not converted on every call.** The pipeline maps the URI
  to the file Nextflow's own image cache uses in `sif_dir` (for example
  `ghcr.io-nextgenusfs-funannotate-1.9.0-rc.1.img`). If that file is missing,
  the first task pulls it once under a file lock; later tasks and Nextflow
  reuse it. (A direct `apptainer exec docker://...` of the funannotate image
  took 885 s and about 15 GB of temporary space on UCR HPCC.)

## Running at another site

The pipeline is portable by construction: process scripts carry **no tools**
and the provisioning is injected per process `label` by the provisioning
profile (see above), so everything except the `ucr_hpcc` (Lmod) axis works
identically at any site. To run elsewhere:

```bash
# single workstation / login node (containers, any linux with apptainer):
nextflow run stajichlab/nf_funannotate1 -profile annotate,local,singularity

# another SLURM cluster (containers on the scheduler):
nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,singularity

# ...or project-local pixi envs instead of containers (any executor):
nextflow run stajichlab/nf_funannotate1 -profile annotate,local,pixi
```

`local` is the portable executor; the `slurm` profile was already generic
(`-N 1 -n 1` clusterOptions moved into it — no SLURM flags leak into other
executors). Non-SLURM schedulers (PBS/LSF/SGE) need only a small executor block
like the `slurm` one in `nextflow.config`.

### One-time site artifacts (checklist)

1. **Container cache** — `sif_dir` resolves from `NXF_APPTAINER_CACHE` /
   `NXF_SINGULARITY_CACHE` / `APPTAINER_CACHE` (fallback documented in
   `nextflow.config`). Export one of those env vars. Nextflow pulls the public
   images on first use: `funannotate:1.9.0-rc.1` (ghcr.io, includes MariaDB),
   `sra_tools:1.3.1` (ghcr.io), and the biocontainers (edirect, prodigal,
   skani, busco, interproscan, setup, braker3). You must build the non-public
   images listed in `conf/provision_singularity.config` yourself: `AAFTF`
   (or use the `skip_fcs` workaround in the quick start),
   `signalp6-fast.sif` + converted GPU weights and `DeepTMHMM-1.0.sif`
   (licensed), and `antismash-standalone-8.0.4-procps.sif`. No `mariadb.sif`
   is needed: MariaDB comes from the funannotate image (see
   [PASA MySQL backend](#pasa-mysql-backend-mariadb)).
2. **Reference databases** — everything the pipeline can auto-download
   (taxonkit taxdump, funannotate DBs, antiSMASH DBs, BUSCO lineages) is
   storeDir-cached under `launchDir` and needs no manual step. What a site
   *must* supply before those processes run: the **eggNOG DB** (`--eggnog_db`:
   eggnog.db, eggnog_proteins.dmnd, eggnog.taxa.db) and **SwissProt fungi
   proteins** (`--proteins`) for `--run_annotate`; a **PASA config**
   (`--pasa_conf_dir`) for the MySQL backend. Point an existing shared install
   at any of them with `--<param> /path` and the setup process no-ops.
3. **Licenses / GPU** — a GeneMark `~/.gm_key` (or `--genemark_path` to a
   host install, vs the braker3-container mode `genemark_container_mode=true`
   on the singularity axis); `signalp6-fast.sif` and `DeepTMHMM-1.0.sif` are
   licensed (no conda/module substitute). Turn `--signalp_gpu` /
   `--deeptmhmm_gpu` off on clusters without GPU nodes — same image, CPU mode.
4. **Conda axis only** — build the frozen envs once into shared storage with
   `environments/conda/build_conda_env.sh`. Point every run at it the same way
   the container cache is pointed: export `CONDA_ENVS_ROOT=/shared/lib/condaenv`
   in the site env/launcher, or pass `--conda_envs_root` per run. Unset, the
   conda axis falls back to a per-user `$HOME/.conda/nf_funannotate1`.
5. **Site config file** — copy `conf/site_template.config` to
   `conf/site_<site>.config`, fill in the shared paths you have, and pass
   `-c conf/site_<site>.config` (or includeConfig it from a site profile).
   `conf/site_ucr_hpcc.config` is the filled-in UCR example.

## Input model (`samples.csv`)

Columns: `SPECIES, STRAIN, ASMID, LOCUSTAG, BUSCO_LINEAGE, TRANSL_TABLE, NCBI_TAXONID, GENOME`

Genome resolution is **dual**:
- a non-empty **`GENOME`** column → use that local FASTA directly (`.fa`/`.fna`,
  gzipped or plain; relative paths resolve against the launch dir);
- otherwise resolve `<source>/<ASMID>/<ASMID>_genomic.fna.gz` from the NCBI_ASM
  `--source` dir.

Useful filters: `--taxon RANK:VALUE`, `--asmid <ASMID>`, `--n_test N`, a
`suppress.txt` ASMID skip-list.

## Throughput, resumability & storage

The pipeline is built to run over thousands of genomes and survive walltime
kills / orchestrator restarts. Four subsystems make this practical; all are on by
default and tunable from `conf/profile_annotate.config` (or `--<param>`).

### Batched genome cleaning (FCS-GX)

Cleaning stages the ~470 GB NCBI **FCS-GX** database into `/dev/shm` (~30 min).
Paying that per genome is wasteful at scale, so by default genomes are grouped
into one SLURM job that stages the DB once and cleans the whole batch
sequentially (`GENOME_CLEAN_BATCH`).

| Param | Default | Effect |
|---|---|---|
| `clean_batch_size` | `1000` | genomes per batch job; `0` → one job per genome (`GENOME_CLEAN`) |
| `skip_fcs` | `false` | bypass FCS-GX entirely (no gxdb / highmem); also forces the per-genome path |

- Already-cleaned genomes are skipped, so a killed batch resumes without redoing
  finished assemblies, and a fully-clean batch is never scheduled (no staging cost).
- Each batch writes a manifest (`clean_batch_*.manifest.tsv`) of what it cleaned.
- **Set `FCS_GX_DB_SRC`** to your gxdb path (see `scripts/setup_fcs_shm.sh`).

### Resumable, persistent prediction

`FUNANNOTATE_PREDICT` computes directly into the durable per-genome dir
(`<target>/<out>/`) and emits a small `<out>.predict.done` marker — there is no
publishDir copy or rsync. funannotate checkpoints into `predict_misc/`, so a job
killed by OOM/timeout resumes completed steps in place on the next run. A current
GBK short-circuits; an RNA-seq / Trinity input newer than the GBK forces a clean
re-predict. Genomes predicted during a run flow straight into the optional
annotate / antiSMASH / InterProScan / SignalP / update steps in the **same** run.

### Too-small / fragmented pre-flight guard

Assemblies that are both small *and* fragmented cannot yield funannotate's 30
required training models and would burn hours before aborting. They are detected
up front (and again from the predict log) and skipped cleanly — flagged in
`<target>/predict_skipped_too_small.tsv` — instead of failing the batch.

| Param | Default | Meaning |
|---|---|---|
| `predict_min_asm_bp` | `8000000` | below this assembled size = "small" (`0` disables the guard) |
| `predict_frag_max_n50` | `10000` | N50 below this = "fragmented" |
| `predict_frag_max_contigs` | `1000` | contig count above this = "fragmented" |

Both the small *and* fragmented gates must trip, so complete small genomes (e.g.
*Malassezia*) are unaffected.

### Standalone ab-initio predictors (GeneMark + optional Prodigal)

Beyond funannotate's own Augustus training, two standalone predictors run ahead of
`FUNANNOTATE_PREDICT` and hand their models over as EVM sources.

**GeneMark-ES/ET (`GENEMARK_RUN`, on by default).** Trains ab-initio models on the
masked genome (ES) or with RNA-seq hints (ET) and passes the `.mod` to predict
via `--genemark_gtf` (EVM source `genemark`, weight 1). Obeys `run_genemark`,
`genemark_mode` (AUTO/ES/ET), and `genemark_container_mode`. The container path
(validated against a real genome in the 2026-08-20 realtest) runs
`gmes_petap.pl` from the public `teambraker/braker3` image, which bundles
GeneMark 4.72 under BRAKER's CC-BY-NC-SA relicense — no `gm_key` needed — plus
Augustus's `bam2hints`/`join_mult_hints.pl` for ET hint prep. The host-licensed
`genemarkESET` module / a private `genemark_path` install is the alternative for
sites that have a key.

**Prodigal (`PRODIGAL_RUN`, off by default).** For short-ORF / intron-poor
special cases — microsporidia (tiny, very-compact genomes where neither Augustus
nor usually GeneMark recover enough training genes for the 30-model floor) —
Prodigal's single-genome mode is strong at short-ORF discovery. Its GFF is
passed to predict as `--other_gff <gff>:<weight>`, which funannotate renames to
EVM source `other_pred1` (`StartWeights['other_pred1']`). `-g` reuses the
assembly's `transl_table`, so alternative / CUG-reassigned tables always agree
with predict's `--table`. Follow-up candidate: Saccharomycotina.

| Param | Default | Meaning |
|---|---|---|
| `run_genemark` | `true` | run the standalone GeneMark step; `false` = `--auto-skip-genemark` (no ab-initio models) |
| `genemark_mode` | `AUTO` | ES / ET / AUTO (ET when a training BAM exists) |
| `genemark_container_mode` | `false` | run GeneMark inside `container_genemark` (braker3) instead of a host licensed install |
| `run_prodigal` | `false` | run `PRODIGAL_RUN` and hand its models to predict via `--other_gff` |
| `prodigal_weight` | `5` | EVM weight for `other_pred1` (calibrated 2026-08-31 on OC4: gene Sn/Sp 0.839/0.817) |
| `prodigal_mode` | `single` | `single` (complete genomes) or `meta` (fragmented/short contigs) |
| `prodigal_lineages` | `[]` | contingency: run `PRODIGAL_RUN` only for genomes whose BUSCO lineage is in this list (e.g. `['microsporidia_odb10']`); empty = all genomes |

Prodigal runs on the same soft-masked genome as the other steps, so it will
under-call inside repeats (acceptable; EVM-compatible) — running it on the
unmasked raw with a mask-interval filter is a documented in-module TODO.

### Compressed storage

Clean and masked genomes are stored gzip-compressed in `input_clean_genomes/`
(`<asmid>.fa.gz`, `<asmid>.masked.fasta.gz`); tools that can't read gzipped FASTA
inflate a local copy on the fly. Completion gating accepts either `.gbk` or
`.gbk.gz`, so finished annotation folders can be archived/compressed without
breaking skip logic on the next run. Legacy uncompressed `.fa` files are still
recognized.

## EarlGrey repeat masking

`earlgrey_mask.nf` builds a curated TE library once per species on the best
representative genome (`> --cutoff_mb`), then applies it to every conspecific
strain with RepeatMasker, writing `<asmid>.masked.fasta.gz` into
`input_clean_genomes/`. funannotate consumes that file in place of its default
tantan mask wherever it exists.

```bash
nextflow run earlgrey_mask.nf -c nextflow.config -profile earlgrey -resume
# restrict to the species that owns one assembly (representative or member):
nextflow run earlgrey_mask.nf -c nextflow.config -profile earlgrey --asmid GCA_XXXXXXXXX.1 -resume
```

EarlGrey runs into a persistent per-species dir (`params.earlgrey_workdir`) so a
walltime-killed run resumes from its checkpoints instead of restarting the
multi-hour discovery; its `-M` memory cap is derived from the SLURM allocation.
Tune `--cutoff_mb`, `--repeat_taxon`, and `--n_test`.

## External / built tools

### Rust helpers (built on deploy, not committed)

Two Rust binaries used by the SRA/RNA-seq steps are **built from source** into
`tools/` (gitignored) rather than checked in (they are dynamically-linked,
platform-specific ELFs). Build them inside the pipeline checkout — for a GitHub
run that is the cached asset dir (`~/.nextflow/assets/stajichlab/nf_funannotate1`):

```bash
module load rust            # Rust toolchain, edition 2024+ (rust >= 1.85)
bash scripts/build_tools.sh
```

| Tool | Source | Used as |
|---|---|---|
| `fix_fastq_header_trinity` | https://github.com/hyphaltip/fix_fastq_header_trinity | `params.fastq_hdr_script` |
| `enforce_seqpair_readlen` | https://github.com/hyphaltip/enforce_seqpair_readlen | `params.readlen_script` |

Revisions are pinned in `build_tools.sh` (override with `FIXHDR_REV` /
`ENFORCE_REV`). Each tool ships a Python fallback (`scripts/enforce_seqpair_readlen.py`,
upstream `fix_fastq_headers.py`) if you can't build the Rust version.

### Site data scripts

- `scripts/clean_genome_fa.py` — min-length contig filter (stdlib only).
- `scripts/setup_fcs_shm.sh` — stages the NCBI **FCS-GX** database into
  `/dev/shm` for `GENOME_CLEAN`. **Set `FCS_GX_DB_SRC`** to your gxdb path.

## Layout

```
nf_funannotate1/           # repo root = pipeline root (runs from GitHub)
  nextflow.config                 # manifest (mainScript=funannotate.nf), shared params, profiles map, singularity block
  funannotate.nf                  # full annotation workflow (labeled processes) — default entry
  earlgrey_mask.nf                # standalone curated repeat masking (EarlGrey)
  conf/
    profile_annotate.config       # params + per-process resources
    provision_{module,pixi,singularity}.config
    profile_earlgrey.config
    test.config                   # self-contained stub profile
  lib/SampleUtils.groovy          # auto-compiled by Nextflow (projectDir/lib)
  scripts/                        # clean_genome_fa.py, setup_fcs_shm.sh, build_tools.sh, *.py fallbacks
  pixi.toml                       # per-label conda envs for the pixi profile
  run_annotate.sh, run_earlgrey.sh
  tools/bin/                      # built Rust helpers (gitignored)
  tests/data/                     # synthetic stub fixtures
```

## Testing

`-profile test -stub-run` exercises the whole graph with synthetic data and no
real tools. The SRA/RNA-seq subgraph is exercised with `--run_sra_fetch true`.
Post-predict steps (antismash/interpro/signalp/annotate) run in the same pass for
genomes predicted in that run, and skip cleanly for genomes already complete
(`predict_results/*.gbk` or `.gbk.gz` present and not stale).
