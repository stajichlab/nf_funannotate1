# nf_euk_genome_annotate

A Nextflow DSL2 framework for **eukaryotic genome annotation** (fungal defaults,
generic-capable) on the UCR HPCC. It runs the full
[funannotate](https://github.com/nextgenusfs/funannotate) workflow — genome
clean → repeat mask → RNA-seq fetch/train → predict → optional
antiSMASH/InterProScan/SignalP → annotate/update — plus a standalone
[EarlGrey](https://github.com/TobyBaril/EarlGrey) curated repeat-masking pipeline
([see below](#earlgrey-repeat-masking)).

The pipeline code lives under [`nextflow/`](nextflow/).

## Quick start

```bash
# graph/dry test (no SLURM, no tools needed)
cd nextflow
nextflow run funannotate.nf -c nextflow.config -profile test -stub-run

# real run on SLURM with environment modules (from the project root)
sbatch nextflow/run_annotate.sh --n_test 1
```

## Orthogonal profiles

Compose one option from each of three axes: `-profile <pipeline>,<executor>,<provisioning>`

| Axis | Options |
|---|---|
| **pipeline** | `annotate` · `earlgrey` · `test` / `stub` |
| **executor** | `slurm` · `local` |
| **provisioning** | `module` (default) · `pixi` · `singularity` |

```bash
nextflow run funannotate.nf -c nextflow.config -profile annotate,slurm,module -resume
nextflow run funannotate.nf -c nextflow.config -profile annotate,local,singularity -resume
```

The `run_annotate.sh` launcher honours `EXECUTOR=` and `PROVISION=` env vars
(default `slurm` / `module`).

Process scripts carry **no `module load`** — provisioning is supplied per process
`label` by the provisioning profile (`conf/provision_*.config`): a `beforeScript`
(module/pixi) or a `container` (singularity).

### Singularity images to build

Public biocontainers are used for edirect/antismash/interproscan/setup. Build
these and point at them with `--container_*` (defaults under
`/bigdata/stajichlab/shared/lib/singularity_cache`):
`funannotate`, `AAFTF` (genome_clean), the SRA multi-tool image, and
`signalp6-gpu`. `mariadb.sif` (PASA) already exists in shared lib.

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
`nextflow/tools/` (gitignored) rather than checked in (they are dynamically-linked,
platform-specific ELFs):

```bash
module load rust            # Rust toolchain, edition 2024+ (rust >= 1.85)
bash nextflow/scripts/build_tools.sh
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
nextflow/
  nextflow.config                 # manifest, shared params, profiles map, singularity block
  funannotate.nf                  # full annotation workflow (labeled processes)
  earlgrey_mask.nf                # standalone curated repeat masking (EarlGrey)
  conf/
    profile_annotate.config       # params + per-process resources
    provision_{module,pixi,singularity}.config
    profile_earlgrey.config
    test.config                   # self-contained stub profile
  lib/SampleUtils.groovy
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
