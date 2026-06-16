# nf_euk_genome_annotate

A Nextflow DSL2 framework for **eukaryotic genome annotation** (fungal defaults,
generic-capable) on the UCR HPCC. It runs the full
[funannotate](https://github.com/nextgenusfs/funannotate) workflow — genome
clean → repeat mask → RNA-seq fetch/train → predict → optional
antiSMASH/InterProScan/SignalP → annotate/update — plus a deferred
[EarlGrey](https://github.com/TobyBaril/EarlGrey) repeat-masking pipeline.

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
| **pipeline** | `annotate` · `earlgrey` (deferred) · `test` / `stub` |
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
  earlgrey_mask.nf                # deferred repeat masking
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
Post-predict steps (antismash/interpro/signalp/annotate) fire on a resumed run
once `predict_results/*.gbk` exist (2-pass production semantics).
