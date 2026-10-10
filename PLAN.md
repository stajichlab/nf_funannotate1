# Plan: Develop `nf_euk_genome_annotate` Nextflow framework

## Context

`nf_euk_genome_annotate/nextflow/` is currently a 3-file stub (`nextflow.config`,
`funannotate.nf` = predict only, `antismash.nf`). A complete, production-tested
fungal annotation pipeline already exists at
`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow` — a 1934-line
`funannotate.nf` (clean → mask → SRA fetch → RNA-seq train → predict →
antismash/interpro/signalp → annotate/update), an `earlgrey_mask.nf`, a layered
`nextflow.config` + `conf/profile_*.config`, `lib/SampleUtils.groovy`, `bin/`,
and `scripts/`.

Goal: migrate that **full** pipeline into `nf_euk_genome_annotate` while adding
what the reference lacks — **orthogonal, selectable provisioning** (module /
pixi / singularity) and **executor** (slurm / local) — and make it
**fungal-by-default but generic-capable** for any eukaryote. EarlGrey masking is
scaffolded but deferred. The reference hardcodes `module load X` inside every
process `script:` block, which prevents switching provisioning; the central
refactor is to lift tool provisioning out of script bodies into per-`label`
`beforeScript`/`container` directives supplied by a provisioning profile (the
reference already proves this pattern works — `conf/test.config` sets
`beforeScript = ':'` per label to disable module loads in stub runs).

Decisions confirmed with user: **full pipeline now**; **orthogonal profiles,
module default**; **fungal defaults but parameterized for generic euk**; **both
input models** (NCBI_ASM source dir by ASMID *and* direct local FASTA path).

## Target layout (`nf_euk_genome_annotate/nextflow/`)

```
nextflow.config                 # manifest, shared params, process defaults, singularity block, profiles map
funannotate.nf                  # FULL migrated annotation workflow (replaces current stub)
earlgrey_mask.nf                # migrated, deferred/optional (masking)
lib/SampleUtils.groovy          # copied verbatim from reference
conf/
  profile_annotate.config       # pipeline params + per-label resources (queue/cpus/mem/time/retry)
  provision_module.config       # per-label beforeScript = module loads (DEFAULT)
  provision_pixi.config         # per-label beforeScript = pixi env activation
  provision_singularity.config  # singularity.enabled=true + per-label container=
  profile_earlgrey.config       # masking params + resources (deferred)
  test.config                   # -stub-run: tiny resources, beforeScript=':', synthetic data
bin/                            # helpers that processes call: clean_genome_fa.py, setup_fcs_shm.sh, enforce_seqpair_readlen, fix_fastq_header_trinity, select_repeat_representatives.py
scripts/                        # tools run by hand or by build steps only (build_tools.sh, audit/backfill helpers)
pixi.toml                       # [workspace] + per-tool environments for pixi provisioning
run_annotate.sh                 # sbatch launcher
run_earlgrey.sh                 # sbatch launcher (deferred)
tests/data/                     # test_samples.csv + tiny inputs for stub
```

## Provisioning architecture (the core new capability)

Three **independent** profile axes, combined on the CLI, e.g.
`-profile annotate,slurm,module` (default) or `-profile annotate,local,singularity`:

1. **Executor** — `slurm` (executor name=slurm, queueSize, submitRateLimit,
   pollInterval, exitReadTimeout) / `local` (executor=local). Lifted from
   reference `nextflow.config` executor block.
2. **Provisioning** — `module` / `pixi` / `singularity`. Each sets the same set
   of `withLabel:` blocks but fills `beforeScript` (module/pixi) **or**
   `container` + `singularity.enabled=true` (singularity).
3. **Pipeline** — `annotate` (params + resources) / `earlgrey` / `test`.

**Refactor required in `funannotate.nf`:** every process gets a `label`, and all
`module load …` lines are **removed** from `script:` bodies. Provisioning is then
the profile's job. `beforeScript` and `script` are concatenated into one job
script by Nextflow, so module/pixi activations persist into the script.

Label → tool/module mapping (extracted from reference grep):

| label | processes | module(s) | singularity container |
|---|---|---|---|
| `funannotate` | MASKREPEAT_TANTAN_RUN (`funannotate mask -m tantan`), RNASEQ_PREPARE, FUNANNOTATE_TRAIN/PREDICT/ANNOTATE/UPDATE | miniconda3, funannotate (+fastp for train) | funannotate `.sif` (user builds) |
| `genome_clean` | GENOME_CLEAN (AAFTF `fcs_gx_purge` + `clean_script`; needs /dev/shm gxdb, highmem) | miniconda3, taxonkit, AAFTF | AAFTF `.sif` (user builds) |
| `setup` | SETUP_TAXONDB | miniconda3, taxonkit, AAFTF | AAFTF/taxonkit img |
| `edirect` | SRA_QUERY, SRA_QUERY_BATCH | ncbi_edirect | biocontainer entrez-direct |
| `sra` | SRA_FETCH, SRA_FETCH_SE | sratoolkit, parallel-fastq-dump, fastp, BBTools, workspace/scratch, seqkit, aria2 | multi-tool img (user builds) |
| `antismash` | ANTISMASH_RUN | miniconda3, antismash | biocontainer antismash |
| `interproscan` | INTERPROSCAN_RUN | interproscan | `docker.io/interpro/interproscan` (per ref interproscan6 profile) |
| `signalp` | SIGNALP_RUN | signalp/6-gpu (GPU) | signalp6 `.sif` (user builds) |

**Special case — PASA MariaDB service:** FUNANNOTATE_TRAIN/UPDATE start a
`mariadb.sif` singularity *instance* for PASA (gated by `params.pasa_mysql`).
This is a runtime service, not tool provisioning, so its `module load singularity
&& singularity instance start` stays inline in the script (guarded by
`params.pasa_mysql`), independent of the chosen provisioning axis. Noted as a
known coupling; `module load singularity` is available on HPCC in all modes.

`nextflow.config` keeps a `singularity {}` block (enabled=false by default;
autoMounts, cacheDir=`/bigdata/stajichlab/shared/lib/singularity_cache`,
`envWhitelist='SCRATCH'`); the `singularity` profile flips `enabled=true`.

## Migration steps

1. **Scaffold config layer.** Write `nextflow.config` (manifest, shared params:
   `samples`, `n_test`, `taxon`, `asmid`, `suppress`, debug; process defaults
   `shell=['/bin/bash','-l']`, `errorStrategy` retry, `maxRetries=2`,
   `cache='lenient'`, `clusterOptions='-N 1 -n 1'`; singularity block; profiles
   map wiring all three axes + test/stub). Base on reference `nextflow.config`.
2. **Copy `lib/SampleUtils.groovy` verbatim** — used for `makeSampleTag`.
3. **Migrate `funannotate.nf` full workflow** from reference (all 17 processes +
   the `workflow{}` orchestration incl. taxon/asmid/suppress filters, SRA branch
   routing, train/predict/annotate/update gating). Apply two transforms:
   - add `label '<x>'` to each process;
   - delete every `module load` from `script:` (keep the gated PASA mariadb
     `singularity instance` block).
   - Generalize hardcoded fungal assumptions into params already present in
     `profile_funannotate.config` (busco via `BUSCO_LINEAGE` column,
     `antismash_taxon`, `proteins`) — keep fungal values as defaults.
4. **Add dual input model** in the `workflow{}` channel builder: if a `GENOME`
   (or `GENOME_PATH`) column is present and non-empty, use that file directly;
   otherwise resolve `${params.source}/${ASMID}/${ASMID}_genomic.fna.gz` as the
   reference does. Both feed the same `GENOME_CLEAN` input tuple.
5. **Write the three `provision_*.config`** with matching `withLabel:` sets per
   the table above. `module` is the default and the only fully-populated one;
   `pixi` references `pixi.toml` envs; `singularity` lists known biocontainer
   URIs and placeholder paths for images the user will build.
6. **Write `profile_annotate.config`** = params block (from reference
   `profile_funannotate.config`) + per-label resources (queue/cpus/mem/time,
   retry escalation for SRA_FETCH/TRAIN, GPU `clusterOptions` for signalp) +
   per-profile `workDir=work/annotate` and `trace/report/timeline` paths. Drop
   the reference's node-pinning `-w h04,h05,h06` (skill: use queue+resources).
7. **Copy supporting scripts** referenced by processes into `bin/` (not `scripts/`;
   see the bin/ rule under Open items):
   `clean_genome_fa.py` (already expected by current stub), `setup_fcs_shm.sh`,
   `enforce_seqpair_readlen`, `fix_fastq_header_trinity`,
   `select_repeat_representatives.py` (for earlgrey).
8. **`pixi.toml`** — `[workspace]` (per global CLAUDE.md, not `[project]`) with
   conda environments per label where feasible (sra, antismash, interproscan,
   edirect, tantan); mark funannotate/signalp6-gpu as best-effort placeholders.
9. **`test.config` + `tests/data/`** — synthetic `test_samples.csv` (2 rows,
   both a local-FASTA row and an ASMID row), tiny FASTA, `beforeScript=':'` per
   label, `executor='local'`, tiny resources — for `-stub-run`.
10. **`run_annotate.sh`** sbatch launcher (load nextflow module, mkdir logs, run
    `-profile annotate,slurm,module -resume "$@"`), modeled on reference.
11. **EarlGrey (deferred):** copy `earlgrey_mask.nf` + `profile_earlgrey.config`
    + `run_earlgrey.sh` with the same label/provisioning refactor, but leave it
    untuned/untested this pass per user ("leave that for later"). It writes
    `input_clean_genomes/<asmid>.masked.fasta` consumed by the tantan storeDir.

## Reference files to reuse (read again at implementation time)

- `…/Fungi_BFD/nextflow/funannotate.nf` (lines 41–1461 processes, 1478–1934 workflow)
- `…/Fungi_BFD/nextflow/nextflow.config`, `conf/profile_funannotate.config`, `conf/test.config`, `conf/profile_interproscan6.config` (singularity wiring), `conf/profile_ANI.config` (singularity profile example)
- `…/Fungi_BFD/nextflow/lib/SampleUtils.groovy`, `earlgrey_mask.nf`, `run_funannotate.sh`, `run_earlgrey.sh`
- `…/Fungi_BFD/scripts/{clean_genome_fa.py?,setup_fcs_shm.sh,enforce_seqpair_readlen,fix_fastq_header_trinity,select_repeat_representatives.py}`

## Verification

1. **Graph/stub (login node, no SLURM):**
   `nextflow run nextflow/funannotate.nf -c nextflow/nextflow.config -profile annotate,local,module -stub-run --n_test 2`
   — confirms channel logic, dual input model, all process stubs fire.
2. **Provisioning axes parse:** repeat `-stub-run` with
   `-profile annotate,local,singularity` and `…,pixi` to confirm the config
   layers load and `withLabel` blocks resolve (stub skips real tools).
3. **Single real sample on SLURM:**
   `sbatch nextflow/run_annotate.sh --n_test 1 --run_repeatmasker true --run_sra_fetch false`
   then inspect `genome_annotation/<sample>/predict_results/*.gbk` and
   `logs/nextflow/annotate_report.html`.
4. **Confirm provisioning swap works for real:** run one sample with
   `-profile annotate,slurm,module` and (once an image exists) one process under
   `…,singularity`, checking the trace shows the container vs module path.

## Open items / notes

- **bin/ vs scripts/ rule.** Put any helper that a process calls or sources in
  `bin/`. Under `-profile singularity` the task runs in `apptainer exec --no-home`
  and Nextflow binds only `${projectDir}/bin` (plus the work dir and the
  configured mounts). A script in `scripts/` is not visible in the container, so
  the task fails with "No such file or directory" (seen 2026-10-10 with
  `setup_fcs_shm.sh` in GENOME_CLEAN_BATCH). Reference these helpers as
  `${projectDir}/bin/<name>`. Keep `scripts/` for tools that processes do not call.

- GENOME_CLEAN depends on FCS-GX (`/dev/shm/gxdb` via `setup_fcs_shm.sh`) and is
  highmem (450 GB) — keep its `highmem` queue resourcing; this is the heaviest
  prerequisite for the SLURM path and won't run under `-profile local` for real.
- `clean_genome_fa.py` (=`params.clean_script`) IS required (used in GENOME_CLEAN
  after AAFTF purge); it lives in `bin/` here (the reference repo keeps it in `scripts/`).
- Singularity images the user must build: funannotate, AAFTF, signalp6-gpu, the
  sra multi-tool image; mariadb (`mariadb.sif`) already exists in shared lib.
