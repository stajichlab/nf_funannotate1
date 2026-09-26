# Changelog

All notable changes to `nf_funannotate1` are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2026-09-25

Changes since the last changelog update (2026-06-26, `1108471`) up to
2026-09-25 (`bdee611`).

### Added
- Flattened pipeline to the repo root so it runs directly from GitHub:
  `nextflow run stajichlab/nf_funannotate1`. Set `manifest.mainScript`.
- `nextflow_schema.json` + nf-schema parameter validation and a schema-driven
  `--help`.
- Repository metadata: `LICENSE` (MIT), `CHANGELOG.md`, `CITATIONS.md`,
  `CODE_OF_CONDUCT.md`, and a GitHub Actions CI workflow (config parse + `-stub-run`).
- Modular DSL2 layout (#2–#10): processes in `modules/local/`; subworkflows
  `INPUT_CHECK`, `SETUP_DBS`, `CLEAN_GENOMES`, `MASK_GENOME`, `FETCH_RNASEQ`,
  `TRAIN_PREDICT`, `ANNOTATE_GENOME` in `subworkflows/local/`; per-sample meta
  map (`SampleUtils.makeMeta`); shared helpers in `lib/FunannotateUtils.groovy`;
  `versions.yml`, `conf/base.config`, `conf/modules.config`.
- `ASM_STATS` module; `SELECT_REPS` made optional.
- EarlGrey masking as an alternative path in `MASK_GENOME`.
- In-run ANI representative-strain selection (BUSCO + skani) and shared
  ab-initio predict reuse.
- GeneMark-ES/ET run from the public `teambraker/braker3` container; standalone
  GeneMark sidecar entrypoint (`genemark_sidecar.nf`, `--genemark_sidecar_dir`)
  with RNA-seq intron hints.
- SignalP 6 and DeepTMHMM steps, with optional GPU.
- `SETUP_ANTISMASH_DB` module.
- Optional Prodigal ab-initio pass-through into `funannotate predict`,
  lineage-gated (`prodigal_lineages`), emitting a gene/mRNA/CDS hierarchy so
  EVM uses it.
- Conda provisioning axis: `conf/provision_conda.config`, env manifests under
  `environments/conda/`, shared aux env `nf_funannotate1-aux`, env root from
  `$CONDA_ENVS_ROOT`.
- funannotate 1.8 release seam: `stop_after_trinity` and `conf/release_1_8.config`.
- Containerized PASA MariaDB backend (`--pasa_mysql`) with an auto-built seed
  datadir (`SETUP_MARIADB_DATADIR`).
- Site override seam: `conf/site_template.config`, `docs/adding_a_site.md`,
  UCR site configs (including the singularity axis and an opt-in preempt config).
- PASA tier thresholds with graceful degrade to ab-initio-only training;
  `--pasa_aligners` to state the PASA aligner list explicitly.
- BUSCO completeness on the final gene set.
- `RNASEQ_PREPARE` keeps the StringTie GTF and splice-junction BED (optional outputs).
- `scripts/audit_rnaseq_community_runs.py` to audit fetched SRA runs for
  community (metatranscriptome) sources.
- Trace fields: attempt, cpus, memory, queue, hostname, submit/start/complete
  timestamps, peak_rss.
- Tests: regression tests for the PASA/TransDecoder PATH and train
  retry-cache-poisoning bugs; `tests/submit_train_pasa_mysql.sh` and
  `tests/check_train_pasa_mysql.sh` (real `FUNANNOTATE_TRAIN` run with PASA on
  MariaDB, singularity and conda arms, one `FUNANNOTATE_VERSION` knob).
- README: quick start for a new site (local genomes + RNA-seq + Singularity),
  list of public vs. build-it-yourself images, PASA MySQL (MariaDB) section.

### Changed
- `run_annotate.sh` / `run_earlgrey.sh` resolve the pipeline by project name
  (`PIPELINE` / `REVISION` overrides) instead of the script's own path, so they are
  safe under `sbatch` (which copies the script to a spool dir).
- Default funannotate is now 1.9.0-rc.1: container
  `ghcr.io/nextgenusfs/funannotate:1.9.0-rc.1`, conda env
  `funannotate-1.9.0-rc.1` (previous defaults moved through 1.9.0-beta.10,
  beta.11 and beta.12).
- UCR HPCC provisioning consolidated into the `ucr_hpcc` profile; hardcoded UCR
  paths removed from core config.
- `container_sra` defaults to the public CI-built image
  `ghcr.io/hyphaltip/sra_tools_container/sra_tools:1.3.1` (was a locally built
  `<sif_dir>/sra_tools.sif`).
- `container_mariadb` defaults to `container_funannotate`, which bundles MariaDB
  (11.8.6 in rc.1); `--container_mariadb` still overrides. The
  `container_funannotate` pin moved to `nextflow.config` so every provisioning
  axis has it. A `docker://` URI used by a direct `apptainer exec` resolves to
  Nextflow's cache file in `sif_dir` and is pulled at most once.
- `SETUP_MARIADB_DATADIR` runs on the host (`container = false`) at every site,
  not only UCR, and requests 32 GB for the one-time image pull.
- Apptainer is the preferred container engine; SLURM `clusterOptions` are
  scoped to the `slurm` profile.
- SRA query skips community-metatranscriptome runs.
- UCR HPCC: `RNASEQ_PREPARE` and `FUNANNOTATE_TRAIN` exclude nodes without
  AVX2/BMI2 (`c01-c30`); `FUNANNOTATE_TRAIN` also excludes `h01-h06`.
- A global queue setting lets jobs move between queues without false failure
  detection.
- nf-schema 2.7.2, `validation.lenientMode`, validation warns instead of failing.

### Fixed
- PASA MariaDB: per-task instance on node-local `$TMPDIR`; self-contained
  setup; readiness polling instead of a fixed sleep; sidecar now starts the
  server when the image's startscript does not (funannotate image);
  `SETUP_MARIADB_DATADIR` points `TMPDIR` inside the task workdir (MariaDB 11.8
  "Read-only file system"); `mysql_install_db` fallback for MariaDB 10.3.
- `SETUP_MARIADB_DATADIR` failed with "no apptainer/singularity binary on PATH"
  outside UCR (it ran inside the gnu-wget `setup` image).
- `FUNANNOTATE_TRAIN` retries: clear `trinity_gg/`, `hisat2/` and `pasa/`
  checkpoints on retry.
- `FUNANNOTATE_TRAIN`: TransDecoder PATH, PASA log preservation, degrade regex,
  evidence gates; `--aligners` on all five train branches; `gmap` alongside
  `minimap2` for 1.8.17.
- Predict re-runs when training evidence changes; `staleGenome()` ported from BFD.
- Conda axis: tasks run in a plain shell so activation survives; `JAVA_HOME`
  unset so the env's java wins; perl deps and `mysql-libs` pins;
  `perl-mce` added to `funannotate-1.9.0-rc.1` (GeneMark-ES/ET needs
  `MCE::Mutex`); 1.8.17 env pins (salmon, goatools/r-base, transdecoder <6).
- SRA: hardened `SRA_QUERY`/`SRA_QUERY_BATCH` against runaway efetch and stale
  metadata hits; bbnorm Java heap capped to the cgroup; `SRA_FETCH` retry
  budget; edirect Perl `Time::HiRes` shim.
- `RNASEQ_PREPARE` Trinity workdir.
- FCS-GX database staging from `/srv` with fan-out copy; `/dev/shm` copy no
  longer fails on chgrp.
- Clean step normalizes FASTA deflines (fixed a funannotate EVM crash).
- Singularity axis: head-job memory and apptainer for image pulls; apptainer
  on PATH per label; login-shell PATH clobber; tantan sandbox-extract disk overflow.
- `ASM_STATS` uses pure awk/gzip and the correct ASMID column.
- `SCRATCH`/`TMPDIR` export strip on containerized tasks, plus a runtime guard
  that falls back to the task workdir.
- nf-schema incompatibility with Nextflow 24.04.0; empty `withName` block
  rejected by Nextflow 25.10.0.

## [0.1.0]

### Added
- Initial framework: full funannotate workflow (clean → mask → RNA-seq fetch/train
  → predict → optional antiSMASH/InterProScan/SignalP → annotate/update) plus a
  standalone EarlGrey repeat-masking pipeline, with orthogonal
  pipeline/executor/provisioning profiles for the UCR HPCC.
