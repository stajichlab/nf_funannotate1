# Changelog

All notable changes to `nf_funannotate1` are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- FUNANNOTATE_PREDICT's own "already complete and current" check now also
  treats the genome's PASA training output, and for reuse siblings the shared
  ab-initio store, as evidence newer than the GBK. Before, it checked only reads
  and Trinity, so a genome the pipeline had flagged as stale exited as a no-op
  and kept its old annotation.
- ANI reuse: re-running the representative pick (e.g. after the sample sheet
  changes) no longer drops GeneMark from the shared store. Its inline backfill
  passes no `.mod`, so `backfill_abinitio_params.py` now reuses the store's
  existing GeneMark model when the store was built by the same representative.
  Before, the rebuilt store lacked GeneMark, its content hash changed, and every
  sibling looked stale.
- ANI reuse: sibling GeneMark and predict tasks now carry a fingerprint of the
  species' shared ab-initio store (`FunannotateUtils.sharedParamsFingerprint`:
  provenance `content_hash`, else size+mtime). Previously a rebuilt store (new
  representative) had the same path, so `-resume` served siblings' old
  GeneMark and predict from Nextflow's cache even when `staleSharedParams`
  flagged them. Seen on the Bd run: 9 pilot strains kept annotations from the
  pilot representative.

### Changed
- NRP fair use (https://nrp.ai/documentation/userdocs/running/jobs/, .../cpu-only/):
  Nextflow now runs as a per-run **Job** (`k8s/run/`) instead of an idle
  `sleep infinity` head Deployment, which NRP prohibits. Staging and the Rust
  helper build are finite Jobs (`k8s/tools/s3-sync`, `build-rust-tools`);
  `k8s/tools/shell` is a 1 h setup/inspection pod. Task pods get
  `priorityClassName: opportunistic` and avoid GPU nodes (`site_nrp.config`).
- k8s run Job: the Nextflow resume cache now runs on pod-local disk
  (`NXF_CACHE_DIR`), with LevelDB memory-mapping off (`-Dleveldb.mmap=false`).
  It is snapshotted to `RUN_DIR/.nextflow-snapshots` atomically (write to a
  temporary directory, fsync, rename) every 5 min and when Nextflow exits, and
  restored on the next launch (`k8s/run/nf-run.sh`). Previously a node lost
  mid-write left the cache on CephFS corrupt and every retry failed.
  A relaunch waits for the previous pod's heartbeat (`RUN_DIR/.nf-run.heartbeat`)
  to go stale, so a `kubectl delete` + `apply` can't copy a cache that is still
  being written. A restored snapshot that won't open is renamed `bad-*`, and the
  Job's retry falls back to the previous one.

### Added
- `train_cleanup` (default off): once a genome's training resolves, and again
  after its predict succeeds (which covers genomes trained earlier),
  `bin/train_cleanup.sh` deletes funannotate-train intermediates that predict
  never reads: `getBestModel/`, the GMAP index, the seqclean copies of Trinity
  and `pasa.step1.gff3`, plus `pasa/` when `run_update` is off. Anything a
  `training/` symlink points at is kept (e.g. `funannotate_train.trinity-GG.fasta`
  -> `trinity.fasta`). On a Bd training dir this took 1.6 GB down to 129 MB.
- Kubernetes execution (test mode): `-profile annotate,k8s` and
  `-profile annotate,nrp` (NRP Nautilus), via `conf/executor_k8s.config` and
  `conf/site_nrp.config`, plus `k8s/` manifests (PVC, head pod), smoke-test
  params and a README. SignalP/DeepTMHMM/antiSMASH and the FCS-GX purge are
  not provisioned on k8s yet.
- Kubernetes setup files in `k8s/`: `head-deployment.yaml` (the head runs as
  a Deployment, since NRP ends controller-less pods after 6 h), `pvc.yaml`,
  `rbac.yaml` (adds `batch/jobs`; not applied yet), `build-tools-pod.yaml`,
  `s3-stage-pod.yaml`, and `params_nrp_test.yaml` for the D. hansenii CDA1
  test. `site_nrp.config` caps tasks at 16 cpus / 32 GB (NRP's limit for pods
  without a controller). Tested end to end on NRP: clean, mask, SRA fetch,
  Trinity/PASA train, GeneMark, predict, BUSCO completeness.
- `params.mariadb_setup_in_image`: SETUP_MARIADB_DATADIR can run
  mariadb-install-db inside its own image (no nested apptainer).
- `params.predict_local_scratch` (on for k8s/nrp, off elsewhere):
  FUNANNOTATE_PREDICT works in node-local `$TMPDIR` and copies only the final
  pruned tree to `params.target`. BUSCO/Augustus training writes thousands of
  small files, which crawled on CephFS.
- `funannotate.nf` stops at startup, with build instructions, when
  `--run_sra_fetch` is on and `fix_fastq_header_trinity` /
  `enforce_seqpair_readlen` are missing (skipped for `-stub-run`).

### Changed
- Container profiles (`singularity`, `k8s`, `nrp`) use `sra_tools:1.4.0`,
  which ships the two Rust read helpers, and call them by name. No
  per-checkout `scripts/build_tools.sh` run is needed there; host-tool
  profiles keep the `tools/bin` paths.
- `GENEMARK_RUN` falls back to ES self-training when GeneMark-ET produces no
  model, instead of failing the genome. ET cannot train on intron-poor
  genomes (e.g. Saccharomycetes).
- `slurm` + `singularity` profiles are now site-neutral for SignalP and
  DeepTMHMM. `conf/provision_singularity.config` no longer sets the UCR
  partitions (`short_gpu`, `epyc`) or `--exclude=gpu13,gpu14`. The generic
  `slurm` profile asks for `--gres=gpu:1` when `signalp_gpu` /
  `deeptmhmm_gpu` is true, with no partition. The UCR partitions and node
  excludes moved to `conf/provision_ucr_hpcc.config` and
  `conf/site_ucr_hpcc_singularity.config`.
- UCR: `SIGNALP_RUN` on the Lmod axis (no container) keeps the `gpu`
  queue. In a container (`-profile ...,ucr_hpcc,singularity`) it goes to
  `short_gpu` with the gpu13/gpu14 exclude in GPU mode, or `epyc` in CPU
  mode, as before. `DEEPTMHMM_ANNOTATION` (container only) gets the same
  container routing under `ucr_hpcc`.

### Fixed
- SRA queries for species with a blank `NCBI_TAXONID` searched
  `txid[Organism:noexp]`, matched nothing, and trained without RNA-seq. They
  now fall back to the species name. `samples.csv` gets CDA1's taxon ID (4959).
- `BUSCO_COMPLETENESS` kept a stale result after a genome was re-predicted
  (storeDir skips on the directory existing). Predict now removes it after
  writing new results, and already-predicted genomes are routed to
  BUSCO_COMPLETENESS so a missing result is filled in.
- PASA's in-image MariaDB failed to start when tasks run as root (k8s):
  `mariadbd` now gets `--user=root` only in that case.
- BUSCO lineage docs and `samples.csv` use odb10: predict's BUSCO training
  needs `lengths_cutoff`, which odb12 datasets lack.
- `BUSCO_COMPLETENESS` did not run for genomes whose `ASMID` differs from
  the `SPECIES_STRAIN` tag. `train_predict.nf` looked for the proteins at
  `genome_annotation/<ASMID>/predict_results/<ASMID>.proteins.fa`, but
  predict writes them under `<SPECIES_STRAIN>/`, so the `exists()` filter
  dropped those genomes. It now uses the tag. The output also moved from
  `genome_annotation/<ASMID>/busco_completeness/` to
  `genome_annotation/<SPECIES_STRAIN>/busco_completeness/<SPECIES_STRAIN>/`,
  next to `predict_results/`.
- `docs/output.md`: output directories and file prefixes under
  `genome_annotation/` and `genome_annotation_training/` are
  `<SPECIES>_<STRAIN>`, not `<ASMID>`. Corrected the RNA-seq paths
  (reads in `rnaseq_reads/`, species-only tag, `.sra_query.csv`), the
  antiSMASH/InterProScan/SignalP/DeepTMHMM locations, and the trace file name.
  Added missing outputs (SRA manifests, GeneMark model, BUSCO, ANI, EarlGrey
  `strains/`, database caches) and corrected the storeDir skip rule.
- `BUSCO_COMPLETENESS` passed `-l <busco_lineages>/<lineage>`, which skips
  the `lineages/` level of a BUSCO download tree. At UCR HPCC every call
  failed (`/srv/projects/db/BUSCO/v10//fungi_odb10 does not exist`), and
  after 3 attempts the global `errorStrategy 'finish'` stopped the run
  from submitting new tasks. It now uses `-l <lineage> --offline
  --download_path <busco_lineages>`, the same convention as `BUSCO_GENOME`.
- `FUNANNOTATE_TRAIN`: funannotate 1.9.0-rc.3 stops train with exit 3 when
  too few sampled RNA-seq reads map to the genome ("RNA-seq concordance
  gate FAILED"). The module treated this as an infra failure and retried
  it, then tripped `errorStrategy 'finish'`. It now writes the
  `.pasa_train_failed` marker (`pasa_tier` = `rnaseq_gate`) and exits 0,
  so predict runs ab initio for that strain. Exit 3 without the gate
  message still hard-fails. Test: `tests/test_train_rnaseq_gate.sh`.
- `tests/test_train_retry_cleanup.sh` did not parse the extracted block:
  it unescaped Groovy `\$` but not `\\`, so line continuations broke the
  rendered shell. It passes again.

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
