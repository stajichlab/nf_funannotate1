# Changelog

All notable changes to `nf_funannotate1` are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Refit EVM weights (#42): `predict_evm_weights` = `augustus:1 hiq:3 genemark:2 snap:1 pasa:4`
  when the genome has a PASA training set, `predict_evm_weights_norna` =
  `augustus:2 hiq:5 genemark:4 snap:1` otherwise; appended to predict's single `-w` group,
  not applied with Prodigal evidence. +3.8 to +4.4 holdout F1 on 38 held-out RefSeq genomes
  (never worse), +9.5 without RNA-seq (BFD DECISIONS D126/D127). `''` keeps funannotate's
  weights.
- funannotate rc.4 train options, passed only when the container supports them:
  `train_pasa_alt_splice` (false), `train_pasa_remove_contained` (`off`),
  `train_pasa_max_isoforms` (0), and `train_pasa_fl_cache` (true): PASA's full-length list is
  cached next to the shared Trinity assembly (`<assembly>.pasa_fl_accs`, md5-checked) when that
  folder is writable, so strains sharing an assembly compute it once.
- Optional `BUSCO_SCORE_LINEAGE` samples column: the BUSCO dataset `BUSCO_COMPLETENESS` scores
  against, when it differs from `BUSCO_LINEAGE`. `BUSCO_LINEAGE` also names the funannotate
  `--busco_db` (e.g. `dikarya`, which BUSCO 6 does not ship as a dataset), so it cannot change
  without changing the gene models. Applied only to that process's input, so other tasks' cache
  hashes do not change. Rows without the column keep using `BUSCO_LINEAGE`.

### Fixed
- FUNANNOTATE_TRAIN/UPDATE (pasa_mysql): mariadbd readiness wait is `MARIADB_START_TIMEOUT`
  seconds (default 180; was 30), and `stop_mysqldb` sends SIGKILL if mariadbd is still running
  60 s after SIGTERM instead of waiting without a limit. In Fungi_BFD a ~78 s container start
  plus an unbounded `wait` on mariadbd hung a TRAIN task for 20 h (Fungi_BFD_runs DECISIONS
  D138); here a slow start failed the task and a mariadbd that ignored SIGTERM would hang it.
- FUNANNOTATE_PREDICT no longer lets a late attempt clobber a delivered result. Nextflow can mark an attempt
  failed while its pod keeps running; the retry delivered, then the old attempt's `sync_back` replaced the
  delivered `predict_results/` and deleted the BUSCO result. With node-local scratch, an attempt whose GBK target
  is newer than its own start now keeps that result and discards its output. The BUSCO invalidation is no longer
  a blind `rm -rf`: BUSCO_COMPLETENESS records the md5 of the proteins it scored (`<id>/proteins.md5`) and predict
  removes the result only when the delivered proteins differ (or there is no record, as before).
- Task ends no longer wait on node-wide disk flushes: FUNANNOTATE_PREDICT's final bare `sync` is now
  `sync <GBK>` (fsync of the delivered file only). On a loaded CephFS node the bare `sync` blocked in
  uninterruptible wait for over 3 h after a task had finished, holding its pod slot.
- BUSCO_COMPLETENESS keeps `hmmer_output/` and `busco_sequences/` (about 99% of its ~3,500 small files per
  genome, read by nothing downstream) as one `run_<lineage>/busco_run_dirs.tar.gz`; summaries,
  `full_table.tsv` and `missing_busco_list.tsv` stay plain files. The storeDir move is a handful of files, which
  also removes the "Directory not empty" failure when a retry met a half-moved directory.
- GENEMARK_RUN runs in node-local scratch under k8s/nrp (`genemark_local_scratch`, default false, true in
  `conf/executor_k8s.config`). gmes_petap.pl wrote ~2,400 files (740 MB) per genome in the task directory on the PVC
  (max 11,800 files, 1.5 GB): 39 of the 41 GB and 95% of the files in `work/` after the bfd_wave1 rc.4 pilot. Only
  the `.gtf`, `.mod` and stdout log are copied back (EXIT trap, every exit path). Not yet run on the cluster.
- FUNANNOTATE_PREDICT: the "Not enough gene models N to train Augustus" catch read the log from the
  local-scratch `RUNDIR`, which `copy_logs_back` had just removed, so it never matched. A tiny genome
  (e.g. a 456 kb assembly that passes pre-flight) failed the task and, after retries, aborted the whole
  run instead of being recorded in `predict_skipped_too_small.tsv`. The catch now reads the copied log.

### Changed
- `container_funannotate` is `ghcr.io/nextgenusfs/funannotate:v1.9.0-rc.6` (was v1.9.0-rc.5).
  rc.6 fixes two crashes this pipeline does not hit: predict's `AUGUSTUS_BASE` error when the
  Augustus config directory is not named `config` (`augustus_config` ends in `/config`), and
  annotate with `--genbank`/`--gff` and no `--table` (ANNOTATE uses `-i`). The unbuilt rc.5 conda
  manifests are renamed to `environments/conda/funannotate-1.9.0-rc.6{,-rust}.yml` (pip pin
  `v1.9.0-rc.6`); `conda_env` default is `funannotate-1.9.0-rc.6`. README and
  `tests/test_interproscan6.sh` use rc.6. The rc.6 image passed version, `check` and unit tests
  on UCR HPCC; the full `funannotate test` run was still in progress.
- `container_funannotate` is `ghcr.io/nextgenusfs/funannotate:v1.9.0-rc.5` (was v1.9.0-rc.4).
  README container references updated to match. Conda envs stay at `funannotate-1.9.0-rc.1`.
  Added conda manifests `environments/conda/funannotate-1.9.0-rc.5{,-rust}.yml` (pip pin
  `v1.9.0-rc.5`); `conda_env` default is `funannotate-1.9.0-rc.5`. `tests/test_interproscan6.sh`
  uses the rc.5 image. Not built or run with rc.5 yet.
  rc.4 turns PASA `--ALT_SPLICE` off by default (its reports are unused; > 2.5 h of a 4.6 h
  train on one genome) and logs transcript/PASA counts and EVM weight sources.
- FUNANNOTATE_TRAIN passes `--aligners minimap2 blat` by default (was `minimap2`): from rc.4,
  `minimap2` alone means no blat, and without blat C. neoformans H99 lost 1.8 holdout F1.
- NRP site config (`conf/site_nrp.config`): right-sized task requests from the bfd_wave1 pilot trace.
  FUNANNOTATE_PREDICT and GENEMARK_RUN ask for 4 cpus / 16 GB (was 8 / 32 GB), BUSCO_COMPLETENESS 2 / 4 GB,
  GENOME_CLEAN (skip_fcs) 1 / 2 GB, ASM_STATS 1 / 2 GB; a failed attempt escalates (8 cpus / 32 GB for
  predict and GeneMark). Cuts the reserved footprint of a 25-pod fan-out roughly in half. Other sites
  are unaffected.
- BFD equivalence for genomes without RNA-seq (stajichlab/nrp-deploy
  `bfd_wave1`): the no-RNA-seq predict path now matches the BFD pipeline's.
  - `container_funannotate` is `ghcr.io/nextgenusfs/funannotate:v1.9.0-rc.3`
    (was 1.9.0-rc.1); same image digest as BFD's rc.3 .sif.
  - FUNANNOTATE_PREDICT passes `--table <TRANSL_TABLE>` (was omitted, so CUG
    yeasts with table 12 were predicted with table 1).
  - Pre-flight guard: BFD's `asm_preflight_stats.py` (adds the `too_small` and
    `no_training_contigs` verdicts and the repeat-masked share) and BFD's
    thresholds (small < 9 Mb, fragmented N50 < 15 kb; were 8 Mb / 10 kb).
    GENEMARK_RUN and FUNANNOTATE_PREDICT skip any verdict other than `ok`
    (were: `small_fragmented` only); a `small_fragmented` genome with Prodigal
    evidence still runs predict.
  - Repeat-aware EVM: at >= 60% soft-masked, predict passes `--repeats2evm
    --evm-partition-interval 1500` and `-w snap:0`
    (`predict_evm_repeat_*` params).
  - Folder names: `SampleUtils.makeSampleTag` uses BFD's `cleanStrain` (`*`
    and shell characters in a strain). Checked against BFD's rule on all
    23,683 BFD samples.csv rows: 0 differences (1 before).
  - `predict_defline_first_word` (default true) makes the defline rewrite
    optional; BFD does not rewrite deflines.
  - predict_misc keeps `weights.evm.txt` and `final_training_models.gff3.gz`
    as BFD does, and no longer fails when `trnascan.no-overlaps.gff3` is absent.
  - `k8s/params_bfd_wave1.yaml`: predict-only BFD settings for NRP, including
    `augustus_config_source` for BFD's staged Augustus config.

### Fixed
- BUSCO_COMPLETENESS now scores the fresh gene set of a genome re-predicted in
  the same run. Its already-complete row no longer arrives first and wins the
  per-genome dedup; already-complete genomes are scored only if this run
  didn't predict them.
- InterProScan re-ran on every run after funannotate annotate had finished:
  annotate gzips `iprscan.xml` and deletes it, and the done-check looked only
  for the plain file. It now accepts `iprscan.xml` or `iprscan.xml.gz`.
- The singularity-axis InterProScan step could not work: it ran
  `interproscan.sh` (IPS5 flags) in `interpro/interproscan:6.0.0`, which has
  no `interproscan.sh`, java or nextflow.
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
- The InterProScan 5 module is now `INTERPROSCAN5_RUN`
  (`--interproscan_engine ips5`).
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
- InterProScan 6 for `--run_interpro` (default `--interproscan_engine ips6`).
  `INTERPROSCAN_RUN` launches the pinned IPS6 workflow (6.0.2.2) as a nested
  `nextflow run`, one task per genome, with IPS6's local executor inside the
  task allocation and its work dir on node-local `$SCRATCH`. It writes
  `annotate_misc/iprscan.xml.gz` (read directly by funannotate annotate) and
  `iprscan.tsv.gz`. The workflow checkout, images and InterPro data are
  shared, set up once by `scripts/setup_interproscan6.sh`. New
  `iprscan6_*` params; UCR paths in `conf/site_ucr_hpcc.config`. Licensed
  IPS6 apps stay off; `assets/interproscan6/licensed_ucr_hpcc.config` has
  correct UCR paths if they are wanted. Real-data test:
  `tests/test_interproscan6.sh` (Ordospora colligata OC4, 1,864 proteins:
  11 min on 16 CPUs, 8.6 GB peak RSS; funannotate's parser extracted
  InterPro terms for 1,348 proteins and GO terms for 1,163).
- InterPro data release pinned to 110.0 (`iprscan6_interpro`; UCR datadir
  `/srv/projects/db/interproscan/6.0.0/110.0`). 110.0 is the Matches API
  release. IPS6 6.0.2.2 does not compare the API release with the local
  data, so with older local data one output mixed two InterPro releases.
- `funannotate.nf` stops at startup when `--run_interpro` is on under the
  `k8s` / `nrp` profiles (skipped for `-stub-run`). Neither InterProScan
  engine can run in a pod, so each genome's task used to fail there instead.
- UCR note, found while testing: the `interproscan6` label uses a non-login
  shell, because a UCR login shell drops the module-loaded apptainer from
  PATH.
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
