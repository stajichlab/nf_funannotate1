# Changelog

All notable changes to `nf_funannotate1` are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Flattened pipeline to the repo root so it runs directly from GitHub:
  `nextflow run stajichlab/nf_funannotate1`. Set `manifest.mainScript`.
- `nextflow_schema.json` + nf-schema parameter validation and a schema-driven
  `--help`.
- Repository metadata: `LICENSE` (MIT), `CHANGELOG.md`, `CITATIONS.md`,
  `CODE_OF_CONDUCT.md`, and a GitHub Actions CI workflow (config parse + `-stub-run`).

### Changed
- `run_annotate.sh` / `run_earlgrey.sh` resolve the pipeline by project name
  (`PIPELINE` / `REVISION` overrides) instead of the script's own path, so they are
  safe under `sbatch` (which copies the script to a spool dir).

## [0.1.0]

### Added
- Initial framework: full funannotate workflow (clean → mask → RNA-seq fetch/train
  → predict → optional antiSMASH/InterProScan/SignalP → annotate/update) plus a
  standalone EarlGrey repeat-masking pipeline, with orthogonal
  pipeline/executor/provisioning profiles for the UCR HPCC.
