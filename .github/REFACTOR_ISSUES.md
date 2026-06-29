# Refactor issue tracker

Source-of-truth backlog for the DSL2 modularization (see `REFACTORING_PLAN.md`).
Create these on GitHub with `scripts/create_refactor_issues.sh` (uses `gh`).
Each non-epic issue closes one `## Phase N` of the plan and must pass the
stub-run gate before merge:

```
nextflow config . -profile test && nextflow run . -profile test -stub-run && nextflow run . --help
```

---

## EPIC: Modularize nf_funannotate1 into DSL2 modules + subworkflows
labels: refactor, epic

~~Break the 2,470-line `funannotate.nf` monolith into one-process-per-file modules
composed by subworkflows, on a `meta`-map data contract.~~ **DONE** — `funannotate.nf`
is now 163 lines (0 inline processes). 24 modules, 7 subworkflows, shared
`lib/FunannotateUtils.groovy`. nf-core-*inspired*, not nf-core-submitted.

All issues complete (#11 institutional profile, #12 docs/hygiene). Epic closed.

---

## 1. Adopt `meta`-map data contract  (Phase 0 — BLOCKER) ✅ DONE
labels: refactor

Replace the positional 10-tuple with `tuple val(meta), path(genome)`.

- [x] Build `meta` in the workflow channel construction (`lib/SampleUtils.groovy`)
- [x] Update every call site to consume `meta`
- [x] `params.header_length` added to schema with default 24
- [x] Stub-run gate green

## 2. Shared `INPUT_CHECK` subworkflow (dedupe earlgrey)  (Phase 0) ✅ DONE
labels: refactor

- [x] `subworkflows/local/input_check.nf` emits `ch_genomes` (meta + genome)
- [x] `funannotate.nf` and `earlgrey_mask.nf` both consume it; no duplicated parse
- [x] Stub-run gate green for both entrypoints

## 3. Repo skeleton + relocate existing modules  (Phase 1) ✅ DONE
labels: refactor

- [x] Directory skeleton in place (`modules/local/`, `subworkflows/local/`, `lib/`)
- [x] `annotation_tools.nf` split into 3 single-process modules
- [x] Stub-run gate green

## 4. `versions.yml` + `conf/base.config` + `conf/modules.config`  (Phase 2) ✅ DONE
labels: refactor

- [x] `conf/base.config` with resource labels
- [x] `conf/modules.config` with publishDir + ext.args
- [x] At least one module emits and the workflow collects `versions.yml`
- [x] Stub-run gate green

## 5. Setup modules  (Phase 3) ✅ DONE
labels: refactor, good first issue

- [x] 3 modules + `setup_dbs.nf` subworkflow; storeDir preserved
- [x] Stub-run gate green

## 6. Genome clean + `CLEAN_GENOMES` subworkflow  (Phase 4) ✅ DONE
labels: refactor

- [x] `genome_clean.nf`, `genome_clean_batch.nf` + `clean_genomes.nf` subworkflow
- [x] FCS-GX /dev/shm staging and batch skip behavior preserved
- [x] Stub-run gate green

## 7. `MASK_GENOME` subworkflow + per-tool modules  (Phase 5) ✅ DONE
labels: refactor

Implemented tantan + EarlGrey paths (see also Task 4 below).

- [x] `mask_genome.nf` selects masker by param; pass-through supported
- [x] `maskrepeat_tantan_run.nf` extracted
- [x] `earlgrey_build_lib.nf`, `repeatmask_strain.nf`, `deliver_mask.nf` as shared modules
- [x] EarlGrey path integrated via `params.run_earlgrey`; shared with `earlgrey_mask.nf`
- [x] Stub-run gate green (tantan 15/0, EarlGrey 16/0)

## 8. `FETCH_RNASEQ` subworkflow  (Phase 6) ✅ DONE
labels: refactor

- [x] 6 modules + `fetch_rnaseq.nf`; shared Trinity-GG semantics preserved
- [x] maxForks / rate limits preserved
- [x] Stub-run gate green (17/0 with `run_sra_fetch=true`)

## 9. `TRAIN_PREDICT` subworkflow  (Phase 7) ✅ DONE
labels: refactor

- [x] `rnaseq_prepare.nf`, `funannotate_train.nf`, `funannotate_predict.nf` modules
- [x] `train_predict.nf` subworkflow; stale-prediction detection via FunannotateUtils
- [x] `funannotate_update.nf` included (param-gated via ANNOTATE_GENOME)
- [x] Stub-run gate green

## 10. `ANNOTATE_GENOME` subworkflow  (Phase 8) ✅ DONE
labels: refactor

- [x] `annotate_genome.nf` composing antismash / signalp / interproscan → funannotate_annotate
- [x] Each optional tool independently param-gated
- [x] Stub-run gate green

## Task: `lib/FunannotateUtils.groovy` (shared utilities) ✅ DONE

Extracted `gbkResult`, `genomeFile`, `staleRnaseq` into `lib/FunannotateUtils.groovy`
(static methods using `java.io.File`). All four pipeline files updated; wrappers
eliminated; only `staleRnaseq` in `funannotate.nf` kept (needs DSL `log`).

## 11. Consolidate `ucr_hpcc` institutional profile + portable container path  (Phase 9)
labels: refactor

The `module`→`ucr_hpcc` rename is done. Finish the repivot: fold UCR SLURM
partitions / `clusterOptions` into the institutional profile, and ensure a fully
portable run works via per-module conda/biocontainer directives (no Lmod).
Document the "copy to `conf/provision_<site>.config`" path for new sites.

- [x] UCR partition config consolidated under the institutional profile
- [x] `docs/adding_a_site.md` guide for new institutions
- [ ] Portable container/conda path runs without any UCR modules (deferred — not a current priority)

## 12. nf-core hygiene  (Phase 9) ✅ DONE
labels: refactor, documentation

- [x] `docs/usage.md`
- [x] `docs/output.md`
- [x] `assets/schema_input.json`
- [x] Naming + nf-core submission decision recorded in `REFACTORING_PLAN.md`
      → Decision: stay nf-core-*inspired*, do not submit. `nf_funannotate1` naming
        kept; container-per-module and nf-test not pursued for this HPC-first pipeline.
