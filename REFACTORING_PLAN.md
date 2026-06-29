# nf_funannotate1 — Modularization Plan (authoritative)

> This is the **single source of truth** for the DSL2 modularization effort. It
> supersedes the earlier `IMPLEMENTATION_SUMMARY.md` and `MODULE_STRUCTURE.txt`
> (removed). Progress is tracked in GitHub issues — see
> `.github/REFACTOR_ISSUES.md` and `scripts/create_refactor_issues.sh`.

## Where we actually are

- `funannotate.nf`: **2,470-line monolith, 20 inline processes.**
- `earlgrey_mask.nf`: separate 394-line pipeline that **duplicates** the
  samplesheet parse + taxonomy/asmid/suppress filtering from `funannotate.nf`.
- Extracted so far: `modules/asm_stats.nf`, `modules/annotation_tools.nf` (2 of 20).
- Strong base already present: nf-schema + `nextflow_schema.json`, orthogonal
  profiles, manifest, CITATIONS/COC/CHANGELOG/LICENSE, stub-run CI.

The goal is an **nf-core-*inspired*** layout (not nf-core submission): adopt the
parts that are pure engineering wins, skip the parts that fight our HPC reality.

---

## Principle 0 — the data contract (do this FIRST; blocks everything else)

Every process is currently wired with a fragile positional 10-tuple:

```groovy
tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, gz, taxonid)
```

Adding an 11th field touches every process. **Replace it with a `meta` map**, the
standard DSL2 idiom. Genome travels as a separate `path`:

```groovy
// Canonical channel element: tuple val(meta), path(genome)
meta = [
    id          : out,          // unique sample tag — used for tag{} and file naming
    asmid       : asmid,
    species     : species,
    strain      : strain,
    locustag    : locustag,
    busco       : busco,        // BUSCO_LINEAGE
    transl_table: transl_table, // default '1'
    taxonid     : taxonid,
]
```

Rules:
- `meta.id` is the **only** field used for naming/`tag`; everything else is payload.
- `header_length` (constant 24) becomes `params.header_length`, **not** a meta field.
- Build `meta` once, in the `INPUT_CHECK` subworkflow (below). No process
  re-parses the samplesheet.
- A module's `input:` declares `tuple val(meta), val(genome)` (genome is an
  absolute-path **string**, kept as `val` so the networked FS isn't re-staged)
  and never positionally unpacks fields it doesn't use.

Until `meta` is adopted, **do not extract more modules** — every module written
against the old tuple is rework.

### Groundwork landed (this is the only safe sub-step)

- `SampleUtils.makeMeta(row)` (`lib/SampleUtils.groovy`) is the **single
  authoritative definition** of `meta`, reproducing the current `jobs`-channel
  cleaning field-for-field so wiring it in is behaviour-preserving. Not yet called.
- `params.header_length` (default 24) added to `nextflow.config` + schema. Still
  threaded through the tuple for now; the conversion removes it from the tuple.

### Conversion recipe (the atomic change, not yet done)

The rest of #3 is **atomic** — source channel, ~40 workflow channel ops, and all
10 carrying processes move together. Recipe:

1. **Source channels** (`jobs` and `postpredict` maps): replace the per-field
   `def`s with `def meta = SampleUtils.makeMeta(row)` and emit `tuple(meta, gz)`
   (jobs) / `meta` (postpredict). `header_length` comes from `params`.
2. **Processes** (shim to keep script bodies intact): change `input:`/`output:`
   tuples to `tuple val(meta), val(genome)` (+ reads paths where present), add
   `tag "${meta.id}"`, and at the top of `script:` add the alias block
   `def out = meta.id; def asmid = meta.asmid; def species = meta.species; …`
   so every existing `${out}`/`${asmid}` interpolation still resolves.
3. **Workflow ops** — translate the position-coupled patterns:
   - `.map { out, asmid, … -> tuple(out, asmid, …) }` re-threads → `.map { meta, genome -> … }`
   - index access: `it[8].exists()` (genome) → `genome.exists()`; reads
     `it[10]/it[12]` → name them in the destructure.
   - slice sentinels: `row[0..8]` / `row[0..-3]` (drop combine/gate tails) →
     destructure `(meta, genome, _sentinel)` explicitly.
   - species-keyed `groupTuple`/`combine`/`join` for RNA-seq → key on
     `meta.species` (compute `species_tag` from `meta.species`), carry `meta`.
   - the reduced 8-tuple in the annotate phase collapses to a single `meta`.
4. **Validate**: `‑profile test ‑stub-run` **and** `‑profile test ‑stub-run
   --run_sra_fetch true` (the default stub skips the RNA-seq subgraph). Stub
   proves graph wiring only — a real-data HPCC run confirms semantics (grouping
   picks the right representative, reads join to the right strain).

---

## Target architecture

```
main.nf                              # thin entrypoint: parse args, call workflow
workflows/
    funannotate.nf                   # wires the subworkflows (was the monolith)
    earlgrey.nf                      # curated-mask entry, reuses shared subworkflows
subworkflows/local/
    input_check.nf                   # samplesheet -> meta channel + taxon/asmid/suppress filters  (SHARED)
    setup_dbs.nf                     # SETUP_TAXONDB / FUNANNOTATE_DB / AUGUSTUS_CONFIG gating
    prepare_genome.nf                # clean -> asm_stats -> mask
    mask.nf                          # selects ONE masker module by params.mask_tool
    rnaseq.nf                        # sra_query -> sra_fetch(_se) -> rnaseq_prepare
    predict.nf                       # train -> predict -> (update)
    annotate.nf                      # antismash|signalp|interpro -> funannotate annotate
modules/local/                       # ONE process per file
    setup_taxondb.nf  setup_funannotate_db.nf  setup_augustus_config.nf
    genome_clean.nf   genome_clean_batch.nf    asm_stats.nf
    mask_tantan.nf    mask_repeatmodeler.nf    mask_repeatmasker.nf  mask_earlgrey.nf
    sra_query.nf      sra_query_batch.nf       collect_sra_query.nf  write_empty_reads.nf
    sra_fetch.nf      sra_fetch_se.nf          rnaseq_prepare.nf
    funannotate_train.nf  funannotate_predict.nf  funannotate_update.nf
    antismash.nf      interproscan.nf          signalp.nf            funannotate_annotate.nf
    select_reps.nf                   # earlgrey representative selection
conf/
    base.config                      # resources by label (process_low/medium/high/...)
    modules.config                   # per-process publishDir + ext.args (nf-core idiom)
```

### Why this and not the old plan

- **One tool = one process = one module file.** The old plan bundled 3 processes
  per file (`sra_query.nf`, `annotation_tools.nf`, `databases.nf`). That is the
  *opposite* of the convention and kills reuse. Group sequences with
  **subworkflows**, which the old plan never mentioned.
- **Masking is a subworkflow, not a mega-process.** A single process with an
  `if/else` over NONE/TANTAN/REPEATMODELER/REPEATMASKER/EARLGREY is an
  anti-pattern. Each masker is its own module; the *selection* lives in
  `subworkflows/local/mask.nf`.
- **Kill the funannotate/earlgrey duplication.** Both entrypoints parse the same
  samplesheet and apply the same filters. Extract that into `input_check.nf` once
  and call it from both. EarlGrey is a whole pipeline (SELECT_REPS + asm_stats +
  representative-per-species), not a masking flavor — it reuses shared
  subworkflows rather than being folded into one module.

---

## Per-process extraction checklist (the stub-run gate)

Apply to **one process per commit/PR**. The monolith must stay runnable at every
commit.

For process `P`:

- [ ] Create `modules/local/<p>.nf` with `process P { ... }`.
- [ ] `input:` uses `tuple val(meta), path(...)` (no positional field unpacking).
- [ ] Add `tag "${meta.id}"` and a resource `label` (`process_low|medium|high`).
- [ ] Keep the existing `stub:` block; keep `storeDir`/`publishDir` behavior.
- [ ] Emit a version: `path "versions.yml", emit: versions` + a `cat <<-END_VERSIONS`
      block capturing the tool version. (Foundational — see Issue 3.)
- [ ] Move the process's resource/`withName` block out of
      `conf/profile_annotate.config` into `conf/modules.config` (name unchanged,
      so existing `withName:` selectors keep matching).
- [ ] In the workflow, replace the inline `process P {}` with
      `include { P } from '../modules/local/p'` and adapt the call site to pass `meta`.
- [ ] **Gate (must pass before commit):**
      ```
      nextflow config .  -profile test
      nextflow run    .  -profile test -stub-run
      nextflow run    .  --help
      ```
- [ ] Commit. One process. Repeat.

---

## Migration order (corrected)

The old plan started with RNA-seq fetch ("least interdependent") — but
`SRA_FETCH` is the **single most complex** process (~270 lines). Prove the
pattern on a leaf first, then attack the hard pieces.

| Phase | Work | Why here |
|------|------|----------|
| 0 | `meta` map contract + `INPUT_CHECK` subworkflow | Blocks all extraction; dedupes earlgrey |
| 1 | Skeleton: `main.nf`, `workflows/`, `subworkflows/local/`, `modules/local/`; move existing `asm_stats` + `annotation_tools` to convention | Establishes layout cheaply |
| 2 | `versions.yml` + `conf/base.config` + `conf/modules.config` | Pattern every later module copies |
| 3 | Setup modules (3 leaf processes) | Easiest real extraction; validates gate |
| 4 | Genome clean + `prepare_genome` subworkflow | Preserve FCS-GX `/dev/shm` staging |
| 5 | `mask` subworkflow + per-tool masker modules | Replaces the mega-process design |
| 6 | `rnaseq` subworkflow (the hard one) | Done *after* pattern is proven |
| 7 | `predict` subworkflow (train/predict/update) | Core |
| 8 | `annotate` subworkflow | Composition of optional tools |
| 9 | nf-core hygiene (docs/usage, docs/output, schema_input, naming, MultiQC) | Stretch |

---

## Provisioning repivot (done)

The `module` provisioning profile is renamed to **`ucr_hpcc`** — an
*institutional* profile in the nf-core/configs sense. The Lmod module names and
`/bigdata` paths only exist at UCR, so the name now says so. Portable runs use
`singularity` (containers) or `pixi`. New sites copy
`conf/provision_ucr_hpcc.config` to `conf/provision_<site>.config` and register a
matching profile.

```
-profile annotate,slurm,ucr_hpcc      # institutional (default on UCR HPCC)
-profile annotate,local,singularity   # portable
```

Issue 11 covers the fuller consolidation (folding UCR SLURM partitions /
`clusterOptions` into the same institutional profile).

---

## Distance from nf-core

| Area | State |
|------|-------|
| Scaffolding (schema, CITATIONS, COC, CHANGELOG, LICENSE, CI) | ~30–40% there |
| `meta` map | none (Phase 0) |
| Structure (`main.nf`/`workflows`/`subworkflows`/`modules`) | monolith |
| `versions.yml` per module + MultiQC | none (mandatory for nf-core) |
| Containers per module | **biggest gap** — relies on Lmod/pixi; nf-core needs conda+biocontainer per process |
| Naming | `nf_funannotate1` violates nf-core naming (underscores/digits) |
| nf-test, `docs/usage.md`+`output.md`, `assets/schema_input.json`, `.nf-core.yml` | missing |

**Verdict:** ~30–40% on peripheral scaffolding, ~0% on the two load-bearing items
(meta-maps + container-per-module). Full `nf-core lint` compliance is a multi-week
rewrite, much of which fights our HPC reality. **Recommendation: nf-core-inspired,
not nf-core-submitted** — adopt meta-maps, one-tool-per-module, subworkflows,
`versions.yml`, `conf/modules.config`; keep `ucr_hpcc` as an institutional profile
but add a real container path so the pipeline is portable.
