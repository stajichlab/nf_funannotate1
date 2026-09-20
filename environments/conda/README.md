# Conda provisioning for the funannotate pipeline

This directory holds **versioned conda env manifests**, one per funannotate
release/variant you want to test under the `-profile conda` provisioning axis.
They are *deployment artifacts* — identical in role to the `.sif` images under
`params.sif_dir`, not the monorepo `pixi.toml`.

Why separate env files instead of extending `pixi.toml`? Master/beta install the
pip package (`pip install git+https://...`) **on top of** pinned conda deps, and
the Rust-EVM variant swaps the EVM backend entirely. A single drifting
`pixi.toml` can't express pip-installed git refs plus a Rust binary build; a
per-release yaml freezes a known-good lock per release and lets you diff what
changed between them.

## Files

| Manifest | funannotate | EVM backend | conda package? |
|---|---|---|---|
| `funannotate-1.8.17.yml` | 1.8.17 (released) | perl `evidencemodeler` | yes |
| `funannotate-master.yml` | git master | perl `evidencemodeler` | no — `pip install git+...` |
| `funannotate-1.9.0-beta.10.yml` | 1.9.0-beta.10 (tag) | perl `evidencemodeler` | no — `pip install git+...` |
| `funannotate-1.9.0-beta.10-rust.yml` | 1.9.0-beta.10 (tag) | Rust EVM + PASA + Trinity (`rust_optimize` forks) | no — `pip install git+...` + source builds (see below) |
| `funannotate-1.9.0-beta.11.yml` | 1.9.0-beta.11 | perl `evidencemodeler` | no — `pip install git+...` |
| `funannotate-1.9.0-beta.11-rust.yml` | 1.9.0-beta.11 | Rust EVM + PASA + Trinity (`rust_optimize` forks) | no — `pip install git+...` + source builds (see below) |
| `funannotate-1.9.0-beta.12.yml` | 1.9.0-beta.12 (tag) | perl `evidencemodeler` | no — `pip install git+...` |
| `funannotate-1.9.0-beta.12-rust.yml` | 1.9.0-beta.12 (tag) | Rust EVM + PASA + Trinity (`rust_optimize` forks) | no — `pip install git+...` + source builds (see below) |
| `nf_funannotate1-aux.yml` | — | shared *peripheral* env | yes |

**1.9.0-beta.12 is the current 1.9 target** (`build_1.9.sh` / `build_1.9_rust.sh`
default to it). The beta.10 and beta.11 manifests are kept for reproducing older
benchmark cells; beta.10's pins have drifted and it is not maintained. beta.12
differs from beta.11 in the pip ref only — it fixes an
`aux_scripts/augustus_parallel.py` crash that killed every parallel Augustus
worker on a predict run with zero protein and zero RNA-seq evidence.

`nf_funannotate1-aux.yml` is the second env the `conda` axis activates: one
shared aux env serving every peripheral tool label (edirect / sra /
genome_clean / skani / busco / prodigal / antismash / interproscan / repeatmask
/ earlgrey / select — see `conf/provision_conda.config`). It mirrors the
per-feature lists in `pixi.toml` plus the two EarlGrey-path labels; build it
exactly like a funannotate release env.

Only **1.8.17** is a pure conda package on bioconda. Master and every 1.9.0-beta
have **no conda release** — those env files install the conda runtime toolchain,
then `pip install` the git ref on top. Confirm with `funannotate check
--show-versions` after creating them and install any newly-required python deps
that drifted in.

## salmon version — which one each env actually runs

Trinity's genome-guided isoform filtering shells out to `salmon` from PATH, via
its bundled `util/support_scripts/salmon_runner.pl`. **salmon 2.x breaks that
call.** `salmon_runner.pl` passes Trinity's supporting-reads file
(`single.fa.SR_supp.fa`) to `salmon quant -r`, and that file is FASTA; salmon 1.x
read FASTA there, salmon 2.7.0 parses `-r` as FASTQ only and aborts, misreporting
the FASTA as a truncated FASTQ:

```
Error: <reads>.fa looks truncated: it holds 305 lines, which is not a whole
       number of 4-line FASTQ records
```

The Trinity run then produces a **0-byte `training/trinity.fasta`** ("0
transcripts derived from Trinity") after burning the full assembly runtime. This
is an input-format regression, not the dropped classic-CLI flags an earlier note
in `funannotate-1.9.0-beta.11-rust.yml` blamed; salmon 2.3.1 still works, 2.7.0
does not. Isolated repro (200-read FASTA + 1-transcript index, no Trinity
involved) and the full write-up live in that yml's `salmon` comment.

So every manifest that installs Trinity pins **`salmon =1.10.3`**, the last
classic 1.x LTS. The pin is NOT optional for the `-rust` manifests either: the
`rust_optimize` Trinity fork uses the same `salmon_runner.pl` and fails the same
way.

### Pinned in the manifest vs. installed in the built env

The pin only takes effect when the env is (re)built. Manifest pins and the
versions actually present under `/bigdata/stajichlab/shared/condaenv`, checked
2026-09-19:

| env | manifest pin | installed in built env |
|---|---|---|
| `funannotate-1.8.17` | `salmon =1.10.3` | 1.10.3 ✅ |
| `funannotate-1.9.0-beta.11` | `salmon =1.10.3` | 1.10.3 ✅ |
| `funannotate-1.9.0-beta.11-rust` | `salmon =1.10.3` | **2.7.0 ❌ — env predates the pin (built 2026-08-30), never rebuilt** |
| `funannotate-1.9.0-beta.12` | `salmon =1.10.3` | not built yet |
| `funannotate-1.9.0-beta.12-rust` | `salmon =1.10.3` | not built yet |
| `funannotate-master`, `funannotate-1.9.0-beta.10` | unpinned (`salmon`) | not built / stale — an unpinned solve resolves to 2.x today |
| `funannotate-1.9.0-beta.10-rust` | `salmon >=1.0` | not built / stale — same, `>=1.0` does not exclude 2.x |

**`funannotate-1.9.0-beta.11-rust` is therefore still broken on disk** and will
keep emitting empty Trinity assemblies until it is rebuilt. Rebuilding it, or
building the beta.12 envs, picks up the correct pin.

Check any env before trusting it:

```bash
/bigdata/stajichlab/shared/condaenv/<env>/bin/salmon --version   # want: salmon 1.10.3
```

## Build (once per release)

Env manifests list channels `bioconda, conda-forge`. Build a **frozen conda env**
into the shared, node-visible root (default `/bigdata/stajichlab/shared/condaenv`
— a per-user `~/conda/envs` is NOT readable by other nodes' SLURM jobs) with the
wrapper script:

```bash
./environments/conda/build_conda_env.sh funannotate-1.8.17    # one env
./environments/conda/build_conda_env.sh nf_funannotate1-aux   # the peripheral env
./environments/conda/build_conda_env.sh --all                  # every manifest (incl. aux)
# verify the tools resolve
/bigdata/stajichlab/shared/condaenv/funannotate-1.8.17/bin/funannotate check --show-versions
```

The wrapper is idempotent (skips already-built prefixes) and uses mamba for fast
solves; large solves (trinity/R-stack) occasionally hit solver dead-ends — bump
memory or retry with `--refresh` if a solve fails. `--refresh` rebuilds over an
existing prefix.

> **Build one env at a time.** Every build writes through the same conda package
> cache (`~/.conda/pkgs`). Two concurrent builds downloading the same package
> race on the `<pkg>.conda.partial` → `<pkg>.conda` rename and one of them dies
> with `[Errno 2] No such file or directory: '.../<pkg>.conda.partial'`. Observed
> 2026-09-19 when `build_1.9.sh` and `build_1.9_rust.sh` were submitted together
> and landed on the same node: the rust build failed after 5 minutes, the
> non-rust one (which won every race) was unaffected. Serialize them —
> `sbatch --dependency=afterany:<first job id> …` — or give each its own
> `CONDA_PKGS_DIRS`.

> **`--sbatch` from inside another SLURM allocation.** Both `build_1.9*.sh`
> self-submit only when `FUNANNOTATE_BUILD_JOB` is unset, a sentinel the sbatch
> call itself exports. Earlier they gated on `SLURM_JOB_ID`, which is set inside
> *any* allocation, so `--sbatch` from an interactive `srun`/`salloc` shell
> silently built inline on that shell's (usually much smaller) allocation
> instead of submitting. Fixed 2026-09-19.

These are heavy, multi-GB, entire-toolchain builds. Build them as **pre-built
envs** (not conda packages — the only non-conda component is the funannotate
python module itself, which is pip, not conda-build), because they:
- are immediately activatable by every compute node (no per-node install),
- are reproducible across runs once frozen (unlike re-solving each time),
- let you re-point the whole pipeline at a release by flipping one string.

For the `-master` and `-1.9.0-beta.*` manifests, the `pip:` block installs
funannotate at build time. They are not fully frozen: run `funannotate check`
and reconcile any version drift.

### Rust build (the `-rust` manifests — EVM + PASA + Trinity)

The Rust variant source-builds the Rust-optimized `rust_optimize` forks of **EVM,
PASA, and Trinity** (plus bowtie2) into the env. The env yaml already pins the
build toolchain (rust, c/cxx-compiler, cmake, …) and each tool's runtime deps,
mirroring `~/projects/funannotate/funannotate-live/pixi.toml`. **One-shot
wrapper** that does the whole end-to-end build (conda env + forks + activate.d +
GLIBC smoke test + `funannotate check`):

```bash
./environments/conda/build_1.9_rust.sh          # build env + Rust forks
./environments/conda/build_1.9_rust.sh --sbatch # submit as one SLURM job instead
```

### Non-rust build (the plain 1.9 manifests — PERL EVM)

The no-Rust 1.9 variant needs no wrapper work: every external is already a conda
package, so building it is just creating the frozen env, then verifying the perl
stack resolves and funannotate pins `FUNANNOTATE_EVM_ENGINE=perl`. One
`build_1.9.sh` does that:

```bash
./environments/conda/build_1.9.sh          # build env, verify perl EVM stack
./environments/conda/build_1.9.sh --sbatch # submit as one SLURM job instead
```

Why it's simpler than the rust build and what it checks:
- **No source builds, no hand-written activate.d** — conda `trinity`/`pasa`/
  `evidencemodeler`/`augustus`/`snap` packages set `TRINITY_HOME`/`PASAHOME`/
  `EVM_HOME`/`AUGUSTUS_*`/`QUARRY_PATH`/`ZOE` through their own
  `etc/conda/activate.d` scripts on activation.
- **Perl-engine guard**: as of 1.9, funannotate auto-detects the EVM engine when
  `FUNANNOTATE_EVM_ENGINE` is unset — a bare `evidence_modeler` binary on PATH
  silently switches it to the **Rust** engine (`predict.py:448`,
  `funannotate-runEVM.py:19`). The conda perl EVM package only installs
  `evidence_modeler.pl`, so auto-detect lands on perl, but the wrapper asserts
  no bare `evidence_modeler` exists and verifies via `FUNANNOTATE_EVM_ENGINE=perl`.
- **No GLIBC smoke test needed**: the solve on this host is
  `__glibc`-self-constraining (Rocky 8 = glibc 2.28; the solver only picks
  packages satisfying the host's virtual `__glibc`), unlike rust builds whose
  cargo/make artifacts bypass the solver — the exact gap `build_1.9_rust.sh`'s
  GLIBC phase exists to catch.

Keep the manifest's `python >=3.6,<3.9` in line with the reference
`pixi.toml` / 1.8.17 recipe (the 1.9 branch is only reference-tested on
`<3.9`) — the `-python` pin is the one thing to watch in
`funannotate-1.9.0-beta.12.yml`.

It uses the *same* helper scripts the reference project's other packaging paths
use (`install_scripts/pixi_install_{bowtie2,trinity,evm,pasa}.sh` +
`pixi_setup_symlinks.sh`), in the Dockerfile build order (bowtie2 → trinity →
evm → pasa), so pixi/Docker/conda builds can't drift apart. Point it at a
different checkout with `FUNANNOTATE_LIVE=/path/to/funannotate-live` if needed;
per-build logs land in `logs/conda_builds/`.

Manual equivalent (what the wrapper automates):

```bash
FUN=~/projects/funannotate/funannotate-live/install_scripts
PREFIX=~/conda/envs/funannotate-1.9.0-beta.12-rust
# copy the helper scripts into the env's bin/ so they run against its CONDA_PREFIX
# (each script clones/checks-out the rust_optimize git branch and MAKE/CARGO-BUILDs
# into $CONDA_PREFIX/opt/..., then symlinks the binaries into $CONDA_PREFIX/bin).
source /rhome/jstajich/.pixi/bin/activate 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$PREFIX"
bash $FUN/pixi_install_evm.sh       # Rust EVidenceModeler   -> $CONDA_PREFIX/bin/evidence_modeler
bash $FUN/pixi_install_trinity.sh   # Rust Trinity           -> $CONDA_PREFIX/opt/trinityrnaseq
bash $FUN/pixi_install_pasa.sh      # Rust PASA              -> $CONDA_PREFIX/opt/pasa
bash $FUN/pixi_install_bowtie2.sh   # bowtie2 source build   (AVX2 fallback fix)
bash $FUN/pixi_setup_symlinks.sh    # fasta -> fasta36
export FUNANNOTATE_EVM_ENGINE=rust
```

> These are heavy source builds and take a while — run them once, not on every
> activation. Each script is idempotent (skips if already built) and tracks the
> upstream `rust_optimize` git branch by default (pin a commit via
> `TRINITY_RUST_COMMIT` / `PASA_RUST_COMMIT` / `EVM_RUST_COMMIT` for
> reproducibility).

Then `funannotate check --show-versions` should report the Rust EVM/PASA/Trinity
rather than perl `evidence_modeler.pl` / conda `trinity` / conda `pasa`. Pairs
with the `-rust` helper scripts' env vars the way the reference pixi env does
(`TRINITYHOME`, `PASAHOME`, `EVM_HOME`, `FUNANNOTATE_EVM_ENGINE=rust`).

The wrapper also installs an `etc/conda/activate.d/funannotate-rust.sh` snippet
into the frozen env that exports exactly those rust env vars at `conda activate`
time (along with `PERL5LIB` for the PASA hooks/PerlLib and `opt/pasa/bin` on
PATH) — `AUGUSTUS_*`, `QUARRY_PATH`, and `ZOE` are NOT needed there because the
conda `augustus`/`codingquarry`/`snap` packages set them via their own
`activate.d`. It finishes with a GLIBC floor smoke test: the source-built forks
are compiled with the conda-forge c/cxx toolchain and can require
`GLIBC_2.38` (cf. `Dockerfile.base`), which a Rocky 8 host (glibc 2.28) cannot
run — the test compares the max `GLIBC_*` symbol each `opt/` binary needs
against the host glibc and fails the build loudly if jobs would die at exec.

## Run

```bash
# pick a release with --conda_env (a subdir of conda_envs_root):
nextflow run <org>/nf_funannotate1 -profile annotate,slurm,conda \
    --conda_env funannotate-1.8.17 -resume
```

The `conda` profile activates `/bigdata/stajichlab/shared/condaenv/<conda_env>`
in each funannotate label's `beforeScript` (see `conf/provision_conda.config`).
Override the root or env, e.g.:
`--conda_envs_root /path/to/envs --conda_env funannotate-1.9.0-beta.12-rust`.

## What this profile does NOT cover

`conda` only provides the pure-funannotate labels cleanly. Tooling that is
licensed or won't install from conda stays on another provider (use the
`ucr_hpcc` or `singularity` profile for those):

- **GeneMark-ES/ET** — licensed, no conda package. Under `conda` set
  `genemark_container_mode=false --genemark_path /path/to/gmes_petap_dir` (host
  install), exactly like the `pixi` profile.
- **SignalP / interproscan / antismash** — licensed or best-effort on conda;
  leave those labels on `singularity`/`ucr_hpcc`.
- **AAFTF (FCS-GX purge)** — IS conda-solvable, but only from the aaftf channel
  (`stajichlab` on anaconda.org, or the local `file://.../dist-conda` build; see
  `nf_funannotate1-aux.yml`). Make sure the channel is in the manifest before
  building the aux env; build with `MAMBA=micromamba` when it is a `file://`
  channel (mamba 1.x cannot solve local channels).

The **PASA MySQL (mariadb)** step already runs via `apptainer exec` independent
of the provisioning axis, so it behaves identically under `conda`.