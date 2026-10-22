# UCR HPCC (example institutional site)

UCR HPCC is the site the pipeline was developed on. Nothing in the pipeline
requires it: the portable route is `slurm` + `singularity` + your own site config
([adding_a_site.md](adding_a_site.md)). This page collects what is specific to
UCR HPCC, as a worked example of a site.

## Profiles and files

| File | How it is loaded | What it is |
|---|---|---|
| `conf/provision_ucr_hpcc.config` | automatically by `-profile ucr_hpcc` | Lmod `module load` per process label, plus the per-process SLURM queue and `--export=ALL,SCRATCH=,TMPDIR=` safety net |
| `conf/site_ucr_hpcc.config` | automatically by `-profile ucr_hpcc` | UCR shared paths (databases, image cache, InterPro data). Its `apptainer.runOptions` binds some personal `/rhome/...` paths; other UCR users must edit them |
| `conf/site_ucr_hpcc_singularity.config` | `-c`, with the container route | UCR partitions and queues, the AAFTF 0.7 image and the `ASM_STATS` container |
| `conf/site_ucr_hpcc_preempt.config` | `-c`, last | opportunistic use of the preemptible partition |

When running by project name (`nextflow run stajichlab/nf_funannotate1`), a `-c`
file needs a path to a local checkout.

## Launching

```bash
# Lmod modules on SLURM, from a launch dir with samples.csv:
nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,ucr_hpcc -resume --n_test 1

# the sbatch launcher, from a local checkout. Default provisioning is containers
# (singularity) layered over ucr_hpcc's SLURM safety net, not Lmod:
sbatch /path/to/nf_funannotate1/run_annotate.sh --n_test 1
# the launcher with Lmod modules:
PROVISION=ucr_hpcc sbatch /path/to/nf_funannotate1/run_annotate.sh --n_test 1
```

## What `run_annotate.sh` does on UCR HPCC

The launcher detects UCR HPCC by the directory `/bigdata/stajichlab`. When it
exists:

- `SITE_PROFILE` defaults to `ucr_hpcc`, and `CONDA_ENVS_ROOT` defaults to
  `/bigdata/stajichlab/shared/condaenv`.
- `ucr_hpcc` stays in the profile list even when `PROVISION=conda` or
  `singularity` is asked for, for example `-profile annotate,slurm,ucr_hpcc,conda`.
  The profile carries the SLURM `SCRATCH`/`TMPDIR` safety net as well as the
  modules. Without it, tasks can exhaust the node-local `/tmp`. Only the
  last-loaded axis's `beforeScript` / `container` wins for each process.

Elsewhere neither default applies. `SITE_PROFILE= sbatch run_annotate.sh` turns
off only the profile; the `CONDA_ENVS_ROOT` default is a separate check of the same
directory. On UCR HPCC it also drops the SCRATCH/TMPDIR safety net, unless you pass
`-c conf/site_ucr_hpcc_singularity.config`.

## UCR-only defaults you may see in the docs

| Setting | UCR value |
|---|---|
| shared container cache (`sif_dir`) | `/bigdata/stajichlab/shared/lib/singularity_cache` |
| shared conda envs (`conda_envs_root`) | `/bigdata/stajichlab/shared/condaenv` |
| shared funannotate database | `/bigdata/stajichlab/shared/lib/funannotate_db` |
| NCBI assemblies (`source`) | `/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM` |
| taxon database (`taxondb`) | `/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/tmp/taxa/` |
| BUSCO lineages (`busco_lineages`) | `/srv/projects/db/BUSCO/v10` |
| eggNOG (`eggnog_db`) | `/srv/projects/db/eggNOG/LATEST` |
| InterPro 6 data (`iprscan6_datadir`; also `iprscan6_pipeline`) | `/srv/projects/db/interproscan/6.0.0/110.0` |
| antiSMASH database | `/bigdata/stajichlab/shared/lib/antismash_db` |
| old MariaDB image | `.../singularity_cache/mariadb.sif` (10.3.9; the funannotate image's MariaDB is the default) |

For your own site, set the equivalents in your site config
(`conf/site_template.config`).

## UCR HPCC and NRP together

UCR HPCC is also the SLURM half of the usual split with NRP. See
[multi_cluster.md](multi_cluster.md).
