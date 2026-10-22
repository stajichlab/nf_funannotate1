# Adding a new HPC site

Most sites need **no new profile**. Use the portable `slurm` + `singularity`
profiles and describe your cluster in one small config file:

```bash
nextflow run stajichlab/nf_funannotate1 \
    -profile annotate,slurm,singularity -c mysite.config -resume
```

## Option A (recommended): `slurm` + `singularity` + a site config

Requirements: Apptainer/Singularity on every compute node, a shared filesystem
for the launch directory and the image cache, internet on the first run.

**1. Copy the template**

```bash
cp ~/.nextflow/assets/stajichlab/nf_funannotate1/conf/site_template.config mysite.config
```

**2. Set what differs at your site.** Everything is optional:

| Concern | Setting | Default if omitted |
|---------|---------|--------------------|
| Container cache | `sif_dir` in your params file, `--sif_dir`, or export `NXF_APPTAINER_CACHE`. Not in this `-c` config (see below) | `<launchDir>/singularity_cache` |
| Partition | `process.queue`, or `withLabel:`/`withName:` `queue` | cluster default partition |
| Account / QOS / extra flags | `process.clusterOptions` (keep `-N 1 -n 1`) | `-N 1 -n 1` |
| Shared databases | `params.*` paths (see the template) | downloaded under `<launchDir>` |
| GPUs (SignalP / DeepTMHMM) | `--signalp_gpu` / `--deeptmhmm_gpu`, plus a GPU `queue` | CPU |

```groovy
// (set sif_dir in your params file, not here: a `-c` file is read after the image paths are built)
process {
    queue          = 'batch'
    clusterOptions = '-N 1 -n 1 --account=mylab'
    withName: '.*:FUNANNOTATE_TRAIN' { queue = { task.attempt <= 2 ? 'medium' : 'large' } }
}
```

Use a closure for `queue` only when it depends on `task.attempt`; otherwise a
plain string works.

**3. Run a stub check, then one genome**

```bash
nextflow run stajichlab/nf_funannotate1 -profile test,slurm,singularity -c mysite.config -stub-run
nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,singularity -c mysite.config --n_test 1 -resume
```

The stub run submits real (empty) SLURM jobs, so it checks your partition and
account settings without needing the tools or databases.

## Option B: an institutional profile (Lmod modules)

Only for sites without a container runtime. UCR HPCC does this with
`conf/provision_ucr_hpcc.config` (Lmod `beforeScript` per process label) and
`conf/site_ucr_hpcc.config` (shared paths); see [site_ucr_hpcc.md](site_ucr_hpcc.md). To add one:

1. `cp conf/provision_ucr_hpcc.config conf/provision_<yoursite>.config`
2. Replace the module names in each `withLabel` `beforeScript`
   (`module load funannotate/dev-1.8.18` becomes whatever your site calls it)
   and the partition names in the `withName` blocks.
3. Register it in `profiles {}` in `nextflow.config`, mirroring the `ucr_hpcc`
   entry (which also sets `process.shell = ['/bin/bash', '-l']` so `module load` works):

```groovy
yoursite { process.shell = ['/bin/bash', '-l']; includeConfig 'conf/provision_yoursite.config' }
```

4. Run: `nextflow run . -profile annotate,slurm,yoursite -resume`

## Other portable provisioning profiles

If you cannot use containers, these also work without site-specific modules:

| Profile | How tools are provided |
|---------|----------------------|
| `conda` | Shared frozen conda envs (see `environments/conda/`) |
| `pixi` | Project-local pixi environments (see `provision_pixi.config`) |

Note: `funannotate` and `signalp` lack public biocontainers; the singularity
profile points at locally-built `.sif` images (paths in `params.sif_dir`).
