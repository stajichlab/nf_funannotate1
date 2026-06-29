# Adding a new HPC site

This pipeline uses an **institutional profile** model (inspired by nf-core/configs).
Site-specific concerns — tool provisioning and SLURM partition names — live in one
file, separate from the portable pipeline logic.

## What is site-specific

Two things differ between HPC sites:

| Concern | Where it lives | Example |
|---------|---------------|---------|
| Tool provisioning (`beforeScript`) | `conf/provision_<site>.config` | `module load funannotate/dev-1.8.18` |
| SLURM partition names (`queue`) | `conf/provision_<site>.config` | `queue = 'epyc'` |

Everything else — resource sizes (cpus/memory/time), pipeline params, process
logic — lives in the portable `conf/profile_annotate.config` and is unchanged
across sites.

## Steps

**1. Copy the UCR HPCC template**

```bash
cp conf/provision_ucr_hpcc.config conf/provision_<yoursite>.config
```

**2. Replace Lmod module names**

Each `withLabel` block has a `beforeScript` that loads tools. Swap the UCR module
names for whatever your site calls them:

```groovy
// UCR HPCC
withLabel: 'funannotate' {
    beforeScript = '''
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load funannotate/dev-1.8.18
        module load fastp
    '''
}

// Your site — replace module names as needed
withLabel: 'funannotate' {
    beforeScript = '''
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load Funannotate/1.8.18-foss-2023a
        module load fastp/0.23.4
    '''
}
```

If your site uses conda/pixi/singularity instead of Lmod, see
`conf/provision_pixi.config` and `conf/provision_singularity.config` for those
patterns — they follow the same `beforeScript`/`container` convention.

**3. Replace SLURM partition names**

The `withName` block at the bottom of the file maps each process to a UCR
partition (`short`, `epyc`, `highmem`, `gpu`). Replace them with your site's
partition names:

```groovy
// UCR HPCC
withName: '.*:FUNANNOTATE_TRAIN' { queue = { task.attempt <= 2 ? 'epyc' : 'highmem' } }

// Your site
withName: '.*:FUNANNOTATE_TRAIN' { queue = { task.attempt <= 2 ? 'medium' : 'large' } }
```

If your SLURM cluster doesn't use named partitions (or uses a default), you can
remove the `queue =` lines entirely — Nextflow will submit to the default partition.

**4. Register the profile in `nextflow.config`**

Add a one-liner under the `provisioning axis` comment in `profiles {}`:

```groovy
yoursite { includeConfig 'conf/provision_yoursite.config' }
```

**5. Run**

```bash
nextflow run . -profile annotate,slurm,yoursite -resume
```

## Portable alternatives

If Lmod modules aren't available at your site, use one of the existing portable
provisioning profiles:

| Profile | How tools are provided |
|---------|----------------------|
| `singularity` | Container image per process label (see `provision_singularity.config`) |
| `pixi` | Project-local pixi environments (see `provision_pixi.config`) |

Note: `funannotate` and `signalp6-gpu` lack public biocontainers; the singularity
profile points at locally-built `.sif` images (paths in `params.sif_dir`).
