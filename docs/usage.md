# Usage

## Quick start

```bash
nextflow run . -profile annotate,slurm,ucr_hpcc --samples samples.csv -resume
```

Profiles compose across three orthogonal axes — pick one from each:

| Axis | Options | Notes |
|------|---------|-------|
| pipeline | `annotate` | funannotate.nf (genome cleaning → masking → training → prediction → annotation) |
| executor | `slurm` \| `local` | `slurm` submits jobs via SBATCH; `local` runs on the head node |
| provisioning | `ucr_hpcc` \| `pixi` \| `singularity` | `ucr_hpcc` = UCR HPCC Lmod modules; see [Adding a site](adding_a_site.md) for other HPC environments |

## Samplesheet (`samples.csv`)

The samplesheet is a CSV with a header row. Required columns:

| Column | Type | Description |
|--------|------|-------------|
| `SPECIES` | string | Binomial species name, e.g. `Aspergillus fumigatus` |
| `ASMID` | string | Assembly ID, e.g. `GCA_000002655.1` or a local slug. Used as the output directory name and primary key. |
| `LOCUSTAG` | string | GenBank locus-tag prefix, e.g. `AFUA` |
| `BUSCO_LINEAGE` | string | BUSCO dataset, e.g. `fungi_odb10` or `saccharomycetes_odb12` |

Optional columns:

| Column | Default | Description |
|--------|---------|-------------|
| `STRAIN` | _(blank)_ | Strain/isolate identifier; first `;`-delimited token is used |
| `TRANSL_TABLE` | `1` | NCBI genetic code table (1 = standard, 4 = Mycoplasma, 12 = alt-yeast) |
| `NCBI_TAXONID` | _(blank)_ | NCBI Taxonomy ID for SRA queries (only needed if `--run_sra_fetch true`) |
| `GENOME` | _(blank)_ | Path to genome FASTA (`.fa`/`.fna`/`.fasta`/`.gz`). Absolute or relative to `launchDir`. If blank, resolved from `params.source/<ASMID>/<ASMID>_genomic.fna.gz`. |
| `PHYLUM`, `SUBPHYLUM`, `CLASS`, `ORDER`, `FAMILY`, `GENUS` | _(blank)_ | Taxonomic rank columns used by `--taxon RANK:VALUE` filtering |

The full column schema is in [`assets/schema_input.json`](../assets/schema_input.json).

### Example

```csv
SPECIES,STRAIN,ASMID,LOCUSTAG,BUSCO_LINEAGE,TRANSL_TABLE,NCBI_TAXONID,GENOME
Aspergillus fumigatus,Af293,GCA_000002655.1,AFUA,eurotiomycetes_odb10,1,746128,
Saccharomyces cerevisiae,S288C,GCA_000146045.2,YAL,saccharomycetes_odb12,12,4932,genomes/S288C.fa.gz
```

## Filtering samples at runtime

```bash
# Process only one assembly
nextflow run . -profile annotate,slurm,ucr_hpcc --asmid GCA_000002655.1

# Process only Ascomycota
nextflow run . -profile annotate,slurm,ucr_hpcc --taxon PHYLUM:Ascomycota

# Smoke-test with the first 2 assemblies
nextflow run . -profile annotate,slurm,ucr_hpcc --n_test 2 -stub-run
```

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samples` | `samples.csv` | Samplesheet path |
| `--target` | `genome_annotation/` | Output directory for annotations |
| `--source` | _(UCR HPCC path)_ | NCBI_ASM source directory for genomes without a `GENOME` column |
| `--skip_fcs` | `false` | Skip FCS-GX contaminant purge (no 470 GB highmem node needed) |
| `--run_earlgrey` | `false` | Use EarlGrey TE discovery instead of tantan for repeat masking |
| `--run_sra_fetch` | `true` | Fetch RNA-seq reads from SRA for funannotate training |
| `--run_annotate` | `false` | Run funannotate annotate (antismash/interpro/signalp must also be enabled) |
| `--run_antismash` | `false` | Run antiSMASH secondary metabolite prediction |
| `--run_interpro` | `false` | Run InterProScan functional annotation |
| `--run_signalp` | `false` | Run SignalP signal-peptide prediction (requires GPU partition) |

Run `nextflow run . --help` for the full schema-driven parameter list.

## Resume a run

```bash
nextflow run . -profile annotate,slurm,ucr_hpcc --samples samples.csv -resume
```

Nextflow's `cache = 'lenient'` (set in `nextflow.config`) handles NFS timestamp
jitter on `/bigdata`-mounted filesystems. Completed processes are skipped; only
changed or failed tasks re-run.

## Stub / dry run

```bash
nextflow run . -profile test -stub-run
```

The `test` profile is self-contained: it sets `beforeScript = ':'` for all
labels so no modules or containers are required. Stubs emit minimal placeholder
files so the full DAG can be traced without running real tools.

## EarlGrey masking path

To use EarlGrey (TE-discovery) instead of tantan repeat masking:

```bash
nextflow run . -profile annotate,slurm,ucr_hpcc \
    --run_repeatmasker false --run_earlgrey true \
    --earlgrey_outdir results/repeatlibrary \
    --masked_dir input_clean_genomes
```

EarlGrey runs once per species (on the most contiguous representative assembly)
to build a curated TE library; RepeatMasker then applies it to every conspecific
strain. Requires the `earlgrey` profile for resource defaults:

```bash
-profile annotate,earlgrey,slurm,ucr_hpcc
```
