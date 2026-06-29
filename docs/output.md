# Output

All paths are relative to `launchDir` (the directory from which `nextflow run` is invoked)
unless noted as params that can be overridden.

## Directory layout

```
launchDir/
├── input_clean_genomes/          # Cleaned + masked genome FASTAs
│   ├── <ASMID>.fa.gz             # FCS-GX cleaned, length-filtered
│   └── <ASMID>.masked.fasta.gz   # Soft-masked (tantan or EarlGrey path)
│
├── genome_annotation/            # Main funannotate output (params.target)
│   └── <ASMID>/
│       ├── predict_results/      # Funannotate predict outputs
│       │   ├── <ASMID>.gbk       # GenBank annotation
│       │   ├── <ASMID>.gff3      # Gene models
│       │   └── <ASMID>.proteins.fa
│       ├── annotate_results/     # Funannotate annotate outputs (run_annotate=true)
│       │   ├── <ASMID>.gbk
│       │   └── <ASMID>.gff3
│       ├── antismash/            # antiSMASH results (run_antismash=true)
│       ├── interproscan/         # InterProScan results (run_interpro=true)
│       └── signalp/              # SignalP results (run_signalp=true)
│
├── genome_annotation_training/   # Funannotate train artifacts (params.training_target)
│   └── <ASMID>/
│       └── training/
│
├── rnaseq_data/                  # Trinity assemblies + normalized reads
│   ├── <species_tag>.trinity-GG.fasta
│   ├── <species_tag>_norm_R1.fastq.gz
│   ├── <species_tag>_norm_R2.fastq.gz
│   └── <species_tag>_norm_SE.fastq.gz
│
├── rnaseq_reads/                 # Fetched + trimmed SRA reads
│   ├── <accession>_R1.fastq.gz
│   └── sra_query/                # Per-species SRA run lists
│       └── <species_tag>_sra.csv
│
├── tables/                       # Assembly statistics (params.tables_dir)
│   └── asm_stats.tsv.gz          # N50, contig count, bp — used by EarlGrey SELECT_REPS
│
├── results/repeatlibrary/        # EarlGrey TE libraries (params.earlgrey_outdir)
│   └── <species_safe>/
│       ├── <rep_ASMID>.families.fa        # Curated TE family library
│       ├── <rep_ASMID>.masked.fasta.gz    # Representative masked genome
│       ├── <species_safe>_RepeatLandscape/
│       └── <species_safe>_summaryFiles/
│
├── funannotate_db/               # Funannotate annotation databases (params.funannotate_db)
├── lib/augustus/3.5/config/      # Writable Augustus config copy (params.augustus_config)
│
└── logs/nextflow/
    ├── annotate_trace.txt        # Per-task resource usage
    ├── annotate_report.html      # Visual execution report
    ├── annotate_timeline.html    # Timeline view
    └── software_versions.yml     # Tool versions collected from all processes
```

## Key outputs

### Genome annotations (`genome_annotation/<ASMID>/`)

The central output of funannotate predict. Each assembly gets its own subdirectory
named by `ASMID`. The `predict_results/` subdirectory always exists after a
successful predict run; `annotate_results/` is added when `--run_annotate true`.

| File | Description |
|------|-------------|
| `predict_results/<ASMID>.gbk` | GenBank flat-file with gene models |
| `predict_results/<ASMID>.gff3` | GFF3 gene models |
| `predict_results/<ASMID>.proteins.fa` | Predicted protein sequences |
| `annotate_results/<ASMID>.gbk` | GenBank with functional annotations added |

### Cleaned and masked genomes (`input_clean_genomes/`)

Intermediate cleaned genomes are written here by storeDir caching. The pipeline
reads back from this directory on resumed runs — the files act as checkpoints.

| File | Description |
|------|-------------|
| `<ASMID>.fa.gz` | FCS-GX cleaned (or length-filtered only if `--skip_fcs true`) |
| `<ASMID>.masked.fasta.gz` | Soft-masked with tantan (default) or EarlGrey + RepeatMasker |

### RNA-seq evidence (`rnaseq_data/`, `rnaseq_reads/`)

RNA-seq data is organized by species tag (`{species}_{strain}` with whitespace
replaced by `_`). Normalized reads are storeDir-cached so the expensive Trinity
assembly is not repeated on `-resume`.

### Assembly statistics (`tables/asm_stats.tsv.gz`)

Generated when `--gen_asm_stats true` (default). Contains N50, contig count, and
assembled bp per assembly. Used by the EarlGrey path to select the most contiguous
representative per species for TE library construction.

### EarlGrey TE libraries (`results/repeatlibrary/`)

Only generated when `--run_earlgrey true`. The curated TE family FASTA
(`<ASMID>.families.fa`) is built once per species and reused by RepeatMasker for
all conspecific strains.

### Logs (`logs/nextflow/`)

Nextflow trace, HTML report, and timeline are written here. `software_versions.yml`
collects the tool versions emitted by each process that outputs a `versions.yml`.

## Persistence model

This pipeline uses `storeDir` (not `publishDir`) for most outputs. This means:

- **First run**: task executes and writes outputs into the storeDir path.
- **Subsequent runs**: if the output file exists and is non-empty in the storeDir path, the task is skipped entirely (no work-dir entry created).
- **Effect**: outputs accumulate across runs in a shared directory tree. Deleting a storeDir file forces that task to re-run.

`publishDir` (copy mode) is used only for: SRA query results, antiSMASH, InterProScan, and SignalP — these are copied from the Nextflow work directory into the output tree after the task completes.
