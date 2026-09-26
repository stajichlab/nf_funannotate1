# Output

All paths are relative to `launchDir` (the directory from which `nextflow run` is invoked)
unless noted as params that can be overridden.

## Names used below

| Name | Value | Example |
|------|-------|---------|
| `<ASMID>` | the `ASMID` column of `samples.csv` | `ANID_A4` |
| `<tag>` | `SPECIES` and `STRAIN` joined with `_`; spaces and `/#[]*?{}` become `_`; quotes are removed; only the first `;`-token of `STRAIN` is used (`lib/SampleUtils.groovy` `makeSampleTag`) | `Aspergillus_nidulans_FGSC-A4` |
| `<species_tag>` | `SPECIES` only, spaces changed to `_` (RNA-seq is shared by all strains of a species) | `Aspergillus_nidulans` |

Two genomes of the same species must have different `STRAIN` values, or they
get the same `<tag>` directory.

## Directory layout

```
launchDir/
├── <samples stem>.rnaseq_sra.csv # merged SRA run list (next to samples.csv; SRA fetch only)
├── rnaseq_se_candidates.csv      # SRA runs to review for single-end use (SRA fetch only)
├── rnaseq_blacklist_candidates.csv  # SRA runs to review for rnaseq_blacklist.csv (SRA fetch only)
├── input_clean_genomes/          # Cleaned + masked genome FASTAs (params.genome_dir)
│   ├── <ASMID>.fa.gz             # FCS-GX cleaned (or length-filtered only with --skip_fcs), gzipped
│   └── <ASMID>.masked.fasta.gz   # Soft-masked (tantan, or EarlGrey path via params.masked_dir)
│
├── genome_annotation/            # Main funannotate output (params.target)
│   ├── predict_skipped_too_small.tsv    # genomes skipped as too small + fragmented
│   ├── predict_blocked_awaiting_representative.tsv  # ab-initio reuse: strains waiting for their representative
│   ├── <tag>/
│   │   ├── predict_results/      # funannotate predict outputs
│   │   │   ├── <tag>.gbk         # GenBank annotation
│   │   │   ├── <tag>.gff3        # Gene models
│   │   │   ├── <tag>.tbl         # NCBI feature table
│   │   │   ├── <tag>.proteins.fa
│   │   │   ├── <tag>.cds-transcripts.fa
│   │   │   ├── <tag>.mrna-transcripts.fa(.gz)
│   │   │   ├── <tag>.scaffolds.fa
│   │   │   └── <tag>.stats.json
│   │   ├── predict_misc/         # funannotate checkpoints (resume in place)
│   │   │   └── ab_initio_parameters/<tag in lower case>.genemark.mod  # GeneMark model
│   │   ├── logfiles/
│   │   ├── antismash_local/      # antiSMASH results (run_antismash=true)
│   │   ├── annotate_misc/        # inputs for annotate (see below)
│   │   │   ├── iprscan.xml       # InterProScan (run_interpro=true)
│   │   │   ├── signalp.results.txt  # SignalP 6 (run_signalp=true)
│   │   │   └── TMRs.gff3         # DeepTMHMM (run_deeptmhmm=true)
│   │   ├── annotate_results/     # funannotate annotate outputs (run_annotate=true)
│   │   │   ├── <tag>.gbk
│   │   │   └── <tag>.gff3
│   │   └── update_results/       # funannotate update outputs (run_update=true)
│   └── <ASMID>/busco_completeness/   # BUSCO on predicted proteins (keyed by ASMID, not <tag>)
│
├── genome_annotation_training/   # funannotate train output (params.training_target)
│   └── <tag>/
│       └── training/             # PASA/Trinity training files, e.g. funannotate_train.pasa.gff3
│
├── rnaseq_reads/                 # RNA-seq reads per species (fetched from SRA or user-supplied)
│   ├── <species_tag>_norm_R1.fastq.gz
│   ├── <species_tag>_norm_R2.fastq.gz
│   ├── <species_tag>_norm_SE.fastq.gz
│   └── sra_query/                # Per-species SRA run lists
│       └── <species_tag>.sra_query.csv
│
├── rnaseq_data/                  # Shared Trinity assembly per species
│   └── <species_tag>.trinity-GG.fasta
│
├── tables/                       # Assembly statistics (params.tables_dir)
│   └── asm_stats.tsv.gz          # N50, contig count, bp — used by EarlGrey SELECT_REPS
│
├── results/
│   ├── genome_stats/BUSCO_genome/     # BUSCO on genomes (params.genome_stats_outdir)
│   │   └── <ASMID>.BUSCO_summary.<BUSCO_LINEAGE>.txt
│   ├── ani/                           # ANI comparisons (run_ani_reuse=true; params.ani_outdir)
│   └── repeatlibrary/                 # EarlGrey TE libraries (params.earlgrey_outdir)
│       └── <species_safe>/
│           ├── <rep_ASMID>.families.fa        # Curated TE family library
│           ├── <rep_ASMID>.masked.fasta.gz    # Representative masked genome
│           ├── <species_safe>_RepeatLandscape/
│           ├── <species_safe>_summaryFiles/
│           └── strains/<ASMID>.masked.fasta.gz  # Other strains masked with the library
│
├── funannotate_db/               # funannotate databases (params.funannotate_db)
├── busco_lineages/               # BUSCO lineage datasets (params.busco_lineages)
├── data/antismash_db/            # antiSMASH databases (params.antismash_db)
├── lib/augustus/3.5/config/      # Writable Augustus config copy (params.augustus_config)
├── work/                         # Nextflow work dir; also work/taxondb (params.taxondb)
│                                 #   and work/mysql_datadir (params.mysql_datadir)
│
└── logs/nextflow/
    ├── annotate_trace.<YYYYMMDD_HHMMSS>.txt  # Per-task resource usage, one file per run
    ├── annotate_report.html      # Visual execution report
    ├── annotate_timeline.html    # Timeline view
    └── software_versions.yml     # Tool versions collected from all processes
```

## Key outputs

### Genome annotations (`genome_annotation/<tag>/`)

The central output of funannotate predict. Each genome gets its own
subdirectory named by `<tag>` (`SPECIES_STRAIN`), not by `ASMID`. The files
inside use the same `<tag>` prefix. `predict_results/` exists after a
successful predict run; `annotate_results/` is added when
`--run_annotate true`, and `update_results/` when `--run_update true`.

| File | Description |
|------|-------------|
| `predict_results/<tag>.gbk` | GenBank flat-file with gene models |
| `predict_results/<tag>.gff3` | GFF3 gene models |
| `predict_results/<tag>.proteins.fa` | Predicted protein sequences |
| `predict_results/<tag>.cds-transcripts.fa` | Predicted CDS sequences |
| `annotate_results/<tag>.gbk` | GenBank with functional annotations added |
| `annotate_misc/iprscan.xml`, `signalp.results.txt`, `TMRs.gff3` | InterProScan, SignalP and DeepTMHMM results that `funannotate annotate` reads |
| `antismash_local/` | antiSMASH output; `antismash_local/<tag>.gbk` is passed to annotate |

Completion checks accept `<tag>.gbk` or `<tag>.gbk.gz`, so you can compress
finished GBK files.

### Cleaned and masked genomes (`input_clean_genomes/`)

Intermediate cleaned genomes are written here by storeDir caching. These files
use `ASMID`. The pipeline reads back from this directory on resumed runs — the
files act as checkpoints.

| File | Description |
|------|-------------|
| `<ASMID>.fa.gz` | FCS-GX cleaned (or length-filtered only if `--skip_fcs true`) |
| `<ASMID>.masked.fasta.gz` | Soft-masked with tantan (default) or EarlGrey + RepeatMasker |

### RNA-seq evidence (`rnaseq_reads/`, `rnaseq_data/`)

RNA-seq is organized by `<species_tag>` (`SPECIES` only, whitespace replaced
by `_`). All strains of one species share the same reads and Trinity
assembly. Reads are storeDir-cached in `rnaseq_reads/`, so SRA download is
not repeated on `-resume`. To use your own reads, put them there under the
same names (see the README Quick Start). `rnaseq_data/` holds the shared
Trinity-GG assembly.

### Assembly statistics (`tables/asm_stats.tsv.gz`)

Generated when `--gen_asm_stats true` (default). Contains N50, contig count, and
assembled bp per assembly. Used by the EarlGrey path to select the most contiguous
representative per species for TE library construction.

### EarlGrey TE libraries (`results/repeatlibrary/`)

Generated by `earlgrey_mask.nf`, or by `funannotate.nf` with
`--run_earlgrey true`. The curated TE family FASTA (`<rep_ASMID>.families.fa`)
is built once per species on a representative genome and reused by
RepeatMasker for all conspecific strains (`strains/`). The masked genomes are
also copied to `input_clean_genomes/<ASMID>.masked.fasta.gz`, where
funannotate reads them.

### Logs (`logs/nextflow/`)

Nextflow trace, HTML report, and timeline are written here. Each run writes a
new time-stamped trace file, so a resumed run does not overwrite the rows of
earlier runs. `software_versions.yml` collects the tool versions emitted by
each process that outputs a `versions.yml`.

## Persistence model

This pipeline uses `storeDir` (not `publishDir`) for most outputs. This means:

- **First run**: task executes and writes outputs into the storeDir path.
- **Subsequent runs**: if all declared outputs of the task exist in the storeDir path, the task is skipped entirely (no work-dir entry created). Nextflow checks that the files exist, not their size, so an empty file also counts.
- **Effect**: outputs accumulate across runs in a shared directory tree. Deleting a storeDir file forces that task to re-run.

Some processes write directly into the output tree from their own script (funannotate predict, train, annotate and update, antiSMASH, SignalP, DeepTMHMM). `publishDir` (copy mode) is used for the SRA query results, the merged `<samples stem>.rnaseq_sra.csv`, InterProScan, the EarlGrey masked-genome delivery (`input_clean_genomes/`) and the ANI results. Files in the launch directory that the pipeline reads back (for example `rnaseq_reads/`) act as caches: delete one to force the step that makes it to run again.
