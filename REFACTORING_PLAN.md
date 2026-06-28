# Nextflow Pipeline Modularization Plan

## Pipeline Workflow Structure

```
┌─────────────────────────────────────────────────────────────────┐
│                    GENOME PREPROCESSING                          │
├─────────────────────────────────────────────────────────────────┤
│  CLEAN              │ MASK                  │ SUMMARY_STATS     │
│  ─────              │ ────                  │ ─────────────     │
│  • FCS_GX           │ • NONE                │ • ASM_STATS       │
│  • sourpurge        │ • TANTAN              │                   │
│  • vecscreen        │ • REPEATMODELER       │                   │
│                     │ • REPEATMASKER        │                   │
│                     │ • EARLGREY            │                   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                    GENE PREDICTION                               │
├──────────────────────────────────┬──────────────────────────────┤
│   RNA-seq Preparation             │  Funannotate Pipeline       │
│   ─────────────────────           │  ──────────────────        │
│   • SRA_QUERY                     │  • TRAIN                   │
│   • SRA_FETCH (PE & SE)           │  • PREDICT                 │
│   • RNASEQ_PREPARE (Trinity)      │  • UPDATE (optional)       │
└──────────────────────────────────┴──────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                    ANNOTATION                                    │
├─────────────────────────────────────────────────────────────────┤
│  • ANTISMASH (secondary metabolites)                            │
│  • SIGNALP (signal peptides)                                    │
│  • INTERPROSCAN (protein domains - when implemented)            │
│  • FUNANNOTATE_ANNOTATE (final annotation)                      │
└─────────────────────────────────────────────────────────────────┘
```

## Completed ✅
✅ **ASM_STATS Module** (`modules/asm_stats.nf`)
- Extracted assembly statistics generation into a separate reusable module
- Used by both `funannotate.nf` and `earlgrey_mask.nf`
- Generates `asm_stats.tsv.gz` with ASMID, total_length_bp, N50_bp, contig_count

✅ **Optional SELECT_REPS** (`earlgrey_mask.nf`)
- Added `--skip_select_reps` flag to skip representative selection
- When enabled, all genomes are processed for EarlGrey masking without size filtering
- Added `--gen_asm_stats` flag (default: true) to auto-generate assembly statistics

✅ **Annotation Tools Module** (`modules/annotation_tools.nf`)
- ANTISMASH_RUN (secondary metabolite detection)
- INTERPROSCAN_RUN (protein domain annotation)
- SIGNALP_RUN (signal peptide prediction)

## Proposed Modularization Structure (User-Approved)

### Phase 1: Genome Preprocessing Modules

#### `modules/AAFTF/asm_stats.nf` (Refactor existing)
**Summary Statistics Generation**
- ASM_STATS: Compute assembly stats (total_length_bp, N50_bp, contig_count)
- *Note: Move existing `modules/asm_stats.nf` here*

#### `modules/AAFTF/FCS_GX.nf` (Extract from GENOME_CLEAN)
**Contamination Screening & Removal**
- FCS_GX contamination detection and removal
- Phylum-aware filtering using NCBI taxonomy

#### `modules/AAFTF/sourpurge.nf` (Extract from GENOME_CLEAN)
**Sourpurge Contamination Detection**
- Source organism contamination screening

#### `modules/AAFTF/vecscreen.nf` (New)
**Vector Contamination Screening**
- NCBI VecScreen vector contamination detection

#### `modules/repeatmasking/masking.nf`
**Repeat Masking Strategy Selection**
- TANTAN: Soft masking (tantan algorithm)
- REPEATMODELER: De novo TE discovery
- REPEATMASKER: Library-based masking (existing library or species)
- EARLGREY: De novo TE discovery + masking (currently in separate earlgrey_mask.nf)
- NONE: Skip repeat masking entirely

### Phase 2: Gene Prediction Modules

#### `modules/rnaseq_fetch/sra_query.nf`
**SRA Discovery & Query**
- SRA_QUERY: Query NCBI SRA for RNA-seq accessions per species
- SRA_QUERY_BATCH: Batched SRA queries to NCBI
- COLLECT_SRA_QUERY: Merge per-species results into manifest

#### `modules/rnaseq_fetch/sra_fetch.nf`
**RNA-seq Download & Normalization**
- SRA_FETCH: Download paired-end RNA-seq, normalize reads
- SRA_FETCH_SE: Download single-end RNA-seq, normalize
- WRITE_EMPTY_READS: Create placeholders for species with no SRA data

#### `modules/rnaseq_fetch/prepare.nf`
**RNA-seq Assembly & Preparation**
- RNASEQ_PREPARE: Trinity assembly and normalization per species
- Output shared Trinity-GG for all strains of a species

#### `modules/funannotate/train.nf`
**Gene Model Training**
- FUNANNOTATE_TRAIN: PASA-based training on representative assembly
- Full training (Trinity + HISAT2 + trimmomatic) for representatives
- PASA-only for non-representative strains

#### `modules/funannotate/predict.nf`
**Gene Prediction**
- FUNANNOTATE_PREDICT: Ab initio and evidence-based gene prediction
- Pre-flight validation (assembly size/fragmentation checks)
- Post-prediction filtering and formatting

#### `modules/funannotate/update.nf` (Optional)
**Prediction Update with RNA-seq**
- FUNANNOTATE_UPDATE: Update predictions with mapped RNA-seq reads
- Optional step for models with available transcriptomics data

### Phase 3: Annotation Modules

#### `modules/annotate/annotation_tools.nf` (Refactor existing)
**Post-prediction Annotation Tools**
- ANTISMASH_RUN: Secondary metabolite cluster detection
- SIGNALP_RUN: Signal peptide prediction
- INTERPROSCAN_RUN: Protein domain annotation (when implemented)

#### `modules/annotate/funannotate.nf`
**Final Funannotate Annotation**
- FUNANNOTATE_ANNOTATE: Functional annotation merging

### Phase 4: Setup & Utilities

#### `modules/setup/databases.nf`
**Database Initialization**
- SETUP_TAXONDB: NCBI taxonomy database (for FCS-GX)
- SETUP_FUNANNOTATE_DB: Funannotate databases (BUSCO, etc.)
- SETUP_AUGUSTUS_CONFIG: Writable Augustus configuration

## Benefits of Modularization

1. **Reusability**: Modules can be composed into different pipelines
2. **Flexibility**: Easy to swap masking strategies or annotation tools
3. **Maintainability**: Smaller, focused files are easier to understand and modify
4. **Testing**: Individual modules can be tested in isolation
5. **Documentation**: Each module documents its inputs, outputs, and dependencies
6. **Scalability**: Easier to add new tools (e.g., new masking strategies)
7. **Git History**: Smaller commits with clear intent

## Implementation Strategy

### Priority Order
1. **Phase 2.1 (RNA-seq Fetch)**: Least interdependent, high reusability
2. **Phase 2.2 (Funannotate Modules)**: Core prediction pipeline
3. **Phase 1 (Genome Preprocessing)**: More complex due to conditional branching
4. **Phase 3 (Annotation)**: Depends on Phase 2 completion
5. **Phase 4 (Setup)**: Last, as foundational

### Implementation Notes
- Each module should be standalone with clear input/output contracts
- Use `include` statements in main workflow files
- Maintain backward compatibility with existing wrapper scripts
- Update nextflow.config to support module-specific params
- Create detailed header comments in each module file
- Use consistent naming conventions: `modules/{category}/{function}.nf`

### Testing & Validation
- Test each module in isolation with stub runs: `-stub-run`
- Verify module reuse works across different pipelines
- Document module interdependencies
- Create example usage in comments
