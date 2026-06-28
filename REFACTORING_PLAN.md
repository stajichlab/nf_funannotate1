# Nextflow Pipeline Modularization Plan

## Completed
✅ **ASM_STATS Module** (`modules/asm_stats.nf`)
- Extracted assembly statistics generation into a separate reusable module
- Used by both `funannotate.nf` and `earlgrey_mask.nf`
- Generates `asm_stats.tsv.gz` with ASMID, total_length_bp, N50_bp, contig_count

✅ **Optional SELECT_REPS** (`earlgrey_mask.nf`)
- Added `--skip_select_reps` flag to skip representative selection
- When enabled, all genomes are processed for EarlGrey masking without size filtering
- Added `--gen_asm_stats` flag (default: true) to auto-generate assembly statistics

## Proposed Modularization Structure

### Phase 1: Annotation Tools (Proposed)
**File:** `modules/annotation_tools.nf`
- ANTISMASH_RUN
- INTERPROSCAN_RUN
- SIGNALP_RUN

### Phase 2: Funannotate Core (Proposed)
**Directory:** `modules/funannotate/`

#### predict.nf
- FUNANNOTATE_PREDICT

#### annotate.nf
- FUNANNOTATE_ANNOTATE

#### update.nf
- FUNANNOTATE_UPDATE

#### train.nf
- FUNANNOTATE_TRAIN

### Phase 3: Genome Preparation (Proposed)
**Directory:** `modules/genome_prep/`

#### clean.nf
- GENOME_CLEAN
- GENOME_CLEAN_BATCH
- MASKREPEAT_TANTAN_RUN

#### rnaseq.nf
- SRA_QUERY
- SRA_QUERY_BATCH
- COLLECT_SRA_QUERY
- WRITE_EMPTY_READS
- SRA_FETCH
- SRA_FETCH_SE
- RNASEQ_PREPARE

#### setup.nf
- SETUP_TAXONDB
- SETUP_FUNANNOTATE_DB
- SETUP_AUGUSTUS_CONFIG

## Benefits of Modularization

1. **Reusability**: Modules can be used independently or in different pipelines
2. **Maintainability**: Smaller, focused files are easier to understand and modify
3. **Testing**: Individual modules can be tested in isolation
4. **Documentation**: Each module documents its inputs, outputs, and dependencies
5. **Git History**: Smaller commits with clear intent

## Implementation Notes

- Each module should be standalone and include all required processes
- Use `include` statements in the main workflow files
- Maintain backward compatibility with existing scripts and workflows
- Update documentation as modules are created

## Recommended Rollout

1. Start with Phase 1 (annotation_tools) - simplest, least interdependent
2. Move to Phase 2 (funannotate core) - incrementally extract predict, then annotate, then update
3. Complete Phase 3 (genome_prep) - most complex due to many interdependencies
