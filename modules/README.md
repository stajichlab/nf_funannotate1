# Nextflow Modules Directory

This directory contains reusable Nextflow process modules organized by functional area. Each module is designed to be independently importable and composable into different workflows.

## Directory Structure

```
modules/
├── README.md (this file)
├── asm_stats.nf
│   └── Currently at root level; move to AAFTF/ in Phase 1
├── annotation_tools.nf
│   └── Currently at root level; move to annotate/ in Phase 3
├── funannotate/ (Phase 2 - Gene Prediction)
│   ├── predict.nf
│   ├── train.nf
│   └── update.nf
├── AAFTF/ (Phase 1 - Genome Preprocessing)
│   ├── asm_stats.nf
│   ├── FCS_GX.nf
│   ├── sourpurge.nf
│   └── vecscreen.nf
├── repeatmasking/ (Phase 1 - Repeat Masking)
│   └── masking.nf (strategies: TANTAN, REPEATMODELER, REPEATMASKER, EARLGREY, NONE)
├── rnaseq_fetch/ (Phase 2 - RNA-seq Preparation)
│   ├── sra_query.nf
│   ├── sra_fetch.nf
│   └── prepare.nf
├── annotate/ (Phase 3 - Annotation)
│   ├── annotation_tools.nf
│   └── funannotate.nf
└── setup/ (Phase 4 - Utilities)
    └── databases.nf
```

## Module Usage

### Including a Module
```groovy
include { PROCESS_NAME } from './modules/category/module.nf'

// In workflow:
PROCESS_NAME(input_channel)
```

### Example: Using Multiple Modules
```groovy
include { ASM_STATS } from './modules/AAFTF/asm_stats'
include { FUNANNOTATE_PREDICT } from './modules/funannotate/predict'
include { ANTISMASH_RUN; SIGNALP_RUN } from './modules/annotate/annotation_tools'

workflow {
    ASM_STATS(samples, genome_dir)
    FUNANNOTATE_PREDICT(genome_channel)
    ANTISMASH_RUN(predict_output)
    SIGNALP_RUN(predict_output)
}
```

## Module Documentation Format

Each module file should include:
1. **Header comment** describing the module's purpose
2. **List of included processes**
3. **Parameter requirements** (e.g., `params.augustus_config`)
4. **Example usage** in comments

### Template
```groovy
/*
 * module_name — Brief description of what this module does
 *
 * Processes:
 *   - PROCESS_1: What it does
 *   - PROCESS_2: What it does
 *
 * Parameters required:
 *   - params.param_name: Description
 *
 * Example usage:
 *   include { PROCESS_1; PROCESS_2 } from './modules/category/module'
 *   PROCESS_1(input_channel)
 */
```

## Development Guidelines

### Before Creating a New Module
1. Check if similar functionality exists
2. Verify the process is truly independent from others
3. Document dependencies clearly

### Module Design Principles
1. **Single Responsibility**: Each module focuses on one functional area
2. **Reusability**: Processes should work in different contexts
3. **Clear Contracts**: Explicit input/output tuples
4. **Minimal Dependencies**: Avoid tight coupling to specific params
5. **Standalone Testing**: Each module can be tested with `-stub-run`

### Parameter Handling
- Prefer `params.param_name` over hardcoded values
- Document all expected params in module header
- Use sensible defaults when possible
- Avoid module-specific param namespacing (e.g., don't use `params.predict_*` in predict.nf)

## Rollout Plan

### Phase 1: Genome Preprocessing (Extract from funannotate.nf)
- `modules/AAFTF/asm_stats.nf` - Move from root
- `modules/AAFTF/FCS_GX.nf` - Extract GENOME_CLEAN
- `modules/AAFTF/sourpurge.nf` - Extract GENOME_CLEAN
- `modules/AAFTF/vecscreen.nf` - New
- `modules/repeatmasking/masking.nf` - Extract + MASKREPEAT_TANTAN_RUN

### Phase 2: Gene Prediction (Extract from funannotate.nf)
- `modules/rnaseq_fetch/sra_query.nf` - SRA processes
- `modules/rnaseq_fetch/sra_fetch.nf` - SRA download processes
- `modules/rnaseq_fetch/prepare.nf` - RNASEQ_PREPARE
- `modules/funannotate/train.nf` - FUNANNOTATE_TRAIN
- `modules/funannotate/predict.nf` - FUNANNOTATE_PREDICT
- `modules/funannotate/update.nf` - FUNANNOTATE_UPDATE (optional)

### Phase 3: Annotation (Extract from funannotate.nf)
- `modules/annotate/annotation_tools.nf` - Move from root + FUNANNOTATE_ANNOTATE
- `modules/annotate/funannotate.nf` - If needed

### Phase 4: Setup & Utilities
- `modules/setup/databases.nf` - SETUP_* processes

## Testing Modules

### Stub Run
```bash
nextflow run funannotate.nf \
  -c nextflow/nextflow.config \
  -profile test,local \
  -stub-run
```

### Unit Test (Single Process)
```bash
nextflow run -c nextflow.config \
  -profile test,local \
  -stub-run \
  --only-module modules/funannotate/predict.nf
```

## Performance Notes

### Module Extraction Impact
- Minimal performance change (include statements are compile-time)
- Negligible memory overhead from modularization
- DAG complexity unchanged

### Caching & Resume
- storeDir and publishDir behavior unchanged
- Resume functionality works across module boundaries
- Workflow checkpoints unaffected

## Links & References

- [REFACTORING_PLAN.md](../REFACTORING_PLAN.md) - Detailed phase breakdown
- Nextflow Module Documentation: https://www.nextflow.io/docs/latest/modules.html
