# Citations

## Pipeline

- **nf_funannotate1** — this pipeline. Stajich Lab, UC Riverside.

## Workflow manager

- [Nextflow](https://doi.org/10.1038/nbt.3820)
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017.
- [nf-schema](https://github.com/nextflow-io/nf-schema) — parameter validation and schema-driven help.

## Core tools

Tools are invoked depending on the steps enabled for a given run.

- [funannotate](https://github.com/nextgenusfs/funannotate) — eukaryotic genome annotation pipeline (clean, mask, train, predict, annotate, update).
- [AAFTF](https://github.com/stajichlab/AAFTF) — assembly QC and contaminant filtering used during genome cleaning.
- [NCBI FCS-GX](https://github.com/ncbi/fcs) — foreign contaminant screening.
- [tantan](https://gitlab.com/mcfrith/tantan) — low-complexity / simple-repeat masking.
- [EarlGrey](https://github.com/TobyBaril/EarlGrey) + [RepeatMasker](https://www.repeatmasker.org/) — curated transposable-element masking (standalone pipeline).
- [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) — ab initio gene prediction / training.
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) — genome-guided RNA-seq assembly for training evidence.
- [PASA](https://github.com/PASApipeline/PASApipeline) — transcript alignment assemblies (uses MariaDB).
- [EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler) — evidence-based gene-model consensus (via funannotate).
- [BUSCO](https://busco.ezlab.org/) — lineage completeness assessment and prediction seeding.
- [fastp](https://github.com/OpenGene/fastp) — read trimming/QC.
- [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (bbnorm) — RNA-seq read normalization.
- [SRA Toolkit](https://github.com/ncbi/sra-tools) — SRA data retrieval.
- [TaxonKit](https://bioinf.shenwei.me/taxonkit/) — NCBI taxonomy queries.
- [antiSMASH](https://antismash.secondarymetabolites.org/) — secondary-metabolite biosynthetic gene cluster prediction.
- [InterProScan](https://www.ebi.ac.uk/interpro/) — protein domain/functional annotation.
- [SignalP](https://services.healthtech.dtu.dk/services/SignalP-6.0/) — signal-peptide prediction.

## Software environments

- [pixi](https://pixi.sh/), [Singularity/Apptainer](https://apptainer.org/), and [BioContainers](https://biocontainers.pro/) for reproducible tool provisioning.
