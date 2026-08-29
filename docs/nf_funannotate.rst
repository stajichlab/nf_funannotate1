.. _nf_funannotate:

nf_funannotate1: a Nextflow front-end for funannotate at scale
================================================================

``nf_funannotate1`` is a Nextflow DSL2 pipeline that drives the
`funannotate <https://github.com/nextgenusfs/funannotate>`_ eukaryotic
(fungal-default, generic-capable) annotation workflow across many genomes at
once: clean → repeat-mask → RNA-seq fetch/train → ``funannotate predict`` →
optional antiSMASH / InterProScan / SignalP / DeepTMHMM → ``funannotate
annotate`` / ``update``. It adds batching, resumability, and per-species
ab-initio parameter reuse on top of the same ``funannotate`` commands
documented elsewhere in this manual (:ref:`predict`, :ref:`annotate`,
:ref:`update`) -- it does not replace them or change their behavior.

The pipeline repo lives at
`github.com/stajichlab/nf_funannotate1 <https://github.com/stajichlab/nf_funannotate1>`_
and runs directly from GitHub, no clone required:

.. code-block:: none

    nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,ucr_hpcc -resume

Nextflow caches the repo under ``~/.nextflow/assets/``; add ``-r
<branch|tag>`` to pin a revision. Outputs and the ``samples.csv`` /
``lib/`` assets are read from your *launch directory*, not the cached
checkout.

Orthogonal profiles
---------------------

Compose one option from each of three independent axes with
``-profile <pipeline>,<executor>,<provisioning>``:

.. list-table::
   :header-rows: 1

   * - Axis
     - Options
     - Notes
   * - pipeline
     - ``annotate`` · ``earlgrey`` · ``test`` / ``stub``
     - ``annotate`` is the full workflow above; ``earlgrey`` is the
       standalone curated-repeat-library entry point (below); ``test``/
       ``stub`` are self-contained dry runs (no SLURM, no real tools).
   * - executor
     - ``slurm`` · ``local``
     - ``slurm`` submits each process as an SBATCH job; ``local`` runs on
       the current node (useful for ``-stub-run`` and small dev runs).
   * - provisioning
     - ``ucr_hpcc`` · ``pixi`` · ``singularity``
     - How each process's tools reach ``$PATH``. ``ucr_hpcc`` is the
       institutional path (UCR HPCC Lmod modules); ``pixi`` and
       ``singularity`` are portable and runnable at any site. See
       :ref:`nf_funannotate_provisioning` below.

.. code-block:: none

    nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,ucr_hpcc -resume
    nextflow run stajichlab/nf_funannotate1 -profile annotate,local,singularity -resume

Process scripts themselves carry **no** ``module load`` lines -- the
provisioning profile supplies each process ``label`` with either a
``beforeScript`` (``ucr_hpcc``/``pixi``) or a container (``singularity``).
This is why the provisioning axis composes cleanly with either executor: the
same ``funannotate.nf`` graph runs unmodified regardless of how a given site
makes ``funannotate``/``augustus``/``gmes_petap.pl``/etc. available.

.. _nf_funannotate_provisioning:

Provisioning: modules, pixi, or containers
---------------------------------------------

``ucr_hpcc``
    Lmod ``module load`` per process label (``conf/provision_ucr_hpcc.config``).
    Institutional -- encodes UCR HPCC's module names/paths, analogous to an
    ``nf-core/configs`` institutional profile. Use this on UCR HPCC.

``pixi``
    Project-local `pixi <https://pixi.sh>`_ environments
    (``conf/provision_pixi.config``, driven by the repo's ``pixi.toml``).
    Portable to any site with ``pixi`` installed; no containers or Lmod
    needed.

``singularity``
    Apptainer/Singularity containers, one per process label
    (``conf/provision_singularity.config``). Most portable -- the same
    images run on HPC and Kubernetes. Despite the profile name
    ``singularity`` (kept for backward compatibility with existing
    invocations), the container *engine* used is Nextflow's ``apptainer {}``
    directive, not ``singularity {}`` -- functionally identical, but it
    avoids ``apptainer`` warning on every task invocation about the legacy
    ``SINGULARITYENV_*`` variable prefix.

GeneMark-ES/ET: host install vs. container
+++++++++++++++++++++++++++++++++++++++++++

``funannotate predict``'s GeneMark-ES/ET step is split out into its own
standalone process, ``GENEMARK_RUN``, upstream of ``FUNANNOTATE_PREDICT``
(so ``FUNANNOTATE_PREDICT`` itself can run fully containerized). How
``GENEMARK_RUN`` gets ``gmes_petap.pl`` depends on the provisioning axis:

.. list-table::
   :header-rows: 1

   * - Provisioning
     - GeneMark source
   * - ``ucr_hpcc``
     - ``module load genemarkESET`` (host install)
   * - ``pixi``
     - ``params.genemark_path``, pointed at your own licensed install --
       there is no conda/pixi package for GeneMark
   * - ``singularity``
     - The public, redistributable ``docker://teambraker/braker3`` image
       (BRAKER3's bundled GeneMark-ES/ET/EP/ETP 4.72), invoked with a
       manually-built ``apptainer exec`` call for just the
       GeneMark/Augustus-aux commands -- see below.

Under ``singularity`` provisioning, ``GENEMARK_RUN`` does **not** get a
plain Nextflow process ``container =`` directive. Its script also calls
``python``/``gzip`` on the bare host shell (for the shared too-small/
fragmented pre-flight guard, ``bin/asm_preflight_stats.py``), which aren't
guaranteed present inside the BRAKER3 image -- a full-script container wrap
would break those calls. Instead the process builds its own ``apptainer
exec`` invocation (``$SING``) at runtime and prefixes only the
GeneMark-family commands (``gmes_petap.pl``, ``bam2hints``,
``join_mult_hints.pl``) with it, with explicit binds for the task work
directory, the pipeline's own ``bin/`` scripts, the resolved license-key
directory, ``$TMPDIR``, and (in ET mode) the RNA-seq training BAM's
directory -- none of which fall under Apptainer's automatic mount set when a
container is launched this way rather than via the native ``container =``
directive.

**A public image does not waive GeneMark's own license.** The
``teambraker/braker3`` image ships the redistributable GeneMark *binary*
(from gatech-genemark's CC-BY-NC-SA-4.0-relicensed release, not the older
key-gated tarball), but ``gmes_petap.pl`` still expects a user-obtained
``~/.gm_key`` regardless of which image runs it -- see
:ref:`genemark-license` and the ``containers`` doc's note on this same point
for the companion ``funannotate-live`` container. ``GENEMARK_RUN`` resolves
that key's real location (it can be a symlink landing outside ``$HOME`` on
some clusters) and binds it into the container defensively, so this holds
even where the key isn't already reachable via a plain ``$HOME`` automount.

Samplesheet (``samples.csv``)
--------------------------------

.. list-table:: Required columns
   :header-rows: 1

   * - Column
     - Description
   * - ``SPECIES``
     - Binomial species name, e.g. ``Aspergillus fumigatus``
   * - ``ASMID``
     - Assembly ID, e.g. ``GCA_000002655.1`` or a local slug. Output
       directory name and primary key.
   * - ``LOCUSTAG``
     - GenBank locus-tag prefix, e.g. ``AFUA``
   * - ``BUSCO_LINEAGE``
     - BUSCO dataset, e.g. ``fungi_odb10``

.. list-table:: Optional columns
   :header-rows: 1

   * - Column
     - Default
     - Description
   * - ``STRAIN``
     - *(blank)*
     - Strain/isolate identifier
   * - ``TRANSL_TABLE``
     - ``1``
     - NCBI genetic code table
   * - ``NCBI_TAXONID``
     - *(blank)*
     - Needed only for ``--run_sra_fetch true``
   * - ``GENOME``
     - *(blank)*
     - Local FASTA path (absolute or relative to ``launchDir``); if blank,
       resolved from ``<source>/<ASMID>/<ASMID>_genomic.fna.gz``
   * - ``PHYLUM``, ``SUBPHYLUM``, ``CLASS``, ``ORDER``, ``FAMILY``, ``GENUS``
     - *(blank)*
     - Used by ``--taxon RANK:VALUE`` filtering

.. code-block:: none

    SPECIES,STRAIN,ASMID,LOCUSTAG,BUSCO_LINEAGE,TRANSL_TABLE,NCBI_TAXONID,GENOME
    Aspergillus fumigatus,Af293,GCA_000002655.1,AFUA,eurotiomycetes_odb10,1,746128,
    Saccharomyces cerevisiae,S288C,GCA_000146045.2,YAL,saccharomycetes_odb12,12,4932,genomes/S288C.fa.gz

Filtering and dry runs
-------------------------

.. code-block:: none

    # process only one assembly
    nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,ucr_hpcc --asmid GCA_000002655.1

    # process only Ascomycota
    nextflow run stajichlab/nf_funannotate1 -profile annotate,slurm,ucr_hpcc --taxon PHYLUM:Ascomycota

    # graph/dry test -- no SLURM, no tools needed
    nextflow run stajichlab/nf_funannotate1 -profile test -stub-run

Key parameters
----------------

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``--samples``
     - ``samples.csv``
     - Samplesheet path
   * - ``--target``
     - ``genome_annotation/``
     - Output directory for annotations
   * - ``--skip_fcs``
     - ``false``
     - Skip FCS-GX contaminant purge (avoids the highmem staging cost)
   * - ``--run_earlgrey``
     - ``false``
     - Use EarlGrey TE discovery instead of tantan for repeat masking
   * - ``--run_sra_fetch``
     - ``true``
     - Fetch RNA-seq reads from SRA for funannotate training
   * - ``--run_annotate`` / ``--run_antismash`` / ``--run_interpro`` / ``--run_signalp``
     - ``false``
     - Optional post-predict functional annotation steps
   * - ``--genemark_mode``
     - ``auto``
     - ``ES`` (self-training), ``ET`` (RNA-seq-hint-seeded when a training
       BAM is available), or ``auto`` (ET when possible, else ES)
   * - ``--genemark_path``
     - *(null)*
     - Licensed GeneMark install dir; used under ``pixi`` provisioning only
       (``ucr_hpcc`` uses the host module, ``singularity`` uses the
       container -- see :ref:`nf_funannotate_provisioning`)
   * - ``--gene_prediction_shared_abinitio``
     - *(disabled)*
     - Shared per-species ab-initio parameter store; enables
       ``--predict_with``-based fast reuse across conspecific strains

Run ``nextflow run stajichlab/nf_funannotate1 --help`` for the full
schema-driven parameter list.

Throughput, resumability & storage
--------------------------------------

Built to run over thousands of genomes and survive walltime kills /
orchestrator restarts:

* **Batched genome cleaning.** FCS-GX's ~470 GB database is staged into
  ``/dev/shm`` once per batch (``--clean_batch_size``, default 1000 genomes)
  rather than once per genome; already-cleaned genomes are skipped on
  resume.
* **Resumable, persistent prediction.** ``FUNANNOTATE_PREDICT`` writes
  directly into the durable per-genome output directory (no ``publishDir``
  copy), so a job killed by OOM/timeout resumes ``funannotate``'s own
  ``predict_misc/`` checkpoints in place.
* **Too-small/fragmented pre-flight guard.** Assemblies that are both too
  small and too fragmented to reach ``funannotate``'s 30-model training
  minimum are detected up front and skipped cleanly (flagged in
  ``predict_skipped_too_small.tsv``) instead of burning hours before
  aborting. Tunable via ``--predict_min_asm_bp`` / ``--predict_frag_max_n50``
  / ``--predict_frag_max_contigs``.
* **Compressed storage.** Clean/masked genomes are stored gzip-compressed;
  completion checks accept either ``.gbk`` or ``.gbk.gz``, so finished
  output can be archived without breaking resume logic.

EarlGrey repeat masking
--------------------------

``earlgrey_mask.nf`` is a separate entry point that builds a curated TE
library once per species (on its most contiguous representative assembly),
then applies it to every conspecific strain with RepeatMasker.
``funannotate predict`` consumes that mask in place of the default tantan
mask wherever it exists:

.. code-block:: none

    nextflow run earlgrey_mask.nf -c nextflow.config -profile earlgrey -resume

See also
----------

* :ref:`predict`, :ref:`annotate`, :ref:`update` -- the underlying
  ``funannotate`` commands this pipeline orchestrates
* :ref:`containers` -- the companion ``funannotate-live`` container this
  pipeline's ``singularity`` provisioning uses for the ``funannotate``
  label, and its GeneMark section (:ref:`genemark-license`)
* ``docs/adding_a_site.md`` in the ``nf_funannotate1`` repo -- porting the
  ``ucr_hpcc`` institutional profile to a different HPC site
