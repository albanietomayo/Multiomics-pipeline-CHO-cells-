# Incremental ChIP-seq generalization

## Purpose

The validated ChIP-seq pilot demonstrated preprocessing, alignment,
duplicate handling, filtering, MACS3 peak calling and signal QC for one
H3K4me3/Input pair.

The incremental-generalization layer separates biological analysis
selection from that pilot-specific execution configuration.

## Analysis planning policy

`workflow/scripts/build_chipseq_analysis_plan.py` consumes the existing:

- ChIP-seq run inventory;
- target and library-role annotations;
- reconstructed experimental conditions;
- Input-control candidate audit;
- curated experimental-eligibility decisions.

The planner is deliberately fail-closed.

An IP is marked `ready` only when:

1. the IP is retained by the current eligibility review;
2. its experimental condition is resolved;
3. exactly one Input candidate exists;
4. the Input has no documented issue;
5. the Input is itself retained by the eligibility review;
6. IP and Input have an exact match across the configured condition key;
7. the library layout is currently supported automatically.

Other states are represented explicitly as:

- `review_required`;
- `blocked_no_control`;
- `excluded`.

The planner therefore does not silently convert ambiguous metadata into
biological pairings.

## Incremental behaviour

The analysis plan is derived from current metadata rather than from fixed
SRR accessions.

This phase generalizes the analysis-planning layer conditional on the
upstream metadata and experimental-eligibility layers already covering
the active ChIP-seq cohort. It does not yet generalize the generation of
experimental-eligibility decisions for previously unseen ENA runs.

Accordingly, a newly discovered ChIP-seq IP/Input pair that has already
been classified by the upstream eligibility layer can enter the analysis
plan deterministically without adding its run accession to source code.

Full ENA-to-analysis incremental execution additionally requires an
incremental eligibility policy. That layer must preserve existing curated
decisions, automatically classify only cases supported by explicit rules
and evidence, and route unsupported or ambiguous new records to
`review_required`.

Ambiguous or incomplete records remain visible but are not executed
automatically.

## Pilot-specific MACS3 parameters

The source-informed `--nomodel --extsize 147` configuration validated for
PRJNA865478 is not treated as a universal ChIP-seq default.

Peak-calling parameter generalization is handled separately after
analysis-pair generalization, allowing assay/study-specific policies to
be introduced without altering the biological pairing logic.

## Incremental experimental eligibility

The second generalization layer replaces the assumption that the active
ChIP-seq catalogue must contain exactly the runs present in the original
curated eligibility review.

Previously reviewed runs retain their curated decision unchanged.

For a previously unseen run, automatic protocol eligibility is deliberately
restricted to the experimentally reviewed envelope. A new run is retained
automatically only when its study, library role, histone target and resolved
condition exactly match a combination already represented by a retained
baseline run from a study explicitly enabled for incremental classification.

Therefore, the incremental eligibility layer supports new sequencing runs
or later ENA depositions belonging to an already reviewed experimental
design without extending the biological interpretation beyond the available
evidence.

A new study, new experimental condition, new histone target, unresolved
label or documented issue is assigned `review_required`.

This distinction is deliberate: automatic catalogue refresh does not imply
automatic scientific approval of previously unseen experimental designs.

## Eligibility does not imply control pairing

Incremental protocol eligibility and IP/Input pairing are intentionally
separate decisions.

A newly discovered run may fall inside the reviewed protocol envelope and
therefore be eligible for processing while still failing automatic control
pairing. For example, if an ENA refresh introduces an additional Input
compatible with an experimental condition that already has an Input, the
control-candidate layer may convert the corresponding pairing from unique
to ambiguous.

Such cases remain `review_required`; incremental eligibility never grants
permission to select arbitrarily among multiple compatible controls.

Conversely, a new IP deposited for an already reviewed condition can reuse
an existing unique compatible Input when the control-candidate rules resolve
that association deterministically.

## Runtime planning integration

The metadata workflow now connects the previously validated incremental
eligibility and analysis-planning layers to the active ENA-derived ChIP-seq
metadata products.

The runtime chain is:

1. active ChIP-seq run inventory;
2. label-derived target and role annotation;
3. reconstructed experimental conditions and control candidates;
4. incremental protocol eligibility;
5. deterministic analysis planning;
6. deduplicated processing-run selection.

Only analyses with `analysis_status = ready` contribute runs to the
processing cohort.

IP runs and shared Input controls are deduplicated before FASTQ acquisition.
A shared Input therefore appears once in the processing-run table even when
it supports several ready analyses.

The current validated catalogue produces 18 ready IP/Input analyses and
22 unique processing runs: 18 IP libraries and four shared Input controls.

This phase does not yet redirect the preprocessing entry point to the
runtime processing-run table. The validated common FASTQ/QC/preprocessing
engine remains unchanged until the runtime planning layer is independently
validated.

## Dynamic preprocessing entry point

The ChIP-seq preprocessing entry point now consumes the dynamic
`chipseq_processing_runs.tsv` generated by the runtime analysis-planning
layer instead of the original two-run pilot selection.

The operational entry point executes two sequential Snakemake stages:

1. `chipseq_metadata_all`, which refreshes the active ChIP-seq metadata,
   incremental eligibility, analysis plan and deduplicated processing cohort;
2. `chipseq_preprocessing_all`, which uses that processing cohort as the
   input run table for the existing generic FASTQ acquisition, FastQC,
   fastp, post-processing QC and MultiQC rules.

This separation preserves modularity: the ChIP-specific planning layer
determines *which* runs are processable, while `common.smk` remains the
generic engine that determines *how* their FASTQ files are acquired and
preprocessed.

Shared Input controls are represented once in the processing cohort and
therefore are downloaded and preprocessed once even when reused by several
ready IP analyses.

The production wrapper no longer assumes exactly two pilot runs. Runtime
validation compares the complete processing-run set against the generated
FASTQ manifest and both structured QC summary tables.

At this stage the preprocessing redirection is implemented but no
full 22-run FASTQ download or preprocessing production execution has yet
been launched.

## Dynamic alignment input planning

The alignment generalization is being introduced in a separate planning
layer before changing the validated Bowtie2 workflow.

`prepare_chipseq_alignment_inputs.py` builds the alignment input cohort
from three runtime artifacts:

1. the dynamic `chipseq_processing_runs.tsv` generated by the analysis
   planning layer;
2. `preprocessing_qc_by_fastq.tsv`, which provides the project-relative
   processed FASTQ path and post-processing sequence count;
3. `preprocessing_qc_by_job.tsv`, which independently records the fastp
   post-processing read count.

The processing-run table remains authoritative for IP/Input role and
cohort membership. The two preprocessing tables must exactly cover the
same run set.

For every automatically alignable run, the planner verifies the
processed FASTQ on disk and records its SHA-256 digest, byte size,
post-fastp read count and project-relative provenance path.

The current automatic alignment envelope is deliberately fail-closed:
only SINGLE-end ILLUMINA ChIP-seq runs are admitted. PAIRED layouts,
unsupported platforms, missing or additional QC rows, inconsistent read
counts, path traversal and missing processed FASTQ files abort planning
rather than being guessed or silently ignored.

The current validated dynamic cohort contains 22 unique runs: 18 IP and
4 Input libraries. The previously validated H3K4me3 pilot
SRR20770297/SRR20770287 is a normal subset of this cohort and is not
special-cased by the planner.

Phase 4A does not yet alter the validated Bowtie2 Snakefile, alignment
parameters, SLURM wrapper or reference-building logic. Those runtime
changes are deferred until the dynamic input contract has been validated
independently.

## Dynamic alignment runtime core

After independent validation of the dynamic alignment-input planner,
the ChIP-seq alignment Snakefile was generalized from the original
two-run pilot to an arbitrary non-empty set of validated run
accessions.

The validated alignment method itself was not changed. Current
SINGLE-end runs continue to use Bowtie2 with local alignment,
the `very-sensitive-local` preset and seed 0, followed by coordinate
sorting with samtools and the existing alignment QC calculations.

The alignment input plan is validated before DAG construction. Run
accessions must be unique and syntactically valid; biological roles
must be `ip` or `input`; staged destinations must follow the canonical
`inputs/<run>.fastq.gz` contract; and, when layout and platform fields
are present, only SINGLE-end ILLUMINA inputs are accepted
automatically.

FASTQ staging was also generalized so that any number of planned runs
can be copied from a preprocessing project, verified against recorded
SHA-256 digests and byte sizes, and recorded in
`inputs/verified.json`.

The alignment Snakefile no longer imports or depends on the original
pilot eligibility workflow. The validated dynamic alignment input plan
is now the execution-level cohort contract, while the reference and
Bowtie2 index remain cohort-independent and therefore reusable across
different validated ChIP-seq cohorts.

Phase 4B deliberately leaves the production alignment submitter and
SLURM wrapper unchanged. Their replacement of pilot-specific
preprocessing discovery and eligibility checks is deferred to the
submission-layer integration phase.

The dynamic alignment runtime is validated with synthetic FASTQ
fixtures and a Snakemake dry-run. No Bowtie2 alignment is executed in
this phase.
