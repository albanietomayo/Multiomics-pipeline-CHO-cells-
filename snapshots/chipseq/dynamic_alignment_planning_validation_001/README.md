# ChIP-seq dynamic alignment input planning validation 001

This snapshot records Phase 4A of the incremental ChIP-seq
generalization.

## Objective

Phase 4A introduces a fail-closed planning layer between the dynamic
preprocessing cohort and the validated ChIP-seq alignment workflow.

The planner does not modify Bowtie2 execution, reference generation,
alignment parameters or the SLURM alignment wrapper.

## Runtime contract

The planner consumes three artifacts:

1. `chipseq_processing_runs.tsv`, which defines the dynamic processing
   cohort and biological IP/Input role;
2. `preprocessing_qc_by_fastq.tsv`, which defines the project-relative
   processed FASTQ path and post-fastp sequence count;
3. `preprocessing_qc_by_job.tsv`, which independently records the
   post-fastp read count.

All three run sets must match exactly.

## Current cohort

The current validated processing cohort contains:

- 22 unique runs;
- 18 IP runs;
- 4 Input runs;
- 22 SINGLE-end libraries;
- 22 ILLUMINA libraries.

The previously validated pilot runs SRR20770287 and SRR20770297 are a
normal subset of this dynamic cohort and are not special-cased by the
planner.

## Fail-closed behavior

Automatic alignment planning currently supports only SINGLE-end
ILLUMINA ChIP-seq runs.

Planning aborts when it encounters unsupported paired layouts,
unsupported platforms, missing or additional preprocessing QC rows,
inconsistent post-fastp read counts, unsafe relative paths, missing
processed FASTQ files or inconsistent cohort metadata.

## Validation scope

Eight unit tests were executed successfully.

The current real 22-run processing cohort structure was also passed
through the planner using synthetic processed FASTQ files with the same
run identities, IP/Input roles, layout and platform as the validated
runtime-planning cohort.

This validates cohort propagation and planner behavior, not biological
or production alignment results.

The complete production preprocessing of these 22 runs has not yet been
executed. Therefore the real 22 processed FASTQ files were not available
for this Phase 4A validation.

No real 22-run alignment was executed and no SLURM job was submitted.

The previously validated ChIP-seq alignment runtime remains unchanged in
Phase 4A.
