# ChIP-seq dynamic filtering runtime validation 001

This snapshot records Phase 5B of the incremental ChIP-seq
generalization.

## Objective

Phase 5B generalizes the duplicate-marking and BAM-filtering runtime
from the original two-run SINGLE-end pilot to the current dynamic
22-run ChIP-seq processing cohort.

## Cohort

The validated synthetic runtime contract contains:

- 22 unique runs;
- 18 IP runs;
- 4 shared Input runs;
- SINGLE-end libraries;
- ILLUMINA platform.

The accessions correspond to the current dynamic ChIP-seq processing
cohort. Synthetic BAM placeholders are used only to construct the
Snakemake dry-run DAG.

## Scientific method

The previously real-data-validated filtering method is preserved:

- Picard MarkDuplicates;
- duplicates marked, not removed during MarkDuplicates;
- duplicate scoring strategy SUM_OF_BASE_QUALITIES;
- optical duplicate parsing disabled;
- nuclear alignments retained;
- minimum MAPQ 30;
- MAPQ 255 excluded;
- exclusion flags 3844;
- mitochondrial accession NC_007936.1 excluded;
- final filtered BAM validation remains part of the workflow.

Nine pre-existing scientific regression tests passed.

## Runtime generalization

The filtering Snakefile no longer requires exactly two pilot runs and
no longer imports the pilot experimental-eligibility workflow.

The input plan is checked fail-closed for schema, layout, platform,
unique accessions, roles, source BAM paths, SHA-256 values, byte sizes
and alignment-count consistency.

Four new dynamic-runtime unit tests passed.

The 22-run Snakemake dry-run contains 112 planned rule jobs:

- 22 input-staging jobs;
- 22 MarkDuplicates jobs;
- 22 filtering jobs;
- 22 filtered-BAM indexing jobs;
- 22 filtered-BAM QC jobs;
- 1 cohort summary job;
- 1 final target.

PAIRED layout and duplicate-run fixtures are rejected.

## Scope boundary

Phase 5B modifies only the filtering runtime core, tests and
generalization documentation.

The filtering configuration, SLURM wrapper and submission layer remain
unchanged and are intentionally still operationally pilot-specific.
Their generalization is reserved for Phase 5C.

No productive 22-run alignment currently exists. No real duplicate
marking or BAM filtering was executed during this validation, and no
filtered BAM was created.
