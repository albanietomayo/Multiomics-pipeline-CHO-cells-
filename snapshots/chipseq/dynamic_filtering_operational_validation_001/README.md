# ChIP-seq dynamic filtering operational validation 001

This snapshot records Phase 5C of the incremental ChIP-seq
generalization.

## Objective

Phase 5C connects the generalized duplicate-marking and filtering
runtime to the dynamic ChIP-seq alignment submission layer.

## Alignment acceptance contract

The filtering submission does not trust the alignment
`latest_submission.txt` pointer alone.

It requires a productive alignment submission containing:

- numeric `job_id.txt`;
- the corresponding `job_<id>` directory;
- `job_status.tsv` with `stage=completed` and `exit_status=0`;
- a successful output-publication validation log;
- `output.sha256`;
- validated reference, input-provenance and alignment-QC outputs.

The small provenance and QC files are SHA-256 checked on the login
node. Full raw-BAM SHA-256 verification is deliberately deferred to
compute-node staging.

## Dynamic cohort

The integrated synthetic validation represents the current cohort:

- 22 unique runs;
- 18 IP runs;
- 4 shared Input runs;
- SINGLE-end;
- ILLUMINA.

The filtering plan is generated automatically from the verified
alignment submission.

## Filtering method

The scientific filtering policy remains unchanged:

- Picard MarkDuplicates;
- duplicate marking without removal at the Picard stage;
- MAPQ >= 30;
- MAPQ 255 excluded;
- exclusion flags 3844;
- mitochondrial accession NC_007936.1 excluded;
- optical duplicate parsing disabled;
- persistent verified filtered BAM and CSI outputs.

## Tests

The validation passed:

- 9 scientific filtering regression tests;
- 4 Phase-5B dynamic-runtime tests;
- 6 Phase-5C operational tests.

The integrated Snakemake dry-run expands to exactly 112 jobs:

- 22 input staging;
- 22 MarkDuplicates;
- 22 BAM filtering;
- 22 filtered-BAM indexing;
- 22 filtered-BAM QC;
- 1 cohort summary;
- 1 final target.

## Safety and storage

`--check` did not call `sbatch`, create a filtering job ID or update the
filtering latest-submission pointer.

No real filtering was executed.

No upstream alignment BAM was deleted. All 22 synthetic raw BAMs
remained present after validation.

Automatic raw-BAM deletion is intentionally not part of Phase 5C.
A separate cleanup action may only be considered after a productive
filtering job has completed, published and revalidated all persistent
filtered outputs.

## Current production boundary

The real project does not yet contain a productive generalized
22-run alignment submission. It therefore correctly fails closed and
cannot start filtering at this stage.
