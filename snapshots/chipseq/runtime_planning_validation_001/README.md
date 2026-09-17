# ChIP-seq runtime planning validation 001

This snapshot validates Phase 3A of the ChIP-seq incremental
generalization work.

## Runtime planning chain

The active metadata workflow now exposes the following planning chain:

1. ENA-derived ChIP-seq run inventory;
2. target and library-role annotation;
3. experimental-condition reconstruction;
4. control-candidate generation;
5. incremental experimental eligibility;
6. deterministic IP/Input analysis planning;
7. deduplicated processing-run selection.

## Current validated cohort

For the current validated catalogue, the analysis plan contains 18
`ready` IP/Input analyses.

These analyses require 22 unique sequencing runs:

- 18 IP runs;
- 4 shared Input runs.

The Input controls are reused as follows:

- ERR868150: 6 analyses;
- ERR868177: 6 analyses;
- SRR20770287: 3 analyses;
- SRR20770295: 3 analyses.

All 22 current runs are SINGLE-end Illumina libraries.

## Compatibility with the existing FASTQ layer

`chipseq_processing_runs.tsv` was passed directly to the existing
`build_validation_manifest.py` implementation.

The resulting FASTQ manifest contained exactly 22 runs and 22 physical
FASTQ files. The processing-run set and FASTQ-manifest run set were
identical.

No modification to `workflow/rules/common.smk` or
`workflow/scripts/build_validation_manifest.py` was required.

## Runtime

Validation used a temporary environment created from the repository
`environment.yml`.

The validated runtime contained:

- Python 3.14.7;
- Snakemake 9.25.2;
- pandas 3.0.5.

The modified metadata workflow parsed successfully with Snakemake and
exposed the new `chipseq_incremental_eligibility`,
`chipseq_analysis_plan` and `chipseq_processing_runs` rules.

## Scope

This snapshot validates runtime planning integration, processing-run
deduplication, compatibility with the existing FASTQ-manifest builder,
and Snakemake parsing.

It does not claim that a complete new ENA refresh has already been run
through the full workflow, and it does not yet redirect the ChIP-seq
preprocessing entry point to the new processing-run table.

That preprocessing redirection is handled separately in Phase 3B.
