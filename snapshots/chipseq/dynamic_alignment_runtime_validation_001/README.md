# ChIP-seq dynamic alignment runtime validation 001

This snapshot records Phase 4B of the incremental ChIP-seq
generalization.

## Objective

Phase 4B generalizes the previously validated two-run ChIP-seq
alignment runtime to the complete dynamic processing cohort.

The Bowtie2 alignment method itself is unchanged.

## Current runtime cohort

The validated cohort contains:

- 22 unique runs;
- 18 IP runs;
- 4 Input runs;
- SINGLE-end libraries only;
- ILLUMINA sequencing only.

## Dynamic runtime

The original two-run restriction was removed from the alignment
Snakefile.

The generic alignment worker is now `chipseq_align_run`.

The runtime input plan is the execution-level cohort contract.

The alignment Snakefile no longer imports the legacy pilot
eligibility workflow and no longer depends on `chipseq_pilot.json` or
`pilot_eligibility.json`.

The reference and Bowtie2 index remain cohort-independent resources.

## Preserved method

Phase 4B preserves the previously validated alignment method:

- Bowtie2 2.5.5;
- local alignment;
- `very-sensitive-local`;
- seed 0;
- SINGLE-end `-U` input;
- coordinate sorting with samtools;
- BAM quickcheck;
- existing alignment QC calculations.

## Validation

Phase 4A regression tests and Phase 4B runtime tests passed.

The real 22-run cohort structure was propagated through synthetic
processed FASTQ fixtures.

All 22 FASTQs were staged and verified through the generalized support
script.

A Snakemake 9.25.2 dry-run produced the expected 72-job DAG:

- 1 nuclear reference;
- 1 mapping reference;
- 1 Bowtie2 index;
- 1 FASTA index;
- 22 alignment jobs;
- 22 BAM index jobs;
- 22 alignment QC jobs;
- 1 alignment summary;
- 1 final aggregation target.

No pilot eligibility rule or pilot alignment worker was present in the
DAG.

No Bowtie2 alignment was executed, no BAM was created and no SLURM job
was submitted.

The production submitter and SLURM wrapper remain unchanged and are
therefore outside the scope of Phase 4B.
