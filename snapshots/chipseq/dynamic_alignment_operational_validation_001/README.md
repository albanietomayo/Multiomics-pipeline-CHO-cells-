# ChIP-seq dynamic alignment operational validation 001

This snapshot records Phase 4C of the incremental ChIP-seq
generalization.

## Objective

Phase 4C connects the previously validated dynamic preprocessing,
alignment-input planning and alignment runtime layers through the
production alignment submission interface.

## Upstream contract

The alignment submitter resolves the preprocessing submission referenced
by `results/chipseq/preprocessing/slurm/latest_submission.txt`.

The latest submission is not trusted merely because it was submitted.

Before alignment planning, the submitter requires:

- a numeric preprocessing `job_id.txt`;
- the corresponding `job_<id>/job_status.tsv`;
- workflow exit status 0;
- final exit status 0;
- `job_<id>/output.sha256`;
- successful revalidation of the preprocessing output manifest;
- the dynamic processing-run table;
- the per-FASTQ preprocessing QC table;
- the per-job preprocessing QC table.

If these requirements are not met, alignment planning fails closed.

## Dynamic alignment plan

The Phase 4A planner generates `config/chipseq_alignment_inputs.json`
from the verified preprocessing project.

The current validated operational cohort contains:

- 22 unique processing runs;
- 18 IP libraries;
- 4 Input libraries;
- SINGLE-end data only;
- ILLUMINA data only.

Every processed FASTQ is bound to the plan by SHA-256 digest and byte
size.

## Submission safety

The `--check` mode was executed against an isolated synthetic successful
preprocessing submission.

A fake `sbatch` sentinel verified that `sbatch` was not called.

The check created neither an alignment `job_id.txt` nor an alignment
`latest_submission.txt` pointer.

## Runtime integration

All 22 processed FASTQ fixtures were staged and SHA-256 verified.

A Snakemake 9.25.2 dry-run generated the expected 72-job alignment DAG:

- 1 nuclear reference job;
- 1 mapping reference job;
- 1 Bowtie2 index job;
- 1 FASTA index job;
- 22 alignment jobs;
- 22 BAM index jobs;
- 22 alignment QC jobs;
- 1 alignment summary job;
- 1 final aggregation target.

No pilot eligibility rule and no pilot alignment worker were present.

## Scope limitation

The synthetic preprocessing submission reproduces the structure and
current 22-run cohort of the dynamic workflow, but its processed FASTQ
files are synthetic fixtures.

The real 22-run dynamic preprocessing workflow has not yet been
executed.

No production Bowtie2 alignment was executed, no BAM was generated and
no SLURM job was submitted during Phase 4C validation.
