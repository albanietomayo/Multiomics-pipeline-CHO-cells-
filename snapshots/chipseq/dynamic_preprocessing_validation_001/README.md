# ChIP-seq dynamic preprocessing validation 001

This snapshot records Phase 3B of the incremental ChIP-seq
generalization work.

## Objective

The ChIP-seq preprocessing entry point was generalized so that it no
longer depends on the original two-run pilot selection.

Instead, preprocessing consumes the dynamic
`chipseq_processing_runs.tsv` produced by the runtime ChIP-seq
analysis-planning layer.

## Validated current cohort

The validated analysis plan currently requires:

- 18 IP runs;
- 4 shared Input runs;
- 22 unique processing runs in total;
- 22 physical FASTQ files;
- 36.434 GiB of FASTQ data;
- 22 SINGLE-end libraries.

Shared Input controls are therefore represented once in the processing
cohort and are not redundantly downloaded or preprocessed for every IP
analysis that reuses them.

## Effective preprocessing DAG

A real Snakemake dry-run using the validated 22-run cohort expanded to
93 jobs:

- 22 `download_validation_fastq`;
- 22 `fastqc_raw`;
- 22 `fastp_single`;
- 22 `fastqc_post`;
- 1 `validation_download_report`;
- 1 `multiqc_raw`;
- 1 `multiqc_post`;
- 1 `summarize_preprocessing_qc`;
- 1 `chipseq_preprocessing_all`.

No `fastp_paired` job was scheduled for the current cohort.

The dry-run created zero FASTQ files.

## Submission-wrapper validation

The production submitter was executed in `--check` mode with an
intercepted `sbatch` command.

The check completed successfully, no `sbatch` invocation occurred and
no `job_id.txt` file was generated.

The isolated submission copy matched the current operational source
files and its source SHA-256 manifest validated successfully.

The submitter itself does not execute the full preprocessing DAG during
`--check`; the Snakemake DAG validation recorded in this snapshot was
performed separately in the VERA environment using the repository
software runtime.

## Architecture

The ChIP-specific planning layer determines which runs enter
preprocessing. The existing generic preprocessing engine remains
unchanged.

In particular, no changes were required to:

- `workflow/rules/common.smk`;
- `workflow/scripts/build_validation_manifest.py`;
- `workflow/scripts/summarize_preprocessing_qc.py`;
- the generic FASTQ download implementation.

## Scope and limitations

This snapshot validates dynamic cohort propagation, FASTQ-manifest
construction, Snakemake wildcard expansion, the complete preprocessing
dry-run and safe submitter check mode.

It does not represent a production preprocessing execution.

No FASTQ files were downloaded during this validation, no production
SLURM job was submitted, and the full 22-run preprocessing workload has
not yet been executed.
