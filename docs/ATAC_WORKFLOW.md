# ATAC-seq through filtered BAM and TSS

## Entry points

The main `Snakefile` already includes `workflow/rules/atacseq.smk`.
`Snakefile.atac` uses the same reference, common and ATAC rules, with a separate
bulk validation list. It does not include metadata acquisition or RNA rules.
The general default `all` is unchanged; ATAC is an explicit modality target.

```bash
snakemake --snakefile Snakefile.atac atacseq_all --cores 8 --use-conda -n
# Equivalent modality target from the general workflow; may refresh metadata:
snakemake --snakefile Snakefile atacseq_all --configfile config/atacseq_bulk.yaml --cores 8 --use-conda -n
```

Remove `-n` to execute. From the main Snakefile, use the ATAC overlay as shown;
otherwise the historical mixed validation manifest is selected. The dedicated
entry point is recommended when using the curated catalogue already on disk.
`config/atacseq_validation_runs.tsv` selects SRR12774931 and SRR12774932.
Every requested run must pass the current shared selection policy. Liver run
SRR12774934 and multiome run SRR29929613 remain blocked. A successful technical
test does not restore eligibility. SINGLE/PAIRED is read structure, not cell
resolution. No currently retained bulk CHO PE library is available in this cohort.

## Workflow products

FASTQ manifest/checkpoint → FASTQ download with integrity checks → fastp →
Bowtie2 against nuclear plus mitochondrial reference → Picard → filtered BAM →
TSS profile and summary → `atacseq_validate_outputs` JSON.

`atacseq_filter_bam` handles both SINGLE and PAIRED; the former rule name
`atacseq_filter_single_bam` is retired. SE keeps the same MAPQ/flag criteria.
PE additionally requires proper pairs and mapped mates, groups candidates by
read name and retains groups with exactly two records, then coordinate sorts.
The validation rule checks one read1/read2 per name and reciprocal mate positions.
Picard and pair sorting use the job's temporary directory. Optional
`hpc_validation_cleanup=true` marks intermediate FASTQ and alignment BAM outputs
as temporary. It is false by default and never marks final BAM/TSS as temporary.

Validation scans the full BAM, checks its index, sort order, nonzero record count,
MAPQ, flags, mitochondrial exclusion and read layout. It checks TSS geometry,
finite nonnegative values, count totals, background, normalization and the maximum
score against the supplied profile. It hashes all four inputs and selection files.
It does not recalculate TSS from BAM or certify biological quality. The report
explicitly records `biological_acceptance: not_assessed`.

## Recheck existing Vera outputs without alignment

```bash
snakemake --snakefile Snakefile.atac atacseq_existing_all --cores 2 --use-conda \
  --config atacseq_existing_job=/cephyr/users/mayoa/Vera/atacseq_bulk_validation_results/launch_20260916T081250Z_2qbkz612/job_10297369 -n
```

This target schedules only two validation rules and their aggregate. It reads
persistent BAM/index and TSS directly, writing new reports under
`results/atacseq/qc/verification_existing`. It does not copy, overwrite or adopt
old BAMs as newly generated outputs and does not invoke alignment or download.
For full scans on Vera, submit `workflow/slurm/atacseq_verify_existing.sbatch`
with `ATAC_PROJECT_HOME`, `ATAC_EXISTING_JOB`, and `ATAC_VERIFY_SAVE`. This worker
uses scratch Conda environments and preserves reports/logs in the requested save
root. The BAMs and source repository are read-only. A pending Slurm job snapshots
code when it starts, so keep the selected worktree unchanged until then.

## Reproducible new bulk execution on Vera

```bash
python3 workflow/scripts/launch_atac_bulk.py --repo "$PWD"
# To submit a new job after reviewing the prepared snapshot:
python3 workflow/scripts/launch_atac_bulk.py --repo "$PWD" --submit
```

The launcher reads the versioned `Snakefile.atac`, rules, configuration and
`workflow/slurm/atacseq_bulk_validation.sbatch`; it does not embed a second copy
of the analysis. It freezes current source bytes, hashes and Git provenance.
It supports the documented two-run SINGLE cohort and refuses an active job or
an existing completed job under its output root. This is intentional: current
results already exist. It must not be used merely to install this integration.
It reserves 64 CPUs/96 GiB for the Vera scratch allocation, limits Snakemake to
16 threads and processes runs sequentially. Persistent copies are verified before
run-specific scratch cleanup. These settings reproduce the validated bulk launch;
they are not a resource benchmark or a portable scheduler profile.

## Evidence and limitations

`benchmarks/atacseq/2026-09-16/` contains the eligible bulk results, profiles,
original checksums, code snapshot and provenance from job 10297369. Large BAMs,
indexes, FASTQ, references and environments are not committed. The dated
historical launcher/worker preserve exactly how that job ran, not the new rules.
The integrated validation also passed on the real persistent BAMs and TSS outputs
for SRR12774931 and SRR12774932 in Vera job 10300242 (2026-09-16, COMPLETED 0:0,
00:00:53), using commit 0c74840255e558e1726b90cb55f5067737db41f1.
Only the two atacseq_validate_existing rules and their aggregate ran; alignment
and TSS calculation were not repeated. Reports, input hashes and execution
evidence are archived in benchmarks/atacseq/2026-09-16/verification_job_10300242.

`benchmarks/atacseq/2026-09-14/paired_technical/` preserves the three supplied PE
working files and patch on top of remote commit 354cf261. The old worker targets
an excluded multiome library and is historical only. Its run reports remain in
Vera; this package does not claim those reports were imported. PE support remains
technical; no eligible bulk CHO PE biological validation is claimed.

TSS scores are the maximum normalized 10-bp bin, not an ENCODE-equivalent score.
The two CHO libraries show high duplication and TSS enrichment near 9.4–9.7.
No universal threshold is applied. Peak calling, FRiP and tracks remain pending.

## Tests

In an environment containing pysam (the ATAC QC environment suffices):

```bash
python -m unittest discover -s workflow/tests -p test_atac_integration.py -v
```

Tests cover valid SE/PE BAMs; duplicate, low-MAPQ and mitochondrial rejection;
even-count orphan pairs; wrong layout; mate inconsistencies; selection exclusions;
actual TSS profiles and corrupted versions; preservation of an earlier report on
failed JSON writing. Additional delivery tests executed the actual filter shells
through Snakemake on synthetic records and checked fresh/existing-only DAGs.
Synthetic PE records are not evidence for eligibility of a biological PE library.
