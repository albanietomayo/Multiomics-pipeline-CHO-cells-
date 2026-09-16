# ChIP-seq pilot duplicate marking and filtering

## Scope and entry point

This stage consumes the two validated SINGLE BAMs from alignment job 10297460:
SRR20770287 (input) and SRR20770297 (H3K4me3 IP). The preceding alignment code
and validation evidence were published at commit 96a13edf8772bf30cc2be588862f42be8642541c.
No realignment or reference download is required.

Run `python3 workflow/slurm/submit_chipseq_filtering.py --check` for a source
snapshot and provenance check; use `--submit` to also submit that snapshot to
SLURM. The submitter records the actual base commit, working tree status, source
files, SHA-256 hashes, settings and generated input plan. Each submission is
isolated; existing alignment outputs are read-only inputs.

The Snakemake entry point is `workflow/rules/chipseq_filtering.smk`, target
`chipseq_filtering_all`. Shared experimental eligibility is a dependency of BAM
staging. Staging verifies each complete BAM against the alignment output
manifest and validates the reference sequence dictionary and read group. The
manifest is anchored to `snapshots/chipseq/alignment_validation_001`.

## Duplicate policy

Picard MarkDuplicates 3.5.0, separate environment with OpenJDK 17, operates on
coordinate-sorted, unfiltered SINGLE BAMs. Duplicates are marked without removing
records, using SUM_OF_BASE_QUALITIES. Input and IP are processed independently.
This is coordinate-based duplicate estimation; it does not establish that every
marked record is a PCR artifact. No UMI-based interpretation is attempted.

Optical duplicate detection is explicitly disabled (`READ_NAME_REGEX=null`):
optical coordinates in these deposited read identifiers have not been validated.
No claim of optical-duplicate rate or reliable library-size estimation is made.
Picard metrics and logs are retained.

## Filtering policy and accounting

Retain nuclear primary alignments with MAPQ >= 30 and < 255, excluding unmapped,
QC-failed and duplicate-marked records. Flags excluded: 3844 = 4 + 256 + 512 +
1024 + 2048. This is the configured pilot analysis policy, not a universal
biological eligibility criterion. MAPQ 255 means unavailable mapping quality;
it is explicitly excluded. Both strands are retained. Mitochondrial reference:
NC_007936.1. Other nuclear scaffolds are retained. No blacklist is applied.

`filtering_flow.tsv` assigns each rejected record to the FIRST applicable reason:

1. Secondary or supplementary alignment.
2. Unmapped.
3. QC failure flag.
4. Mitochondrial alignment.
5. MAPQ unavailable (255).
6. MAPQ below 30.
7. Duplicate mark.

Input count = sum of sequential discards + retained count. `filtering_qc.json`
also reports independent criterion counts, which can overlap and MUST NOT be
summed. Thus the duplicate count in the final sequential step can be lower than
Picard's total duplicate count. Picard PERCENT_DUPLICATION is a fraction (0-1),
whereas `retained_pct_of_input` is a percentage (0-100).

The marked BAM is scanned fully. Total, mapped and nuclear MAPQ>=30 counts must
match the validated upstream alignment reports; marked mapped-primary duplicate
counts must match Picard. The filtered BAM is indexed with CSI and independently
read in full to check record count and filtering invariants. Zero retained reads
cause failure before downstream analysis. QC includes flagstat, idxstats and
samtools stats; signal quality remains `not_evaluated`.

## Storage, execution and outputs

SLURM requests 4 CPUs, 24 GiB RAM, four hours. Picard Java heap is 12,000 MiB;
the scheduler permits one duplicate-marking job at a time. Conda environments,
copied input BAMs, duplicate-marked BAMs and sorting scratch use node TMPDIR.
Copied inputs and duplicate-marked BAMs are temporary Snakemake outputs.

Only filtered BAMs, CSI indexes, reports, parameters, reference provenance and
hash-verified input descriptions are published under:
`results/chipseq/filtering/slurm/submission_*/job_*/outputs/`.
The original raw BAMs remain in the alignment submission. Filtered output copies
are SHA-256 checked before the job is marked completed. If execution fails,
available small reports are preserved; large intermediate BAMs are not copied.
No automatic deletion of existing persistent results is performed.

## Validation and limitations

Local regression tests use small synthetic BAMs, including overlapping discard
criteria, MAPQ 29/30/255, both strands, mitochondria, duplicate-count mismatch,
changed reference dictionary, upstream count mismatch and altered BAM hashes.
Run `python workflow/tests/test_chipseq_filtering.py` in the filtering environment.
On 2026-09-16, all nine regression tests passed with pysam 0.24.0. A
Snakemake 9.25.2 dry-run against synthetic placeholder inputs scheduled all 13
jobs, including shared eligibility, staging, duplicate marking, filtering,
indexing, final BAM validation and the summary. Python and Bash syntax checks
also passed. Picard itself was not executed locally.

Local testing does not establish real-data duplicate rates or successful VERA
Conda resolution. Real execution, result review, validation snapshot, commit and
push remain pending until confirmed by a successful VERA job.

This is a pilot-specific downstream entry point. It validates frozen curated
eligibility and previously generated inputs; it does not yet create missing
alignment or preprocessing outputs automatically. Incremental eligibility rules
for new samples, paired-end processing, a unified ENA-to-peaks DAG, peak calling
and signal QC remain pending. No ATAC-specific Tn5 shift is applied.

## Primary documentation

- https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard
- https://www.htslib.org/doc/samtools-flags.html
- https://www.htslib.org/doc/samtools-view.html
