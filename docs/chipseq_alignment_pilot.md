# ChIP-seq pilot alignment

Scope: one reviewed SINGLE-END IP/input pilot from `config/chipseq_pilot.json`.
Entry point: `workflow/rules/chipseq_alignment.smk`.
Launcher: `python3 workflow/slurm/submit_chipseq_alignment.py --check|--submit`.

The launcher snapshots source and configuration and checks the existing eligibility
validator. It verifies upstream QC report hashes and completed job status. On the
compute node, processed FASTQs are copied from the selected preprocessing attempt
and their original SHA-256 hashes are verified. No FASTQ download is scheduled.

The workflow reuses `fetch_reference_genome.py` and its environment to acquire the
configured NCBI assembly. A separate rule streams the nuclear FASTA and a versioned
mitochondrial sequence into a mapping reference, checks duplicate sequence IDs,
and records hashes. An existing mitochondrial ID causes a stop for explicit review.
The original fetcher checks assembly identity and NCBI MD5 values. It records the
configured annotation release but does not independently verify that release;
annotation is not used by these alignment rules.

Bowtie2 2.5.5 runs in very-sensitive-local mode with fixed seed 0. This allows
soft clipping while retaining the complete read sequence in the BAM. The choice
is recorded and does not establish that the observed composition bias is harmless.
Samtools 1.24 sorts, creates a CSI index and reports flagstat, idxstats and stats.
The custom QC pass reads the complete BAM, checks primary read totals against
the processed FASTQ QC totals and reports mitochondrial mapping, MAPQ>=30 counts,
and 5'/3' soft clipping relative to the original read orientation. MAPQ is measured,
not used to filter this BAM. No duplicate marking, BAM filtering, Tn5 shift, peak
calling, enrichment evaluation or biological replication inference is performed.

## Pilot preprocessing review

SRR20770287 input: 43,830,994 reads; SRR20770297 H3K4me3 IP: 41,724,399.
FastQC base-quality modules pass; final-position mean quality is approximately Q34.
About 97% of reads have length 82-84 nt. Both samples have composition bias at
the start and A=0% at position 84. The raw and processed composition profiles are
identical. The cause of the terminal bias remains unresolved. The fastp command
uses -Q -L -G: quality/length filtering and polyG trimming are disabled; adapter
trimming is enabled, with no reported trimming in these two runs. Zero low-quality
filter counts are not evidence that an enabled quality filter was passed.
Proceed to pilot alignment without additional fixed trimming; evaluate mapping
and clipping before deciding downstream filtering or drawing signal conclusions.

## Storage and provenance

SLURM requests 8 CPUs, 32 GiB RAM and 6 hours. Reference, index, environments and
FASTQ copies live in node-local TMPDIR. Preserve the sorted BAMs/CSI, QC/logs,
the exact combined reference FASTA compressed with gzip, its FAI and provenance,
source snapshot and explicit Conda package lists. Each saved output is hashed
before copy and verified after copy. Job success is reported only after copying.
The original nuclear SHA manifest is retained as provenance; its annotation paths
refer to the temporary build and are not included in the exported artifact set.
The output SHA manifest covers the actual saved files. The Bowtie2 index is not
persisted in this pilot: new attempts rebuild it. Shared reference/index caching
remains future work, subject to storage capacity and exact reference hashes.

Each submission is a distinct attempt. --check does not schedule work. Repeating
--submit creates another job, so use latest_submission.txt to monitor an attempt.
Copy failures mark the job failed; partial output may remain in its own job folder.
Successful scheduler execution does not imply acceptable biological ChIP signal.

## Pending generalization

This entry point deliberately rejects paired or multi-file inputs and only uses
the selected pilot. It does not implement automatic eligibility for new ENA runs,
incremental catalog updates or production-wide ChIP processing. The rule-based
eligibility generalization requested in docs/chipseq_pending_work.md remains pending.

References: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
and https://www.htslib.org/doc/samtools.html .

## Development validation (2026-09-16)

Python syntax and Bash syntax checked. A Snakemake 9.25.2 dry-run with synthetic
inputs scheduled 13 jobs, including the shared eligibility rule, reference and
index preparation, both alignments and QC. Targeted tests verified forward/reverse
soft-clip orientation, mitochondrial/unmapped counts, upstream report hashes,
corrupt FASTQ rejection, streaming reference assembly and duplicate-ID rejection.
No real CHO reference download or BAM alignment was executed in this development
environment. VERA run validation and scientific interpretation remain pending.
