# ChIP-seq P10/P50/P90 resource benchmark

## IMPLEMENTED / VALIDATED WITHOUT EXECUTION

This directory defines a deterministic upstream resource benchmark for the 22
unique SINGLE-end Illumina runs in the production ChIP-seq cohort. The fixed
P10/P50/P90 representatives span the observed compressed-input distribution:

| Class | Run | Study | Compressed bytes |
|---|---|---|---:|
| P10 | ERR868176 | PRJEB9291 | 1,273,744,805 |
| P50 | ERR868152 | PRJEB9291 | 1,678,499,646 |
| P90 | SRR20770294 | PRJNA865478 | 2,193,948,334 |

`benchmark_selection.tsv` contains only the fixed class/accession choices.
`chipseq_benchmark_manifest.tsv` is generated from `config/samples.tsv` by
`workflow/scripts/chipseq_benchmark.py generate-manifest`; URLs, byte counts,
MD5 values, layout, platform, study, source hash, and Git commit are inherited
and then checked byte-for-byte against a fresh deterministic derivation.

One explicitly selected run is processed per SLURM job through restartable
boundaries: acquisition with size/MD5 validation, fastp with the production
SINGLE settings (`4` threads, `-Q -L -G`), one-thread post-fastp FastQC,
Bowtie2/Samtools alignment against the harmonized CriGri-PICRH-1.0 /
GCF_003668045.3 plus mitochondrial reference, and the validated ChIP duplicate
marking/nuclear filtering policy. Shared reference acquisition and index setup
are recorded separately and are not attributed to a run-size stage.

Each successful stage gets an atomic completion record containing output sizes
and SHA-256 values. A resumed/requeued worker validates those records before
skipping work. Raw FASTQ, processed FASTQ, raw BAM, duplicate-marked BAM,
filtered BAM, QC outputs, logs, and inventories are retained; the harness has
no automatic data cleanup.

Per-stage output includes UTC timestamps, wall seconds, GNU-time CPU and MaxRSS
values, exit status, input/output bytes, and sampled scratch use. Software
versions, source commit, manifest hash, source hashes, reference setup costs,
and stage status are retained. `collect_chipseq_benchmark_sacct.sh` requests
JobID, JobName, State, Elapsed, TotalCPU, AllocCPUS, MaxRSS, AveRSS, and ExitCode;
the parser writes `not_available` for empty optional values and rejects missing
columns or malformed input.

The safe submitter defaults to check-only behavior. It requires one class and
optionally accepts the matching accession as a redundant check. Check mode
cannot invoke `sbatch`; `--submit` is mandatory for submission. Submission
requires the expected branch, a clean worktree, an exact regenerated manifest,
adequate storage for at least the known compressed input, and no pre-existing
class/accession output directory. `--development-dirty-check` is limited to
non-submitting validation.

The storage preflight reports known compressed-input bytes and available output
and scratch bytes. It rejects either filesystem if it cannot hold even that
known minimum. Downstream expansion remains explicitly unknown until this
benchmark runs; no unsupported free-space multiplier or resource conclusion is
encoded. Operators must review the reported space before real submission,
especially because preserving all intermediates may exceed the previously
observed roughly 14 GiB in the VERA area.

## Intentionally separate

PhantomPeakQualTools and MACS3 are excluded from this upstream run-size
benchmark. They require paired biological analyses and study-specific fragment
handling: PRJEB9291 requires its own valid PPQT estimate with no fallback,
whereas PRJNA865478 uses the publication-supported fixed 147 bp value. They
must be benchmarked later as a separate downstream analysis-level exercise so
these policies cannot contaminate input-size scaling.

The RNA-seq benchmark contributed useful engineering patterns—isolated classes,
scratch sampling, GNU time, inventories, provenance, and machine-readable
metrics. Its paired/mixed layout assumptions, STAR/counting target, cleanup,
measured resource conclusions, and all RNA-derived scaling models were not
reused.

## PENDING REAL SLURM BENCHMARK

No resource row is complete and no resource value is claimed. After human
review and commit, run a check from the repository root (this performs no
submission):

```bash
python workflow/slurm/submit_chipseq_benchmark.py \
  --benchmark-class P10 \
  --run-accession ERR868176
```

Repeat with the matching P50/P90 selection. Only after reviewing storage and
the provisional bootstrap envelope may an authorized operator add `--submit`.
Do not submit all three at once until storage has been reviewed. After each job
finishes, collect accounting read-only:

```bash
bash workflow/slurm/collect_chipseq_benchmark_sacct.sh JOB_ID OUTPUT_DIRECTORY
```

Actual stage and `sacct` observations will be aggregated to define defensible
CPU, memory, walltime, and scratch classes for the 22-run cohort. PPQT/MACS3
resource classes will be derived independently from analysis-level evidence.
