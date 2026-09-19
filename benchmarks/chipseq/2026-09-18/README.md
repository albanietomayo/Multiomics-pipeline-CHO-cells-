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

One explicitly selected run is processed per SLURM job: acquisition with
size/MD5 validation, fastp with the production SINGLE settings (`4` threads,
`-Q -L -G`), one-thread post-fastp FastQC, Bowtie2/Samtools alignment against
the harmonized CriGri-PICRH-1.0 / GCF_003668045.3 plus NC_007936.1 reference,
and the validated ChIP duplicate-marking/nuclear-filtering policy.

## Exact artifact and duplication audit

| Artifact | Created by | Node TMPDIR copy | Persistent successful-run copy | Duplicate persistent copies | Retention reason | Regenerability |
|---|---|---:|---:|---:|---|---|
| nuclear FASTA, GFF3, GTF, sequence report, metadata | reference provisioning outside benchmark | none | shared root only | zero per run | immutable reference and annotation identity | regenerable from fixed NCBI assembly/provenance |
| combined nuclear/mitochondrial FASTA + FAI | shared reference provisioning | none | shared root only | zero per run | exact mapping reference and verified nuclear span | regenerable from the verified nuclear FASTA and versioned mitochondrial accession |
| six Bowtie2 large-index components | shared reference provisioning | none | shared root only | zero per run | alignment prerequisite | regenerable from the hashed combined FASTA and recorded command |
| downloaded raw FASTQ | acquisition | one | none | zero | stage input; expected bytes and MD5 remain in manifest | public ENA URL/accession makes it regenerable |
| processed FASTQ | fastp | one | none | zero | alignment input | regenerable from raw FASTQ, exact fastp parameters and environment |
| alignment-stage FASTQ | old harness alignment wrapper | none in lean harness | none | zero | eliminated implementation-only copy | not applicable |
| raw sorted BAM + CSI | alignment | one | none | zero | filtering input; size/hash and alignment QC retained | regenerable from processed FASTQ, reference/index and parameters |
| filtering-stage raw BAM copy | filtering input validation | one | none | zero | temporary validation input | regenerable from raw alignment BAM during a clean rerun |
| duplicate-marked BAM | Picard MarkDuplicates | one | none | zero | filtering input; byte/hash and Picard metrics retained | regenerable from raw BAM and exact Picard settings |
| filtered BAM + CSI | filtering | one until verified copy | exactly one | zero | reusable scientific product | retained |
| fastp/FastQC/alignment/filtering reports | their respective stages | one while running | exactly one allowlisted copy | zero | scientific interpretation and benchmark evidence | retained |
| stage timing/resource records and logs | harness | no large scratch copy | direct persistent write | zero | wall/CPU/MaxRSS/scratch/input-output evidence | retained |
| source snapshot and submission plan | submitter | copied into TMPDIR only as executable source | one submission snapshot | zero per job result | exact code/config/manifest/Git provenance | retained |

The old worker recursively copied the entire work tree after every stage and
again on exit. That made the raw and processed FASTQs, alignment FASTQ copy, raw
BAM, filtering input copy, duplicate-marked BAM, reference/annotation and
Bowtie2 index persistent implementation conveniences. The lean worker never
recursively persists its work tree. Every successful stage first records each
declared output's logical bytes, SHA-256, retention class and regenerability in
a compact persistent stage record. The scheduler may then discard TMPDIR.

## Persistent result contract

Only the final filtered BAM/CSI, fastp JSON/HTML, post-fastp FastQC HTML/ZIP,
alignment reports, Picard/filtering reports, generated input-provenance plans,
stage metrics/logs, software versions and explicit Conda environment exports,
shared-reference inventory, logical persistent-footprint measurement, source
snapshot/manifest/commit, job status and final completion record persist.
The final completion record is atomic and is created only after every
allowlisted copy has passed size and SHA-256 validation. A failure writes no
completion record and preserves small diagnostics already written to the job
directory when practical.

Per-stage metrics contain UTC times, wall seconds, GNU-time user/system CPU and
MaxRSS, exit status, input/output logical bytes, and scratch
baseline/peak/increment. These are captured before transient artifacts
disappear. `collect_chipseq_benchmark_sacct.sh` separately retains scheduler
JobID, state, elapsed time, TotalCPU, allocated CPUs, MaxRSS, AveRSS and exit
code.

## Shared reference prerequisite

No canonical complete ChIP reference/index was found in the inspected project
area, so the submitter requires `--shared-reference-root`; there is no default.
Before either check or submit it requires nonempty nuclear FASTA, GFF3, GTF,
sequence report, reference metadata and SHA manifest; combined FASTA,
mitochondrial FASTA, FAI and ChIP provenance; and all six
`genome_plus_mt.*.bt2l` components. It validates CriGri-PICRH-1.0 /
GCF_003668045.3 / GCA_003668045.2, annotation release 104, NC_007936.1,
the FAI-derived nuclear span, file sizes and SHA-256 values. The immutable
inventory includes exact execution paths and regeneration commands. The worker
recomputes it against the submission plan before processing.

The safe submitter defaults to check-only behavior and cannot invoke `sbatch`
without `--submit`. Submission also requires the expected branch, a clean
worktree, exact regenerated manifest, adequate output/scratch minimums and no
pre-existing class/accession output. A validated completed job is never silently
overwritten.

## Restart and storage consequences

Large-stage checkpoint reuse is intentionally removed. A node failure or clean
rerun starts again from the public ENA object; this is scientifically
reproducible from accession, URL, expected bytes/MD5, reference and index
hashes, exact software environments, parameters and Git commit. Small failure
diagnostics remain, while the submission/run collision guards require an
operator to choose a new isolated attempt instead of overwriting evidence.

For P10, the two committed historical pilot observations give final-filtered
BAM/source-FASTQ byte ratios of 0.8814 and 0.8449. Treating their midpoint only
as an empirical engineering estimate (not a biological scaling law) gives a
1.026 GiB filtered BAM for the 1,273,744,805-byte P10 input. Adding measured
small QC/source evidence and modest metadata/index overhead gives a 1.04 GiB
central persistent estimate. Applying the larger observed ratio plus a
conservative evidence allowance gives 1.15 GiB. Against 9.847 GiB free, the
corresponding margins are approximately 8.807 GiB and 8.697 GiB. Scratch
expansion remains unknown until a real P10 run and is not represented by these
persistent estimates.

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
  --shared-reference-root /ABSOLUTE/PATH/TO/CriGri-PICRH-1.0 \
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
