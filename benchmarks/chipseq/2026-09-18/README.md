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
| compressed combined nuclear/mitochondrial FASTA | historical reference run | none | existing external object; linked, not copied | zero per run | hashed archival source and index regeneration | public source is regenerable; exact compressed and content hashes are pinned |
| combined-reference FAI | historical reference run | one 24,707-byte TMPDIR copy for filtering | existing external object; linked, not copied | zero per run | exact contig dictionary and verified nuclear span | regenerable from the pinned compressed FASTA |
| reference provenance + immutable inventory | historical run / one-time provisioning | none | one small canonical copy plus existing provenance link | zero per run | assembly, source, build and component identity | regenerable from validated evidence, except that hashes are retained as the identity contract |
| six Bowtie2 2.5.5 large-index components | one-time provisioning | none during benchmark | shared root only | zero per run | alignment prerequisite | regenerable once from the pinned compressed FASTA; never rebuilt per run |
| nuclear FASTA, mitochondrial FASTA, sequence report, metadata and source-hash manifest | historical reference acquisition | none | existing evidence only | zero | provenance only; not consumed by benchmark stages | regenerable from fixed accessions and hashes |
| GFF3 and GTF | historical reference acquisition | none | not required by this benchmark contract | zero | no acquisition, QC, alignment, duplicate, filtering or metric consumer | regenerable from fixed assembly/annotation provenance |
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

The benchmark has no default reference path. Its minimal canonical root contains
only `genome_plus_mt.fa.gz`, `genome_plus_mt.fa.fai`,
`reference_provenance.json`, `reference_inventory.json`, and the six
`bowtie2_index/genome_plus_mt.*.bt2l` files. The three historical source files
may be symlinks to the validated external objects; every check follows them and
verifies exact bytes and SHA-256 values. The uncompressed combined FASTA is not
persistent. GFF3, GTF, standalone nuclear/mitochondrial FASTAs, sequence report,
metadata and the old hash list are not operational inputs.

Before check or submission, and again in the worker before biological
processing, validation requires the exact CriGri-PICRH-1.0 /
GCF_003668045.3 / GCA_003668045.2 identity, configured annotation release 104,
NC_007936.1, a FAI-derived nuclear span of 2,366,634,374 bp, mitochondrial
length 16,284 bp, the three validated historical source hashes, all six index
components, Bowtie2 2.5.5 large-index build metadata, and an exact match to the
immutable inventory. The submission plan then pins the inventory itself and all
execution paths, sizes and hashes for the worker recheck. Annotation release 104
is retained as configured provenance and is not misrepresented as independently
verified.

The safe submitter defaults to check-only behavior and cannot invoke `sbatch`
without `--submit`. Submission also requires the expected branch, a clean
worktree, exact regenerated manifest, adequate output/scratch minimums and no
pre-existing class/accession output. A validated completed job is never silently
overwritten.

### Existing validated evidence

The reusable source directory is
`/cephyr/users/mayoa/Vera/TFM_multiomics_pipeline_chipseq/results/chipseq/alignment/slurm/submission__niyvyr2/job_10297460/outputs/reference`.
Direct inspection gave:

| File | Bytes | SHA-256 | Contract role |
|---|---:|---|---|
| `genome_plus_mt.fa.gz` | 882,813,762 | `f5e1effa4d063b9005eeaa0f245e7a09e84d970d6efaeb80f442cfbbf6216923` | linked archival source |
| `genome_plus_mt.fa.fai` | 24,707 | `2dfb28d82459be6cbddb6e4be79e2c328c436654cc37f2e703c7abb9c50a0d57` | linked runtime dictionary |
| `reference_provenance.json` | 710 | `3646c60b547d946814704af383464807dda890cfd371ff78f24edf8307fd582f` | linked runtime provenance |
| `mitochondrial.fa` | 16,580 | `b7ccf6b1c6981c2b4a0c9e57245bd6712a861e685571db56d274c04fc30ed37c` | provenance only |
| `sequence_report.jsonl` | 195,537 | `e96372c787e1422cabd813efd8a2c2a7981b5f0136605b7531813a9beb318332` | provenance only |
| `reference_metadata.tsv` | 494 | `eb752199fe08accc16b561199a2092c1b4dae4778ee24183167649edf0456e07` | provenance only |
| `original_nuclear_resources.sha256` | 601 | `669b4844f084ad51c5863931bae7fbe4585f3fee2dfbfaa300c583aec5d2464f` | provenance only |

The relevant project trees contain no surviving Bowtie2 component and no
standalone nuclear FASTA. The latter remains represented inside the validated
compressed combined source; its original hash, plus the GFF3/GTF and sequence
report hashes, remains in `original_nuclear_resources.sha256`.

### One-time provisioning plan (do not execute on a login node)

This is a compute workload: the historical build used about 6.7 GB maximum RSS
and roughly 15 minutes just for the forward/reverse index build. Provision it as
one explicitly authorized SLURM job, not as a benchmark stage and not once per
benchmark run. The job should:

1. Refuse an existing canonical target and create a staging directory on the
   same shared filesystem.
2. Symlink the exact existing `genome_plus_mt.fa.gz`, FAI and provenance files
   from the historical output listed above.
3. Verify their recorded hashes, decompress the FASTA only to node-local
   `$TMPDIR`, and verify the decompressed SHA-256
   (`dfd445e136c4bc9c11f9616d251b86732d1aab0cbdbf9d1ababb8093f7430e94`).
4. Run Bowtie2 2.5.5 exactly once with `--large-index --threads 8`, writing all
   build output under `$TMPDIR`.
5. Copy each completed component to a `.part` name in the shared staging tree,
   verify size/hash, then rename it. Generate `reference_inventory.json` with
   `chipseq_benchmark.py create-reference-inventory`.
6. Run `check-reference` against the staging tree, remove write permissions from
   the payload/inventory, atomically rename the complete staging directory to
   the previously absent canonical root, and run `check-reference` once more.
7. Remove only node-local temporary files through the scheduler's TMPDIR cleanup.

The whole-directory final rename prevents benchmark jobs from observing a
partial index. No permanent uncompressed FASTA and no annotation copy are
created.

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
1.026 GiB filtered BAM for the 1,273,744,805-byte P10 input. The conservative
successful-run persistence allowance remains 1.15 GiB (1,234,803,098 bytes when
rounded upward).

On 2026-09-19, the final Ceph snapshot reported 23,005,142,471 bytes used of a
32,212,254,720-byte home quota: 9,207,112,249 bytes (8.574791 GiB) free.
The historical reference
output already counted in usage is 883,052,391 bytes; the three objects reused
by the operational contract account for 882,839,179 of those bytes. No index
component survives in either relevant ChIP alignment-results tree.

The historical Bowtie2 log records exact sizes of 796,690,238 bytes for each of
`.1.bt2l` and `.rev.1.bt2l`, and 1,182,112,388 bytes for each of `.2.bt2l` and
`.rev.2.bt2l`: 3,957,605,252 bytes (3.685807 GiB) known new data. It does not
record `.3.bt2l` or `.4.bt2l` sizes, and the new inventory size depends on the
new component hashes. Therefore the exact new-reference total and exact final
quota margins remain unknown until the one-time build. Based only on the known
four files, free space after provisioning is strictly less than 5,249,506,997
bytes (4.888984 GiB), and after the conservative P10 allowance strictly less
than 4,014,703,899 bytes (3.738984 GiB). These are upper bounds on remaining
space, not forecasts. Scratch expansion remains unknown until a real P10 run.

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
