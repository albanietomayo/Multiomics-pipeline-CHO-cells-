# ChIP-seq Phase 6C operational validation

Validation date: 2026-09-18 on VERA. Starting checkpoint:
`314093950e3110a53822f4982249c687884624c5` (`Implement dynamic ChIP-seq
Phase 6B core`) on `chipseq-incremental-generalization` with a clean worktree.
No SLURM job, real PhantomPeakQualTools job, MACS3 job, or real 22-run workflow
was started.

## VALIDATED NOW

### Biological and reference contracts

The 18 `ready` rows in
`snapshots/chipseq/incremental_planning_validation_001/chipseq_analysis_plan.tsv`
and the generated peak plan have exactly identical `(analysis_id, study, IP,
Input)` tuples. The plan contains 18 IP runs, four shared Input runs, and 22
unique runs: 12 PRJEB9291 and six PRJNA865478 analyses; eight narrow and ten
broad analyses. Six PRJNA865478 rows use fixed 147 bp and all 12 PRJEB9291
rows require analysis-specific PhantomPeakQualTools estimates. No PRJEB9291
row contains a fixed or 147-bp value.

The runtime validation retains the single GCF_003668045.3 /
CriGri-PICRH-1.0 reference contract, requires NC_007936.1, and derives the
MACS3 `-g` value from the verified FAI after excluding that mitochondrial
contig. This remains a nuclear-reference-span approximation, not a
mappability-derived effective genome size.

### VERA software and environment resolution

The bare shell provided Python 3.9.25. Available and executed modules were:

- Miniforge3/24.1.2-0: Conda 24.1.2 and Mamba 1.5.7.
- snakemake/8.27.0-foss-2024a: Snakemake 8.27.0 and Python 3.12.3.
- R/4.4.2-gfbf-2024a: R and Rscript 4.4.2.
- SAMtools/1.21-GCC-13.3.0: samtools/htslib 1.21.

VERA also advertised Snakemake 8.4.2 and R 4.3.2/4.4.1. MACS3, bedtools,
PhantomPeakQualTools, and `run_spp.R` were not present in the bare shell or
found as targeted modules. They are resolvable through the repository Conda
environments.

The unchanged `workflow/envs/chipseq_phantompeakqualtools.yaml` solved and was
created on linux-64 with PhantomPeakQualTools 1.2.2, samtools/htslib 1.24,
R 4.4.3, spp 1.16.0, and Rsamtools 2.22.0. The unchanged peak environment
solved and was created with Python 3.10.21, MACS3 3.0.4, samtools 1.24, and
bedtools 2.31.1. The main environment also solved and was created with its
pinned Python 3.14.7, Snakemake 9.25.2, and Conda 26.7.2. All validation
environments were under `/tmp`, outside the repository.

Snakemake 8.27.0 can parse the workflow, but its `--use-conda` mode rejects
the Miniforge module's Conda 24.1.2 because it requires Conda 24.7.1 or newer.
The production-pinned Snakemake 9.25.2/Conda 26.7.2 combination succeeds.
The SLURM wrapper now sets strict Conda channel priority to remove the warning
and avoid user-configuration-dependent channel mixing.

### Installed `run_spp.R` interface

The two independently created PPQT environments contained byte-identical
`run_spp.R` scripts (SHA-256
`871441dd7a94f703578217dcceeb7799c2d7a9efa9a42390b64adce8084c342f`).
This installed 1.2.2 script and the upstream project documentation agree on:

- `-c=<ChIP BAM>` for input, `-p=<nodes>` for parallel workers,
  `-savp=<PDF>` for the cross-correlation plot, and `-out=<file>` for the
  tab-delimited result;
- 11 result fields: filename, read count, `estFragLen`, `corr_estFragLen`,
  phantom peak and correlation, minimum-correlation shift and value, NSC,
  RSC, and quality tag;
- up to three comma-separated fragment candidates within 90% of the maximum,
  ordered by decreasing correlation, with correlations in the same order;
- the first candidate as the predominant fragment-length estimate in almost
  all cases. The installed source sorts by decreasing correlation and assigns
  the first candidate to its primary cross-correlation peak.

Therefore Phase 6B's first-candidate selection is consistent with the formal
tool behavior. NSC, RSC, and quality tag remain descriptive and do not create
new pipeline PASS/FAIL thresholds. `run_spp.R --help` is not a supported flag
in this build (it exits 1); its embedded usage text and installed source were
inspected instead. Cross-correlation is called with `accept.all.tags=TRUE`;
there is no separate duplicate-removal flag in this invocation. That is
consistent with the pipeline's already-deduplicated BAM contract.

Because `-out` appends and an existing `-savp` output can abort the tool, the
Snakemake rule now writes both outputs to per-process temporary paths, checks
that both are non-empty, and replaces only Snakemake-owned outputs after a
successful run. Published job outputs retain collision refusal.

### Snakemake, preflight, and submission safety

Using a synthetic runtime with the exact real 18-analysis plan, both
Snakemake 8.27.0 and the pinned Snakemake 9.25.2 parsed and dry-ran the graph.
The pinned stack also succeeded with `--use-conda` after Snakemake-only
creation of both rule environments. The DAG contained 98 jobs: 12 PPQT, 12
PPQT parsers, 18 parameter, 18 MACS3, 18 FRiP, 18 per-analysis QC, one cohort
summary, and one all-target job. No biological command executed.

The real configured pointer,
`results/chipseq/filtering/slurm/latest_submission.txt`, is absent. The direct
preflight and `submit_chipseq_peak_calling.py --check` therefore exit nonzero
with `No real productive dynamic filtering pointer exists; peak calling
remains fail-closed`. Regression tests also verify rejection of an incomplete
or failed job, missing output manifest, checksum changes, missing BAM or CSI,
wrong run roles, and mismatched reference provenance.

Submission requires the explicit, mutually exclusive `--submit` action.
`--check` returns before the sole `sbatch` call site, and failed local/preflight
validation occurs before creation of a submission directory. Imports have no
submission side effect. The source bundle is checksummed before execution;
base commit and worktree status are recorded. Each submitted job uses a unique
directory, MACS3 refuses pre-existing result files, publication refuses an
existing output directory/manifest, and Snakemake uses per-analysis outputs
with `--rerun-incomplete` restart boundaries.

### Storage and resource evidence

The repository occupied 36 MiB; no real 22-run BAM/CSI or FASTQ cohort was
present. The authoritative processing plan records 39,120,498,415 FASTQ bytes
(39.12 GB decimal; 36.43 GiB) for the 22 unique runs, but this does not predict
filtered BAM, PPQT scratch, bedGraph, or peak-output sizes. At inspection time
the shared `/cluster` filesystem reported 6.3 PB free and login-node `/tmp`
reported 710 GB free. These are filesystem observations, not a verified user
quota or compute-node `$TMPDIR` guarantee.

Read-only accounting for recorded historical ChIP jobs showed only metadata,
eligibility, and a two-run alignment pilot. The alignment job completed in
01:16:20; its `/usr/bin/time` report recorded a 6,705,104 KiB maximum RSS.
There is no PPQT or MACS3 accounting evidence. Consequently the current
four-CPU, 32-GiB, four-hour peak-calling request is uncalibrated and must not
be represented as validated.

The DAG already separates each of 12 PPQT estimates and each of 18 MACS3 calls
into restartable files, and shared controls are referenced rather than copied
per analysis. However, all stages currently run inside one SLURM allocation.
A monolithic production job is not operationally justified until representative
PPQT and MACS3 resource measurements establish walltime, RAM, and node-local
scratch needs. Stage-specific allocations can be considered after that
benchmark without changing biological pairings or fragment policy.

## PENDING REAL-HPC VALIDATION

- Produce and validate the completed productive dynamic filtering submission
  for all 22 runs and its `latest_submission.txt` pointer.
- Verify the real BAM/CSI/reference hashes and FAI-derived nuclear span through
  the real preflight.
- Run representative PPQT and MACS3 benchmarks in an authorized compute
  allocation; PPQT execution on the real PRJEB9291 BAMs was intentionally not
  attempted on the login node.
- Measure walltime, peak RSS, node-local scratch use, and generated bedGraph/
  peak storage before approving resource requests or a monolithic allocation.
- Confirm project quota and compute-node `$TMPDIR` capacity for the measured
  workload.
- Review real PPQT plots and any multi-candidate estimates. The documented
  first candidate is tool-consistent, but the upstream tool documentation
  recommends visual review when profiles contain multiple strong peaks.
- Execute an authorized real-data dry-run after the productive filtering
  dependency exists. Do not submit the 18-analysis production workflow until
  all preceding STOP conditions are cleared.

Current operational disposition: `SAFE_FOR_REAL_SLURM=NO`. The code,
regression coverage, resolved environments, and synthetic dry-run are coherent
for a Phase 6C checkpoint after human review; no real-data execution claim is
made.
