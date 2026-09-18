# ChIP-seq Phase 6B dynamic core

The biological pairing source remains
`snapshots/chipseq/incremental_planning_validation_001/chipseq_analysis_plan.tsv`.
Only its 18 `ready` rows enter peak planning, and downstream code does not
reconstruct or optimize IP/Input relationships. Four Inputs are deliberately
shared, yielding 22 unique referenced runs.

The generated plan is `config/chipseq_peak_calling_plan.tsv`. It contains 8
narrow and 10 broad analyses. The six PRJNA865478 analyses use the
publication-supported fixed 147 bp policy. The 12 PRJEB9291 analyses require
their own formal PhantomPeakQualTools result and have no 147 bp fallback.
Formal single- and multi-candidate `estFragLen` records are supported. Candidate
fragment lengths and correlations must have equal, non-zero cardinality; the
documented first candidate is selected without applying a new quality rule,
while both complete lists and descriptive PhantomPeakQualTools fields are
retained. Each parsed record binds the analysis ID, IP run, exact BAM path and
SHA-256, and raw table path and SHA-256. These bindings and the raw table are
revalidated when parameters are generated and immediately before MACS3.

`prepare_chipseq_peak_calling_runtime.py --check` accepts only a successfully
completed dynamic filtering job referenced by
`results/chipseq/filtering/slurm/latest_submission.txt`. It verifies the
filtering manifest, BAM/CSI hashes, QC identities and counts, authoritative
plan and policy hashes, and mandatory reference provenance. The filtering copy
of the upstream reference provenance must exactly match the published copy and
its established nuclear, mitochondrial, combined-reference mapping, upstream
manifest, and FAI hashes are required. The FAI must contain formal five-field
records. The workflow derives MACS3 `-g` from
the validated FAI after excluding exactly one `NC_007936.1` contig. This is a
nuclear-reference-span approximation for MACS3 genome size, not a true
mappability-derived effective genome size.

Results are keyed by the stable upstream `analysis_id`. Narrow analyses
require narrowPeak and summit outputs. Broad analyses require broadPeak and
gappedPeak outputs and never depend on a summit file. Both use
`--keep-dup all`, because duplicate-marked reads were removed upstream.
Immediately before execution, the full MACS3 contract is reconstructed from
the runtime analysis, verified reference FAI, fixed study policy or bound
Phantom record, and BAM identities. A changed parameter record fails before an
external MACS3 process can start. Snakemake-owned JSON records use temporary
files and atomic replacement on legitimate reruns; MACS3 result directories
and external publication destinations retain collision refusal.
FRiP and Input-overlap values are descriptive and have no universal pass/fail
threshold. No CHO blacklist or validated replicate/IDR grouping is available,
so both limitations are recorded explicitly.

At this checkpoint no completed real 22-run dynamic filtering pointer exists.
Consequently real-cohort runtime preparation and any downstream production
execution correctly fail closed. Synthetic tests exercise the contracts; they
do not claim real-data validation.
