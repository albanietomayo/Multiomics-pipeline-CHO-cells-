# Clean-clone reproduction validation

Date: 2026-09-30. This is evidence for the prospective, unstaged v1.0.1 patch,
not a release or replacement of historical scientific outputs.

Only **one of the 18 PASS ChIP-seq analyses** was re-executed from FASTQ:
SRR20770297 IP versus SRR20770287 Input, PRJNA865478, fixed 147 bp.
**The remaining 17 PASS analyses were not re-executed from FASTQ.** All
reproduction findings below are restricted to this single validation analysis.

For this analysis, historical peak count, FRiP and peak/summit coordinates
match. Strict signal identity fails. Substituting only this analysis's
reproduced products into the frozen downstream aggregation preserves its locus
SPMR and coordinate-derived features and yields a byte-identical final
37 x 524 ChIP context table. The ChIP values supplied to the 37 x 588 multi-omic
master matrix are unchanged. Ancillary Input QC fields differ in condition and
study tables; these fields are excluded from those final products.

`single_analysis_downstream_validation.json` records the exact downstream
comparisons for this one-analysis substitution. Its historical baseline
reconstruction matched all 592 relevant SPMR values. The temporary study table
has 74 differing ancillary QC cells, while the mark-consensus and final context
tables are byte-identical. `single_analysis_treatment_spmr_forensics.json`
records the six differing treatment-signal bases and the descriptive one-read
test. Neither result relaxes strict reproduction identity or extends
re-execution claims to the other 17 analyses.

The saved `clean_clone_dry_run.log` and `clean_clone_dry_run.json` demonstrate
22 physical runs, 18 planned analyses, eight narrow and ten broad, with no
production jobs executed during the dry-run. Planning all 18 analyses does
not mean they were all re-executed.

## Selection and execution

The selected analysis is **SRR20770297**, PRJNA865478, H3K4me3, narrow peaks,
with its frozen shared Input **SRR20770287**. Selection used compressed IP
FASTQ volume only: 1,909,331,888 bytes, the smallest of six eligible IP runs.
The Input is 1,992,978,505 bytes. `selection_candidates.tsv` records the ranking.
The fragment policy is the historical fixed 147 bp policy.

A disposable Git clone of the unchanged starting HEAD was overlaid only with
the reviewed prospective patch. `commands.txt` records clone and execution
commands; `tested_sources.tsv` binds the tested scientific code, configuration
and planning inputs to the public worktree. The job ran outside the production
repository. Both FASTQs were retrieved and validated against ENA sizes/MD5s.
The original frozen reference files/index were reused as read-only inputs,
with their content/inventory/FAI hashes verified before processing. Reference
acquisition was not rerun. All new large data were confined to disposable
scratch or external validation storage, not the public repository.

The scientific environment was installed from the selected analysis's exact
historical explicit package record: all 198 package URL/hash entries match.
Snakemake 9.25.2 orchestrated the run in a separate environment. No software
version, scientific setting, pairing or comparison tolerance was changed.

## HISTORICAL PRODUCTION

Recorded source commit: `48d80c961acfd188033db9e1cf2972222fbc48a5`.
The portable implementation uses the reviewed `290c35...` production helper;
it is not claimed to be the byte-identical source set for every historical run.
Per-analysis source identities remain in `production_2026-09-23/`.
Historical SRR20770297 has 59,186 peaks, 31,707,556 usable IP reads and
FRiP 0.2082672344724393. Its Input has 33,103,220 usable reads.

## V1.0.1 CLEAN-CLONE REPRODUCTION

Job **10413311 completed successfully**, with SLURM state `COMPLETED` and exit
code `0:0`, at 2026-09-30 16:49:20 cluster time. The saved workflow log records
all four steps complete, and `job-status.tsv` records stage `completed`, exit
status 0. The intentionally cancelled dependent comparison job 10413990 was
not restarted. The completed analysis was not rerun.

The scientific comparison is **FAIL**. Usable IP reads match exactly at
31,707,556. The Input has 33,103,222 usable reads versus historical 33,103,220:
this remains **MISMATCH**, with no tolerance or exception. Total Picard
duplicates match, but filtering removes 2,658,712 duplicates versus historical
2,658,714. The cause of the two-read difference remains unresolved.

Peak count matches exactly at **59,186**, and IP FRiP matches exactly at
**0.2082672344724393** (6,603,645 overlapping IP reads). All narrowPeak and
summit coordinates match. Full scientific peak content does not match:
1,876 narrowPeak rows and 1,222 summit rows differ in numeric values.
Input peak-overlap reads match at 551,627, but the Input overlap fraction and
IP/Input overlap ratio differ because the Input denominator differs.

FASTQ MD5/size, reference identity and hashes, fixed 147 bp fragment size,
narrow peak mode and scientific MACS3 parameters match. Decoded SPMR treatment
and control-lambda checksums differ from history; these are scientific
mismatches, not acceptable gzip packaging differences. Full decoded validity
checks and exact hashes are recorded in `comparison.tsv` and
`product_checksums.tsv`. Historical control lambda can only be compared with
its preserved lossless-compression checksum because its file was removed.
Direct treatment-track inspection confirms different values at the same
genomic interval: `NC_048604.1:92740176-92740177` (0-based) is historical
0.25231 versus reproduced 0.28384. This is not just different row segmentation.

The unchanged peak count, FRiP and coordinates therefore do not establish
scientific reproduction. Numeric peak content and both SPMR tracks differ
alongside the Input count discrepancy. This association does not prove that
the two additional Input reads caused every output difference: identical IP
read counts alone do not establish identical IP alignments or fragment content.

`input_mismatch_investigation.json` preserves the initial discrepancy.
The historical control BAM was removed during validated production cleanup on
2026-09-21, preventing a read-level comparison. Ordering-dependent duplicate
representative selection is a hypothesis, not a confirmed explanation.
No rerun, tolerance relaxation or scientific parameter change was used to seek
a match. Commit/release readiness remains blocked by the scientific mismatches.

`comparison.tsv` classifies checks as `MATCH`,
`EXPECTED_NONIDENTICAL_PACKAGING` or `MISMATCH`; `scientific_comparison.json`
summarizes the result. New execution paths, dates and manifest packaging are
expected differences only where scientific fields are checked separately.
The MACS3 XLS header contains a changed scientific control-tag count, and its
scientific rows differ; its whole-file difference is classified `MISMATCH`.

## Evidence and limits

Strict equality is required for scientific metrics and deterministic text/
peak products. Execution dates, paths, manifest packaging and gzip containers
are classified separately; decoded signal content remains a scientific check.
The historical control-lambda file was deleted; its recorded lossless-
compression checksum is available for comparison. Treatment signal remains
available for direct checksum verification.

Setup-only failures and an adapter issue with Snakemake-created empty output
directories were corrected before the first alignment. The adapter correction
removes only known empty directories and has regression tests. It changes no
scientific processing. One Input and one IP analysis were then processed.

Raw ENA XML and its archive are excluded from Git, GitHub assets and Zenodo.
Metadata replay uses only the sanitized scientific snapshot. No frozen public
product, historical production output, tag or citation metadata was changed.
Nothing was staged, committed, pushed, tagged or released.
