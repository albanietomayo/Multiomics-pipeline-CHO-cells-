# Eligible CHO bulk ATAC validation — 2026-09-16

Job 10297369 completed 0:0 in 00:39:44. Both runs are CHO-K1 SINGLE.

| Run | Filtered alignments | Bowtie2 alignment | Picard duplication | TSS score |
|---|---:|---:|---:|---:|
| SRR12774931 | 7068304 | 98.63% | 77.2402% | 9.693333 |
| SRR12774932 | 5858817 | 98.40% | 74.4819% | 9.405718 |

Both used 42148 of 42215 TSS and peak bin [-60,-50). The score is the maximum
10-bp bin normalized by the mean of 100-bp flanks in the ±2000-bp window.
No ENCODE or universal acceptance threshold is applied. High duplication reduces
retained depth; these technical results do not complete atlas quality assessment.
Peak calling and FRiP are pending. Single-end is not single-cell.

Persistent BAMs and indexes (not included in Git):
/cephyr/users/mayoa/Vera/atacseq_bulk_validation_results/launch_20260916T081250Z_2qbkz612/job_10297369/runs/<RUN>/artifacts/results/atacseq/filtered_bam/<RUN>/filtered.bam

Per-run SHA256SUMS.txt retains the original artifact-relative paths and hashes,
including BAM/index. It is not a checksum manifest for this flattened report folder.
The persistent verification logs record the original successful copy verification.
The separate evidence_SHA256SUMS.txt covers files included here.

Historical code is deliberately retained byte-for-byte; its hardcoded old hashes
and paths are provenance, not the current launcher. Use workflow/scripts/launch_atac_bulk.py.
The new validation rules were not part of job 10297369. Their execution on the original
BAMs remains to be run using atacseq_existing_all; local tests do not substitute for it.
