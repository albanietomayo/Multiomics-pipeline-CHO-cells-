# Integrated ATAC verification on real persisted results

Vera job 10300242 completed on 2026-09-16, exit 0:0, elapsed 00:00:53.
Source commit: 0c74840255e558e1726b90cb55f5067737db41f1.
Two atacseq_validate_existing rules and atacseq_existing_all completed.
No alignment or TSS calculation was repeated; the inputs remain those of job 10297369.

| Run | BAM records | Technical checks | TSS score |
|---|---:|---|---:|
| SRR12774931 | 7068304 | passed | 9.693333 |
| SRR12774932 | 5858817 | passed | 9.405718 |

The imported JSON reports match their persistent SHA256SUMS.txt. The four input
hashes in each report match the original BAM, index, TSS profile and summary
hashes already archived in Git. Selection hashes match the validation source commit.
The original BAMs were not transferred or read again by this archiving script.
Biological acceptance remains not_assessed. Peak calling and FRiP are pending.
The original CRLF/attribute-only differences in the worktree are preserved in
source_changes.patch; no RNA files were modified by this archival step.
SHA256SUMS.txt is the original report checksum list. evidence_SHA256SUMS.txt
covers the imported/generated evidence in this directory. Execution log is gzip
compressed without altering its decompressed contents.
