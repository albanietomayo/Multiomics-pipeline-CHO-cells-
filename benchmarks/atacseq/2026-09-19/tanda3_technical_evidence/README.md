# ATAC-seq Tanda 3 technical evidence

This capsule preserves the final technical and developmental evidence required
before removal of historical ATAC-seq worktrees and technical validation
directories.

## Integrated historical worktree

Branch: `atacseq-integrated-20260916`.

Its branch HEAD is preserved remotely. The remaining working-tree differences
were verified with `git diff --ignore-cr-at-eol --quiet` and contained no
substantive content changes; they represented CRLF/LF line-ending differences.

## Paired-end historical worktree

Branch: `atacseq-paired-filter`.

Its branch HEAD is preserved remotely.

The uncommitted experimental state was substantive, but the complete binary Git
diff was already captured in validation job 10280253 as
`pipeline_changes.patch`. The SHA-256 of that saved patch was verified to be
identical to the current worktree diff.

The submitted paired-end validation script stored by job 10280253 was likewise
verified to be byte-identical to the current worktree launcher.

This experimental state is historical technical evidence and is not the
canonical production ATAC-seq workflow.

## Paired-end biological interpretation

SRR29929613 was used only for paired-end technical validation. It derives from
a multiome/single-cell context and is excluded from the eligible bulk CHO ATAC
atlas. Therefore its outputs demonstrate technical behavior rather than
biological validation of bulk paired-end ATAC-seq.

## Other technical directories

The small filter/integration validation, pre-integration snapshot, and
existing-output verification directories are preserved here in full.

No cleanup is performed by creation of this capsule.
