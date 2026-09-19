# Original ATAC-seq pilot evidence

This directory contains compact original evidence copied from historical
ATAC-seq validation jobs before cleanup of large pilot outputs.

Canonical scientific results remain those from production job 10318071.

Preserved historical stages:

- 10297369: bulk CHO alignment/filtering/TSS pilot.
- 10302460: MACS3 peak-calling and FRiP pilot.
- 10304910: reproducibility/IDR pilot.
- 10308134: validated BigWig pilot.

Each copied file is recorded in source_manifest.tsv together with the original
path, byte size, source SHA-256, copied SHA-256 and verification status.

Large BAM, tagAlign, bedGraph, peak-output intermediates and BigWig pilot files
are intentionally not copied here because the final production workflow
supersedes them and they are regenerable from the preserved pipeline,
parameters and source data.

The IDR diagnostic PNG files were intentionally excluded from Git because they
are regenerable graphical outputs rather than unique provenance evidence.

No cleanup is authorized by this directory.
