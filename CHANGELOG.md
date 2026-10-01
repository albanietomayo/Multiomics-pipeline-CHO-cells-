# Changelog

## v1.0.1 - ChIP-seq reproducibility patch (2026-10-01)

Adds the validated upstream ChIP-seq processing implementation, stage-specific
software environments, portable execution interface and production provenance
that were unintentionally omitted from the original public v1.0.0 release.

Validation from a clean clone re-executed **one of the 18 PASS analyses** from
FASTQ: SRR20770297 IP versus SRR20770287 Input, PRJNA865478, with the fixed
147 bp fragment policy. For this analysis, peak count, FRiP and peak/summit
coordinates matched historical production. Strict signal-level identity was
not obtained: the Input retained-read count differed by two reads, and small
localized numerical signal differences remained. The cause remains unresolved.

Substituting the reproduced products of this single analysis into the frozen
downstream aggregation preserved its locus-level SPMR and peak-coordinate
features, produced a byte-identical final 37 x 524 ChIP context table, and left
the ChIP values entering the 37 x 588 multi-omic master matrix unchanged.
Ancillary Input QC fields differed in the condition and study tables; these
fields are excluded from the final ChIP context and master matrix.
**The other 17 PASS analyses were not re-executed from FASTQ.**

Scientific results, benchmark labels, frozen ChIP-seq outputs and final
multi-omic products are unchanged. Metadata replay uses a minimal sanitized
scientific snapshot; raw ENA XML is excluded from Git, GitHub release assets
and Zenodo.

This is a code-availability/reproducibility patch, not a correction of
scientific results. The single-analysis strict reproduction remains FAIL;
stable downstream values do not establish full processing identity.
