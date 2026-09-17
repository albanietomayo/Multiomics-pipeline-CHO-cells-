# ChIP-seq H3K4me3 peak-calling validation

Pilot comparison for PRJNA865478, day-3 H3K4me3 ChIP-seq.

Treatment/IP:
- SRR20770297

Matched Input/control:
- SRR20770287

## Final configuration

The production pilot uses MACS3 3.0.4 with:

- BAM single-end mode
- q-value 0.01
- `--keep-dup all`
- `--scale-to small`
- `--nomodel --extsize 147`
- `-B --SPMR`
- CHO-specific effective genome size: 2,366,634,374 bp

The 147-bp extension is source-informed and reproduces the peak-calling
extension used for these data in the source study.

## Sensitivity analysis

An initial automatic MACS3 model estimated a fragment length of 127 bp.
MACS3 also emitted a warning because the predicted fragment length was
smaller than twice the 80-bp tag length.

The automatic and fixed configurations were therefore compared as a
technical sensitivity analysis.

AUTO-127:
- 65,689 peaks
- FRiP 0.20812781
- median peak width 340 bp

FIXED-147:
- 59,186 peaks
- FRiP 0.20826724
- median peak width 402 bp

Spatial concordance:
- 92.334% of AUTO peaks overlap FIXED
- 96.342% of FIXED peaks overlap AUTO
- genomic Jaccard = 0.920933
- 96.859% of AUTO peak bases are covered by FIXED
- 94.929% of FIXED peak bases are covered by AUTO
- median nearest-summit distance = 17 bp (AUTO to FIXED)
- median nearest-summit distance = 14 bp (FIXED to AUTO)

These results support FIXED-147 as the primary pilot configuration while
showing that the inferred H3K4me3 landscape is robust to the fragment-length
strategy.

## Limitations

This is a single IP/Input pilot pair. Replicate concordance has not been
evaluated and blacklist filtering was not performed. The comparison therefore
demonstrates technical robustness of peak calling, not full biological
reproducibility.

Large bedGraph outputs are intentionally not stored in this Git snapshot.
Peak/XLS/summit identities are preserved through SHA-256 hashes.

## Repository normalization

The copied `peak_qc.tsv` files were normalized from CRLF to LF line endings
before Git archival to ensure portable SHA-256 validation across Linux
checkouts. Numerical and textual TSV content was not otherwise modified.
MACS3 log files were preserved verbatim as primary execution evidence,
including whitespace emitted by MACS3.
