# ChIP-seq filtering validation 001

## Scope

Real-data validation of duplicate marking and BAM filtering for the
ChIP-seq pilot selected from PRJNA865478.

Pilot condition:

- histone mark: H3K4me3
- time point: day 3
- library layout: SINGLE-END
- Input: SRR20770287
- IP: SRR20770297

## HPC execution

- SLURM job: 10297738
- workflow branch: chipseq-metadata-audit
- base repository commit before filtering implementation:
  96a13edf8772bf30cc2be588862f42be8642541c

The job completed successfully with exit code 0:0.

## Filtering strategy

Picard MarkDuplicates was used to mark duplicate reads.

Filtering retained nuclear primary alignments satisfying the configured
mapping-quality threshold (MAPQ >= 30) while excluding unmapped reads,
QC-failed reads, mitochondrial alignments, unknown MAPQ values and
duplicate-marked reads.

Removal counts in `filtering_flow.tsv` are sequential and mutually
exclusive. Independent diagnostic counts in `filtering_qc.json` may
overlap and therefore must not be summed.

Blacklist filtering was not performed at this stage.

Optical duplicate detection was not performed because read-name optical
duplicate parsing was disabled.

Signal quality was intentionally not evaluated at this stage and will be
assessed downstream using peak calling and enrichment-based metrics.

## Final results

### Input — SRR20770287

- input records: 43,830,994
- sequentially removed: 10,727,770
- retained reads: 33,103,224
- retained fraction: 75.5247 %
- Picard percent duplication: 10.1687 %

### IP — SRR20770297

- input records: 41,724,399
- sequentially removed: 10,016,844
- retained reads: 31,707,555
- retained fraction: 75.9928 %
- Picard percent duplication: 12.0830 %

## Validation

For both runs:

- complete BAM traversal: PASS
- sequential filtering arithmetic: PASS
- filtering-step continuity: PASS
- final flagstat count agreement: PASS
- filtered BAM present: PASS
- CSI index present: PASS

The persistent copied outputs were additionally verified using SHA-256.

Large BAM and index files are intentionally not stored in this Git
snapshot.
