# ChIP-seq incremental eligibility validation 001

This snapshot validates the second ChIP-seq generalization layer.

## Baseline preservation

All 120 previously curated ChIP-seq runs retain their original
`retained_by_rules` decision. Incremental classification does not rewrite
the curated baseline.

## Incremental eligibility

A synthetic new IP was introduced using a condition already represented by
the reviewed PRJNA865478 H3K4me3 day-3 protocol.

The new IP was not added to any hard-coded accession list.

The incremental eligibility layer classified it as `retained_by_rules`
because it matched the reviewed study/role/target/condition envelope.

## Integration with real control-candidate logic

The existing `build_chipseq_control_candidates.py` candidate-selection
function was used directly.

With the original catalogue of Inputs, the synthetic new IP had exactly one
compatible control: SRR20770287. After eligibility classification, the
Phase-1 analysis planner therefore classified the synthetic IP analysis as
`ready`.

## Ambiguity test

A second synthetic Input with the same experimental condition was then added
to the candidate set.

The existing candidate-selection logic returned two compatible Inputs and
the status `multiple_candidates_needs_review`.

No automatic Input selection was made.

This demonstrates that incremental protocol eligibility and automatic
IP/Input pairing remain separate fail-closed decisions.

## Scope limitation

The snapshot validates incremental eligibility and analysis planning after
metadata/annotation/condition generation.

It does not yet demonstrate automatic downstream execution through FASTQ
download, preprocessing, alignment, duplicate handling, BAM filtering and
MACS3 peak calling. That integration is handled in the next development
phase.
