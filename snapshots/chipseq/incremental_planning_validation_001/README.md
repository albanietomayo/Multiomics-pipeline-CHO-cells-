# ChIP-seq incremental planning validation 001

This snapshot validates the first generalization layer of the ChIP-seq
workflow.

## Scope

The planner consumes the already generated ChIP-seq run inventory,
target/role annotations, reconstructed experimental conditions,
control-candidate audit and experimental-eligibility decisions.

It does not yet implement incremental generation of experimental-
eligibility decisions for previously unseen ENA records.

## Current catalogue

The validated catalogue contains 120 ChIP-seq runs:

- 114 IP libraries;
- 6 Input libraries.

The planner produces one analysis-plan record per IP.

## Validation result

The current cohort resolves to:

- 18 `ready` analyses;
- 6 `review_required` analyses;
- 90 `blocked_no_control` analyses.

A `ready` analysis requires one unique eligible Input with an exact
condition match.

The six ambiguous PRJEB9291 analyses associated with candidate Inputs
ERR868163 and ERR868170 remain `review_required`; no automatic choice is
made.

## Pilot recovery

The previously validated pilot is independently reconstructed by the
general pairing rule:

- IP: SRR20770297;
- Input: SRR20770287;
- target: H3K4me3;
- condition: 3 days;
- pairing rule: `unique_exact_condition_input`.

No SRR accession from the pilot is encoded in the analysis-plan builder.

## Incremental interpretation

Synthetic regression testing demonstrates that a newly introduced
IP/Input pair can enter the plan automatically once the upstream
metadata, annotation, condition and eligibility layers contain the new
records.

This snapshot does not claim full ENA-to-analysis automation. Incremental
generation of experimental-eligibility decisions remains a separate
upstream requirement.
