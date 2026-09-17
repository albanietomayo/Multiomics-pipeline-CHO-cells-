# Incremental ChIP-seq generalization

## Purpose

The validated ChIP-seq pilot demonstrated preprocessing, alignment,
duplicate handling, filtering, MACS3 peak calling and signal QC for one
H3K4me3/Input pair.

The incremental-generalization layer separates biological analysis
selection from that pilot-specific execution configuration.

## Analysis planning policy

`workflow/scripts/build_chipseq_analysis_plan.py` consumes the existing:

- ChIP-seq run inventory;
- target and library-role annotations;
- reconstructed experimental conditions;
- Input-control candidate audit;
- curated experimental-eligibility decisions.

The planner is deliberately fail-closed.

An IP is marked `ready` only when:

1. the IP is retained by the current eligibility review;
2. its experimental condition is resolved;
3. exactly one Input candidate exists;
4. the Input has no documented issue;
5. the Input is itself retained by the eligibility review;
6. IP and Input have an exact match across the configured condition key;
7. the library layout is currently supported automatically.

Other states are represented explicitly as:

- `review_required`;
- `blocked_no_control`;
- `excluded`.

The planner therefore does not silently convert ambiguous metadata into
biological pairings.

## Incremental behaviour

The analysis plan is derived from current metadata rather than from fixed
SRR accessions.

This phase generalizes the analysis-planning layer conditional on the
upstream metadata and experimental-eligibility layers already covering
the active ChIP-seq cohort. It does not yet generalize the generation of
experimental-eligibility decisions for previously unseen ENA runs.

Accordingly, a newly discovered ChIP-seq IP/Input pair that has already
been classified by the upstream eligibility layer can enter the analysis
plan deterministically without adding its run accession to source code.

Full ENA-to-analysis incremental execution additionally requires an
incremental eligibility policy. That layer must preserve existing curated
decisions, automatically classify only cases supported by explicit rules
and evidence, and route unsupported or ambiguous new records to
`review_required`.

Ambiguous or incomplete records remain visible but are not executed
automatically.

## Pilot-specific MACS3 parameters

The source-informed `--nomodel --extsize 147` configuration validated for
PRJNA865478 is not treated as a universal ChIP-seq default.

Peak-calling parameter generalization is handled separately after
analysis-pair generalization, allowing assay/study-specific policies to
be introduced without altering the biological pairing logic.

## Incremental experimental eligibility

The second generalization layer replaces the assumption that the active
ChIP-seq catalogue must contain exactly the runs present in the original
curated eligibility review.

Previously reviewed runs retain their curated decision unchanged.

For a previously unseen run, automatic protocol eligibility is deliberately
restricted to the experimentally reviewed envelope. A new run is retained
automatically only when its study, library role, histone target and resolved
condition exactly match a combination already represented by a retained
baseline run from a study explicitly enabled for incremental classification.

Therefore, the incremental eligibility layer supports new sequencing runs
or later ENA depositions belonging to an already reviewed experimental
design without extending the biological interpretation beyond the available
evidence.

A new study, new experimental condition, new histone target, unresolved
label or documented issue is assigned `review_required`.

This distinction is deliberate: automatic catalogue refresh does not imply
automatic scientific approval of previously unseen experimental designs.

## Eligibility does not imply control pairing

Incremental protocol eligibility and IP/Input pairing are intentionally
separate decisions.

A newly discovered run may fall inside the reviewed protocol envelope and
therefore be eligible for processing while still failing automatic control
pairing. For example, if an ENA refresh introduces an additional Input
compatible with an experimental condition that already has an Input, the
control-candidate layer may convert the corresponding pairing from unique
to ambiguous.

Such cases remain `review_required`; incremental eligibility never grants
permission to select arbitrarily among multiple compatible controls.

Conversely, a new IP deposited for an already reviewed condition can reuse
an existing unique compatible Input when the control-candidate rules resolve
that association deterministically.
