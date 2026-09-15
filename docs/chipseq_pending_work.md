# Pending ChIP-seq work

## Incremental eligibility rules — requested, not implemented

Generalize eligibility to new samples and studies using versioned rules,
following the approach used for RNA-seq.

Requirements:
- Preserve decisions for existing runs when relevant evidence is unchanged.
- Classify new or changed runs as retained_by_rules, excluded or review_required.
- Record the rule, reason and evidence supporting every decision.
- Keep unresolved cases outside the approved processing set.
- Do not block unchanged approved runs merely because new catalog rows appear.
- Version metadata snapshots, rules, decisions and approved catalogs.
- Keep protocol eligibility, control pairing and technical QC separate.
- Do not approve all new runs solely because their study was previously reviewed.

The current implementation validates a frozen reviewed cohort.
Incremental updates must be implemented and validated before being advertised
as supported.

## Processing

- Interpret preprocessing QC for pilot job 10280363.
- Preserve the pilot QC evidence and its interpretation.
- Implement and validate alignment, filtering, peak calling and signal QC.
