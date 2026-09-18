# Project guardrails

- Treat `snapshots/chipseq/incremental_planning_validation_001/chipseq_analysis_plan.tsv` as the authoritative biological pairing layer; do not infer or rewrite IP/Input pairings downstream.
- Fail closed on missing, malformed, ambiguous, mismatched, or unsupported inputs and evidence.
- Never leak the PRJNA865478 fixed 147 bp fragment policy into another study; PRJEB9291 requires its own valid PhantomPeakQualTools estimate with no fallback.
- Use the single harmonized CriGri-PICRH-1.0 / GCF_003668045.3 reference and derive the nuclear reference span from its verified FAI.
- Do not submit or start production SLURM jobs without explicit authorization.
- Do not push automatically.
- Run relevant tests and dry-runs before production execution.
