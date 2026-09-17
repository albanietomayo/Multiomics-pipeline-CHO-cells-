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

## Runtime planning integration

The metadata workflow now connects the previously validated incremental
eligibility and analysis-planning layers to the active ENA-derived ChIP-seq
metadata products.

The runtime chain is:

1. active ChIP-seq run inventory;
2. label-derived target and role annotation;
3. reconstructed experimental conditions and control candidates;
4. incremental protocol eligibility;
5. deterministic analysis planning;
6. deduplicated processing-run selection.

Only analyses with `analysis_status = ready` contribute runs to the
processing cohort.

IP runs and shared Input controls are deduplicated before FASTQ acquisition.
A shared Input therefore appears once in the processing-run table even when
it supports several ready analyses.

The current validated catalogue produces 18 ready IP/Input analyses and
22 unique processing runs: 18 IP libraries and four shared Input controls.

This phase does not yet redirect the preprocessing entry point to the
runtime processing-run table. The validated common FASTQ/QC/preprocessing
engine remains unchanged until the runtime planning layer is independently
validated.

## Dynamic preprocessing entry point

The ChIP-seq preprocessing entry point now consumes the dynamic
`chipseq_processing_runs.tsv` generated by the runtime analysis-planning
layer instead of the original two-run pilot selection.

The operational entry point executes two sequential Snakemake stages:

1. `chipseq_metadata_all`, which refreshes the active ChIP-seq metadata,
   incremental eligibility, analysis plan and deduplicated processing cohort;
2. `chipseq_preprocessing_all`, which uses that processing cohort as the
   input run table for the existing generic FASTQ acquisition, FastQC,
   fastp, post-processing QC and MultiQC rules.

This separation preserves modularity: the ChIP-specific planning layer
determines *which* runs are processable, while `common.smk` remains the
generic engine that determines *how* their FASTQ files are acquired and
preprocessed.

Shared Input controls are represented once in the processing cohort and
therefore are downloaded and preprocessed once even when reused by several
ready IP analyses.

The production wrapper no longer assumes exactly two pilot runs. Runtime
validation compares the complete processing-run set against the generated
FASTQ manifest and both structured QC summary tables.

At this stage the preprocessing redirection is implemented but no
full 22-run FASTQ download or preprocessing production execution has yet
been launched.

## Dynamic alignment input planning

The alignment generalization is being introduced in a separate planning
layer before changing the validated Bowtie2 workflow.

`prepare_chipseq_alignment_inputs.py` builds the alignment input cohort
from three runtime artifacts:

1. the dynamic `chipseq_processing_runs.tsv` generated by the analysis
   planning layer;
2. `preprocessing_qc_by_fastq.tsv`, which provides the project-relative
   processed FASTQ path and post-processing sequence count;
3. `preprocessing_qc_by_job.tsv`, which independently records the fastp
   post-processing read count.

The processing-run table remains authoritative for IP/Input role and
cohort membership. The two preprocessing tables must exactly cover the
same run set.

For every automatically alignable run, the planner verifies the
processed FASTQ on disk and records its SHA-256 digest, byte size,
post-fastp read count and project-relative provenance path.

The current automatic alignment envelope is deliberately fail-closed:
only SINGLE-end ILLUMINA ChIP-seq runs are admitted. PAIRED layouts,
unsupported platforms, missing or additional QC rows, inconsistent read
counts, path traversal and missing processed FASTQ files abort planning
rather than being guessed or silently ignored.

The current validated dynamic cohort contains 22 unique runs: 18 IP and
4 Input libraries. The previously validated H3K4me3 pilot
SRR20770297/SRR20770287 is a normal subset of this cohort and is not
special-cased by the planner.

Phase 4A does not yet alter the validated Bowtie2 Snakefile, alignment
parameters, SLURM wrapper or reference-building logic. Those runtime
changes are deferred until the dynamic input contract has been validated
independently.

## Dynamic alignment runtime core

After independent validation of the dynamic alignment-input planner,
the ChIP-seq alignment Snakefile was generalized from the original
two-run pilot to an arbitrary non-empty set of validated run
accessions.

The validated alignment method itself was not changed. Current
SINGLE-end runs continue to use Bowtie2 with local alignment,
the `very-sensitive-local` preset and seed 0, followed by coordinate
sorting with samtools and the existing alignment QC calculations.

The alignment input plan is validated before DAG construction. Run
accessions must be unique and syntactically valid; biological roles
must be `ip` or `input`; staged destinations must follow the canonical
`inputs/<run>.fastq.gz` contract; and, when layout and platform fields
are present, only SINGLE-end ILLUMINA inputs are accepted
automatically.

FASTQ staging was also generalized so that any number of planned runs
can be copied from a preprocessing project, verified against recorded
SHA-256 digests and byte sizes, and recorded in
`inputs/verified.json`.

The alignment Snakefile no longer imports or depends on the original
pilot eligibility workflow. The validated dynamic alignment input plan
is now the execution-level cohort contract, while the reference and
Bowtie2 index remain cohort-independent and therefore reusable across
different validated ChIP-seq cohorts.

Phase 4B deliberately leaves the production alignment submitter and
SLURM wrapper unchanged. Their replacement of pilot-specific
preprocessing discovery and eligibility checks is deferred to the
submission-layer integration phase.

The dynamic alignment runtime is validated with synthetic FASTQ
fixtures and a Snakemake dry-run. No Bowtie2 alignment is executed in
this phase.

## Dynamic alignment operational integration

The alignment submission layer resolves its upstream data from the
latest submitted dynamic ChIP-seq preprocessing workflow, but does not
assume that the latest submission completed successfully.

Before an alignment submission can be prepared, the operational layer
requires the preprocessing submission to contain a numeric `job_id.txt`,
a corresponding `job_<id>/job_status.tsv` with both workflow and final
exit status equal to zero, and `job_<id>/output.sha256`.

The complete preprocessing output manifest is revalidated against the
preprocessing submission `project/` directory. The dynamic processing
run table and both preprocessing QC tables are then required explicitly.

Only after those checks does the Phase 4A planner construct
`config/chipseq_alignment_inputs.json`. The generated plan therefore
binds the alignment submission to one verified preprocessing project,
its dynamic run cohort and the SHA-256-verified processed FASTQ files.

The legacy pilot-specific `support.prepare()` interface, frozen
`submission_ouZEJU` paths, `chipseq_pilot.json` dependency and pilot
eligibility preflight are no longer part of the operational alignment
path.

The alignment SLURM wrapper now validates the generated dynamic plan,
stages and SHA-256-verifies all planned FASTQs, performs a full Snakemake
dry-run and only then enters the real alignment stage.

Operational validation of Phase 4C uses an isolated synthetic successful
preprocessing submission. The production 22-run preprocessing workflow
has not yet been executed, so the real project is expected to fail
closed until a successful dynamic preprocessing submission exists.

## Dynamic duplicate marking and filtering runtime

Phase 5B generalizes the duplicate-marking and BAM-filtering runtime
from the original two-run pilot to an arbitrary non-empty validated
SINGLE-end ChIP-seq cohort.

The runtime consumes `config/chipseq_filtering_inputs.json` and no
longer imports the pilot experimental-eligibility rule or requires
exactly two accessions.

The filtering input plan is fail-closed. At workflow parsing time it
requires:

- schema version 1;
- SINGLE-end layout;
- ILLUMINA platform;
- unique ENA run accessions;
- roles restricted to `ip` and `input`;
- at least one IP and one Input;
- canonical upstream BAM paths of the form
  `outputs/<run>/raw.sorted.bam`;
- valid SHA-256 digests and positive BAM byte sizes;
- internally consistent input, mapped and nuclear MAPQ counts.

The validated scientific filtering method is unchanged. Picard
MarkDuplicates marks duplicates without removing them at the duplicate
marking stage. The subsequent filtering step retains nuclear primary
records with MAPQ >= 30 and excludes MAPQ 255, unmapped, non-primary,
QC-failed, mitochondrial and duplicate-marked records according to the
previously validated policy.

The current generalized cohort contains 22 unique runs: 18 IP and
4 shared Input libraries. This phase validates runtime expansion only;
the real 22-run alignment outputs do not yet exist and no real
duplicate marking or filtering is executed here.

The existing filtering submitter, filtering SLURM wrapper and the
pilot-specific upstream `prepare()` interface remain intentionally
unchanged during Phase 5B. Their replacement by the verified dynamic
alignment-submission contract is reserved for Phase 5C.

## Dynamic alignment-to-filtering operational integration

Phase 5C connects the generalized filtering runtime to the dynamic
alignment submission layer.

The filtering submission interface no longer consumes a frozen pilot
alignment job, archived pilot eligibility or control-validation
artifacts.

Instead it resolves
`results/chipseq/alignment/slurm/latest_submission.txt` and fails closed
unless that submission has a numeric `job_id.txt`, a corresponding
`job_<id>` directory, `job_status.tsv` with `stage=completed` and
`exit_status=0`, a successful output-publication log and an
`output.sha256` manifest.

The latest-submission pointer alone is never treated as evidence of
successful alignment.

Small alignment provenance and QC outputs are SHA-256 revalidated before
the filtering submission is created. The alignment input provenance,
per-run alignment QC and reference provenance must agree on run
accessions, IP/Input roles, SINGLE-end layout, ILLUMINA platform,
diagnostic MAPQ and mitochondrial accession.

For each run, the filtering plan records the upstream
`outputs/<run>/raw.sorted.bam` SHA-256 digest, byte size, input-read
count, mapped-read count and nuclear MAPQ>=30 count. Full BAM SHA-256
verification remains deliberately deferred to the compute-node staging
step to avoid reading all large BAMs on the login node.

The filtering `--check` mode validates and snapshots this contract but
does not call `sbatch`, create `job_id.txt` or update the filtering
`latest_submission.txt` pointer.

The filtering SLURM wrapper no longer executes pilot eligibility.
It validates the dynamic filtering plan, performs the Snakemake dry-run,
runs duplicate marking/filtering when productively submitted, publishes
the verified filtered outputs and independently rechecks the final
`output.sha256` manifest before setting `stage=completed`.

### Upstream BAM retention policy

Phase 5C does not automatically delete alignment `raw.sorted.bam`
files.

The staged filtering copies and duplicate-marked BAMs remain temporary
node-local/Snakemake intermediates. The persistent scientific products
are the verified filtered BAMs, CSI indexes, QC reports and provenance.

Alignment raw BAMs become eligible for a separate storage-cleanup step
only after a productive filtering job has:

1. reached `stage=completed` with `exit_status=0`;
2. published all filtered BAMs and CSI indexes;
3. validated each filtered BAM by a complete read;
4. generated `output.sha256`; and
5. successfully revalidated that output manifest.

Cleanup is therefore intentionally separated from the filtering job
itself. A failed filtering run never deletes its alignment upstream.

The current project still has no productive dynamic 22-run alignment
submission, so Phase 5C validation uses synthetic alignment-output
fixtures and executes no real duplicate marking or filtering.
