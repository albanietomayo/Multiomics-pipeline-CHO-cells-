import csv

configfile: "config/chipseq_metadata.yaml"

CHIP = config["chipseq_metadata"]
PLAN = CHIP["plan_dir"]
SNAPSHOT = CHIP["snapshot_dir"]

INCREMENTAL_ELIGIBILITY = (
    f"{SNAPSHOT}/eligibility/incremental_eligibility.json"
)
ANALYSIS_DIR = f"{SNAPSHOT}/analysis_plan"
ANALYSIS_PLAN = f"{ANALYSIS_DIR}/chipseq_analysis_plan.tsv"
ANALYSIS_SUMMARY = f"{ANALYSIS_DIR}/analysis_plan_summary.json"
PROCESSING_RUNS = f"{ANALYSIS_DIR}/chipseq_processing_runs.tsv"
PROCESSING_SUMMARY = f"{ANALYSIS_DIR}/processing_run_summary.json"

# Legacy curated eligibility validation remains available as an
# explicit target but is no longer the default metadata-all gate.
# ChIP eligibility integration v1
CHIP_ELIGIBILITY_SAMPLES = CHIP["samples"]
CHIP_ELIGIBILITY_CONDITIONS = f"{SNAPSHOT}/control_audit/chipseq_conditions.tsv"
CHIP_ELIGIBILITY_REPORT = f"{SNAPSHOT}/eligibility/experimental_eligibility.json"
CHIP_ELIGIBILITY_PILOT = None


def chipseq_plan(wildcards):
    return str(checkpoints.chipseq_metadata_plan.get().output.plan)


def chipseq_record_files(wildcards):
    with open(chipseq_plan(wildcards), encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        f"{SNAPSHOT}/records/{row['accession']}/{name}"
        for row in rows for name in ("record.xml", "receipt.json")
    ]


rule chipseq_metadata_all:
    input:
        INCREMENTAL_ELIGIBILITY,
        ANALYSIS_PLAN,
        ANALYSIS_SUMMARY,
        PROCESSING_RUNS,
        PROCESSING_SUMMARY,
        f"{SNAPSHOT}/xml_inventory.tsv",
        f"{SNAPSHOT}/annotations/chipseq_target_annotations.tsv",
        f"{SNAPSHOT}/annotations/chipseq_label_evidence.tsv",
        f"{SNAPSHOT}/annotations/chipseq_label_summary.json",
        f"{SNAPSHOT}/control_audit/chipseq_conditions.tsv",
        f"{SNAPSHOT}/control_audit/chipseq_control_candidates.tsv",
        f"{SNAPSHOT}/control_audit/control_candidate_summary.json"


checkpoint chipseq_metadata_plan:
    input:
        samples=CHIP["samples"],
        code="workflow/scripts/build_chipseq_metadata_plan.py"
    output:
        runs=f"{PLAN}/chipseq_runs.tsv",
        plan=f"{PLAN}/ena_request_plan.tsv",
        provenance=f"{PLAN}/request_plan_provenance.json"
    params:
        outdir=PLAN
    log:
        f"{PLAN}/logs/request_plan.log"
    shell:
        "python3 {input.code:q} --samples {input.samples:q} "
        "--outdir {params.outdir:q} > {log:q} 2>&1"


rule chipseq_fetch_xml:
    input:
        plan=chipseq_plan,
        code="workflow/scripts/fetch_chipseq_ena_xml.py"
    output:
        xml=f"{SNAPSHOT}/records/{{accession}}/record.xml",
        receipt=f"{SNAPSHOT}/records/{{accession}}/receipt.json"
    params:
        outdir=SNAPSHOT
    log:
        f"{SNAPSHOT}/logs/{{accession}}.log"
    shell:
        "python3 {input.code:q} --plan {input.plan:q} "
        "--outdir {params.outdir:q} --accessions {wildcards.accession:q} "
        "> {log:q} 2>&1"


rule chipseq_xml_inventory:
    input:
        plan=chipseq_plan,
        records=chipseq_record_files,
        code="workflow/scripts/summarize_chipseq_xml.py",
        validator="workflow/scripts/fetch_chipseq_ena_xml.py"
    output:
        f"{SNAPSHOT}/xml_inventory.tsv"
    params:
        snapshot=SNAPSHOT
    log:
        f"{SNAPSHOT}/logs/inventory.log"
    shell:
        "python3 {input.code:q} --plan {input.plan:q} "
        "--snapshot {params.snapshot:q} --output {output:q} > {log:q} 2>&1"

rule chipseq_annotate_labels:
    input:
        runs=lambda wildcards: str(
            checkpoints.chipseq_metadata_plan.get().output.runs
        ),
        inventory=f"{SNAPSHOT}/xml_inventory.tsv",
        rules="config/chipseq_label_rules.json",
        code="workflow/scripts/annotate_chipseq_labels.py",
        validator="workflow/scripts/fetch_chipseq_ena_xml.py"
    output:
        annotations=f"{SNAPSHOT}/annotations/chipseq_target_annotations.tsv",
        evidence=f"{SNAPSHOT}/annotations/chipseq_label_evidence.tsv",
        summary=f"{SNAPSHOT}/annotations/chipseq_label_summary.json"
    params:
        snapshot=SNAPSHOT,
        outdir=f"{SNAPSHOT}/annotations"
    log:
        f"{SNAPSHOT}/logs/annotations.log"
    shell:
        "python3 {input.code:q} --runs {input.runs:q} "
        "--inventory {input.inventory:q} --snapshot {params.snapshot:q} "
        "--rules {input.rules:q} --outdir {params.outdir:q} "
        "> {log:q} 2>&1"

rule chipseq_control_candidates:
    input:
        runs=lambda wildcards: str(
            checkpoints.chipseq_metadata_plan.get().output.runs
        ),
        inventory=f"{SNAPSHOT}/xml_inventory.tsv",
        annotations=f"{SNAPSHOT}/annotations/chipseq_target_annotations.tsv",
        label_summary=f"{SNAPSHOT}/annotations/chipseq_label_summary.json",
        rules="config/chipseq_condition_rules.json",
        study_evidence="config/chipseq_study_evidence.json",
        code="workflow/scripts/build_chipseq_control_candidates.py",
        validator="workflow/scripts/fetch_chipseq_ena_xml.py"
    output:
        conditions=f"{SNAPSHOT}/control_audit/chipseq_conditions.tsv",
        candidates=f"{SNAPSHOT}/control_audit/chipseq_control_candidates.tsv",
        summary=f"{SNAPSHOT}/control_audit/control_candidate_summary.json"
    params:
        snapshot=SNAPSHOT,
        outdir=f"{SNAPSHOT}/control_audit"
    log:
        f"{SNAPSHOT}/logs/control_candidates.log"
    shell:
        "python3 {input.code:q} --runs {input.runs:q} "
        "--snapshot {params.snapshot:q} --rules {input.rules:q} "
        "--study-evidence {input.study_evidence:q} "
        "--outdir {params.outdir:q} > {log:q} 2>&1"


rule chipseq_incremental_eligibility:
    input:
        runs=lambda wildcards: str(
            checkpoints.chipseq_metadata_plan.get().output.runs
        ),
        annotations=(
            f"{SNAPSHOT}/annotations/"
            "chipseq_target_annotations.tsv"
        ),
        conditions=(
            f"{SNAPSHOT}/control_audit/"
            "chipseq_conditions.tsv"
        ),
        baseline_review=(
            "config/chipseq_experimental_eligibility.json"
        ),
        baseline_conditions=(
            "snapshots/chipseq/control_validation_001/"
            "chipseq_conditions.tsv"
        ),
        study_evidence="config/chipseq_study_evidence.json",
        policy=(
            "config/chipseq_incremental_eligibility_policy.json"
        ),
        code=(
            "workflow/scripts/"
            "build_chipseq_incremental_eligibility.py"
        )
    output:
        INCREMENTAL_ELIGIBILITY
    log:
        f"{SNAPSHOT}/logs/incremental_eligibility.log"
    shell:
        "python3 {input.code:q} "
        "--runs {input.runs:q} "
        "--annotations {input.annotations:q} "
        "--conditions {input.conditions:q} "
        "--baseline-review {input.baseline_review:q} "
        "--baseline-conditions {input.baseline_conditions:q} "
        "--study-evidence {input.study_evidence:q} "
        "--policy {input.policy:q} "
        "--output {output:q} "
        "> {log:q} 2>&1"


rule chipseq_analysis_plan:
    input:
        runs=lambda wildcards: str(
            checkpoints.chipseq_metadata_plan.get().output.runs
        ),
        annotations=(
            f"{SNAPSHOT}/annotations/"
            "chipseq_target_annotations.tsv"
        ),
        conditions=(
            f"{SNAPSHOT}/control_audit/"
            "chipseq_conditions.tsv"
        ),
        candidates=(
            f"{SNAPSHOT}/control_audit/"
            "chipseq_control_candidates.tsv"
        ),
        eligibility=INCREMENTAL_ELIGIBILITY,
        policy="config/chipseq_analysis_policy.json",
        code="workflow/scripts/build_chipseq_analysis_plan.py"
    output:
        plan=ANALYSIS_PLAN,
        summary=ANALYSIS_SUMMARY
    params:
        outdir=ANALYSIS_DIR
    log:
        f"{SNAPSHOT}/logs/analysis_plan.log"
    shell:
        "python3 {input.code:q} "
        "--runs {input.runs:q} "
        "--annotations {input.annotations:q} "
        "--conditions {input.conditions:q} "
        "--candidates {input.candidates:q} "
        "--eligibility {input.eligibility:q} "
        "--policy {input.policy:q} "
        "--outdir {params.outdir:q} "
        "> {log:q} 2>&1"


rule chipseq_processing_runs:
    input:
        plan=ANALYSIS_PLAN,
        samples=CHIP["samples"],
        code="workflow/scripts/prepare_chipseq_processing_runs.py"
    output:
        runs=PROCESSING_RUNS,
        summary=PROCESSING_SUMMARY
    log:
        f"{SNAPSHOT}/logs/processing_runs.log"
    shell:
        "python3 {input.code:q} "
        "--analysis-plan {input.plan:q} "
        "--samples {input.samples:q} "
        "--output {output.runs:q} "
        "--summary {output.summary:q} "
        "> {log:q} 2>&1"


include: "chipseq_eligibility.smk"
