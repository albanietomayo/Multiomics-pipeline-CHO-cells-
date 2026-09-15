import csv

configfile: "config/chipseq_metadata.yaml"

CHIP = config["chipseq_metadata"]
PLAN = CHIP["plan_dir"]
SNAPSHOT = CHIP["snapshot_dir"]

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
        CHIP_ELIGIBILITY_REPORT,
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


include: "chipseq_eligibility.smk"
