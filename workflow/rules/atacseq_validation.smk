# Technical checks are workflow products, separate from biological acceptance.
import csv
from selection_gate import require as require_atac_eligibility


def atacseq_selected_runs():
    path = config["atacseq"].get("validation_runs", "config/atacseq_validation_runs.tsv")
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    runs = [row["run_accession"] for row in rows]
    if not runs or len(runs) != len(set(runs)):
        raise ValueError("ATAC validation list is empty or has duplicate runs")
    for row in rows:
        require_atac_eligibility(row["run_accession"], omics="ATAC-seq")
        if row["omics"] != "ATAC-seq" or row["fastq_structure"] not in ("SINGLE", "PAIRED"):
            raise ValueError("Unsupported ATAC validation metadata")
    return rows


def get_atacseq_verified_outputs(wildcards):
    return [config["atacseq"].get("verification_dir", "results/atacseq/qc/verification")
            + "/" + row["run_accession"] + ".json" for row in atacseq_selected_runs()]


def atacseq_existing_artifact(wildcards, kind):
    require_atac_eligibility(wildcards.run_accession, omics="ATAC-seq")
    root = config.get("atacseq_existing_job")
    if not root:
        raise ValueError("Set atacseq_existing_job to the completed job directory")
    run = wildcards.run_accession
    suffix = {"bam": f"filtered_bam/{run}/filtered.bam",
              "bai": f"filtered_bam/{run}/filtered.bam.bai",
              "profile": f"qc/{run}/tss_profile.tsv",
              "summary": f"qc/{run}/tss_enrichment_summary.tsv"}[kind]
    return str(Path(root) / "runs" / run / "artifacts/results/atacseq" / suffix)


def atacseq_existing_structure(wildcards):
    rows = [r for r in atacseq_selected_runs() if r["run_accession"] == wildcards.run_accession]
    if len(rows) != 1:
        raise ValueError("Run absent from ATAC selection list")
    return rows[0]["fastq_structure"]


def get_atacseq_existing_outputs(wildcards):
    return ["results/atacseq/qc/verification_existing/" + row["run_accession"] + ".json"
            for row in atacseq_selected_runs()]


rule atacseq_all:
    input:
        get_atacseq_verified_outputs


rule atacseq_existing_all:
    input:
        get_atacseq_existing_outputs


rule atacseq_validate_outputs:
    input:
        bam=get_atacseq_eligible_filtered_bam,
        bai=config["atacseq"]["filtered_bam_dir"] + "/{run_accession}/filtered.bam.bai",
        profile=config["atacseq"]["tss_profile"],
        summary=config["atacseq"]["tss_enrichment_summary"],
        manifest=lambda wc: checkpoints.validation_fastq_manifest.get().output.manifest,
        script="workflow/scripts/validate_atac_outputs.py",
        selection=list(SELECTION_FILES)
    output:
        report=config["atacseq"].get("verification_dir", "results/atacseq/qc/verification") + "/{run_accession}.json"
    params:
        structure=lambda wc: get_atacseq_run_structure(wc.run_accession),
        mapq=config["atacseq"]["min_mapq"],
        mt=config["atacseq"]["mitochondrial_accession"],
        window=config["atacseq"]["tss_window_bp"],
        bins=config["atacseq"]["tss_bin_size_bp"],
        background=config["atacseq"]["tss_background_bp"]
    resources:
        mem_mb=2000
    conda:
        "../envs/atacseq_qc.yaml"
    shell:
        """
        python {input.script:q} --run {wildcards.run_accession:q} --layout {params.structure:q} \
          --bam {input.bam:q} --bai {input.bai:q} --profile {input.profile:q} --summary {input.summary:q} \
          --min-mapq {params.mapq} --mitochondrial {params.mt:q} \
          --window {params.window} --bin-size {params.bins} --background {params.background} \
          --tmpdir {resources.tmpdir:q} --output {output.report:q}
        """


rule atacseq_validate_existing:
    input:
        bam=lambda wc: atacseq_existing_artifact(wc, "bam"),
        bai=lambda wc: atacseq_existing_artifact(wc, "bai"),
        profile=lambda wc: atacseq_existing_artifact(wc, "profile"),
        summary=lambda wc: atacseq_existing_artifact(wc, "summary"),
        script="workflow/scripts/validate_atac_outputs.py",
        runs=config["atacseq"].get("validation_runs", "config/atacseq_validation_runs.tsv"),
        selection=list(SELECTION_FILES)
    output:
        report="results/atacseq/qc/verification_existing/{run_accession}.json"
    params:
        structure=atacseq_existing_structure,
        mapq=config["atacseq"]["min_mapq"],
        mt=config["atacseq"]["mitochondrial_accession"],
        window=config["atacseq"]["tss_window_bp"],
        bins=config["atacseq"]["tss_bin_size_bp"],
        background=config["atacseq"]["tss_background_bp"]
    resources:
        mem_mb=2000
    conda:
        "../envs/atacseq_qc.yaml"
    shell:
        """
        python {input.script:q} --run {wildcards.run_accession:q} --layout {params.structure:q} \
          --bam {input.bam:q} --bai {input.bai:q} --profile {input.profile:q} --summary {input.summary:q} \
          --min-mapq {params.mapq} --mitochondrial {params.mt:q} \
          --window {params.window} --bin-size {params.bins} --background {params.background} \
          --tmpdir {resources.tmpdir:q} --output {output.report:q}
        """
