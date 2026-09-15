configfile: "config/config.yaml"

# Legacy QC summary uses these paths relative to its own project copy.
for section, key, expected in (
    ("fastq", "validation_manifest", "results/metadata/validation_fastq_manifest.tsv"),
    ("fastq", "validation_output_dir", "data/raw/fastq"),
    ("preprocessing", "fastp_dir", "results/preprocessing/fastp"),
    ("preprocessing", "reports_dir", "results/preprocessing/fastp/reports"),
    ("qc", "raw_fastqc_dir", "results/qc/raw/fastqc"),
    ("qc", "post_fastqc_dir", "results/qc/post/fastqc"),
):
    if config[section][key] != expected:
        raise ValueError(f"QC summary path mismatch: {section}/{key}")

config["fastq"]["validation_runs"] = "config/chipseq_pilot_runs.tsv"

include: "common.smk"

# ChIP eligibility integration v1
CHIP_ELIGIBILITY_SAMPLES = "config/samples.tsv"
CHIP_ELIGIBILITY_CONDITIONS = "snapshots/chipseq/control_validation_001/chipseq_conditions.tsv"
CHIP_ELIGIBILITY_REPORT = "results/chipseq/eligibility/pilot_eligibility.json"
CHIP_ELIGIBILITY_PILOT = "config/chipseq_pilot.json"

rule chipseq_pilot_selection:
    input:
        eligibility=CHIP_ELIGIBILITY_REPORT,
        pilot="config/chipseq_pilot.json",
        samples="config/samples.tsv",
        code="workflow/scripts/prepare_chipseq_pilot_runs.py"
    output:
        "config/chipseq_pilot_runs.tsv"
    shell:
        "python {input.code:q} --pilot {input.pilot:q} "
        "--samples {input.samples:q} --output {output:q}"

rule chipseq_preprocessing_all:
    input:
        CHIP_ELIGIBILITY_REPORT,
        config["fastq"]["validation_download_report"],
        config["preprocessing"]["qc_by_fastq"],
        config["preprocessing"]["qc_by_job"],
        f"{config['qc']['raw_multiqc_dir']}/multiqc_report.html",
        f"{config['qc']['post_multiqc_dir']}/multiqc_report.html"
    default_target: True

include: "chipseq_eligibility.smk"
