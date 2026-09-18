import sys
from pathlib import Path
sys.path.insert(0, str(Path("workflow/scripts").resolve()))
from selection_gate import decision as selection_decision, FILES as SELECTION_FILES

import pandas as pd

rule fastq_manifest:
    input:
        samples="config/samples.tsv",
        script="workflow/scripts/build_fastq_manifest.py"

    output:
        manifest=config["fastq"]["manifest"]

    conda:
        "../envs/fastq_io.yaml"

    shell:
        """
        python {input.script} \
            --samples {input.samples} \
            --output {output.manifest}
        """


checkpoint validation_fastq_manifest:
    input:
        manifest=config["fastq"]["manifest"],
        validation_runs=config["fastq"]["validation_runs"],
        script="workflow/scripts/build_validation_manifest.py",
        selection=list(SELECTION_FILES)

    output:
        manifest=config["fastq"]["validation_manifest"],
        selection=config["fastq"]["validation_manifest"] + ".selection.tsv"

    shell:
        """
        python {input.script} \
            --manifest {input.manifest} \
            --validation-runs {input.validation_runs} \
            --output {output.manifest}
        """

def get_validation_fastq_outputs(wildcards):

    manifest_path = (
        checkpoints.validation_fastq_manifest
        .get()
        .output.manifest
    )

    manifest = pd.read_csv(
        manifest_path,
        sep="\t",
        dtype=str
    ).fillna("")

    return [
        (
            f"{config['fastq']['validation_output_dir']}/"
            f"{run_accession}/{filename}"
        )
        for run_accession, filename in zip(
            manifest["run_accession"],
            manifest["filename"]
        )
    ]

rule download_validation_fastq:
    input:
        manifest=lambda wildcards: (
            checkpoints.validation_fastq_manifest
            .get()
            .output.manifest
        ),
        script="workflow/scripts/download_fastq.py"

    output:
        fastq=(
            f"{config['fastq']['validation_output_dir']}/"
            "{run_accession}/{filename}"
        )

    shell:
        """
        python {input.script} \
            --manifest {input.manifest} \
            --run-accession {wildcards.run_accession} \
            --filename {wildcards.filename} \
            --output {output.fastq}
        """


rule validation_download_report:
    input:
        manifest=lambda wildcards: (
            checkpoints.validation_fastq_manifest
            .get()
            .output.manifest
        ),
        fastqs=get_validation_fastq_outputs,
        script="workflow/scripts/download_validation_subset.py"

    output:
        report=config["fastq"]["validation_download_report"]

    params:
        output_dir=config["fastq"]["validation_output_dir"]

    shell:
        """
        python {input.script} \
            --manifest {input.manifest} \
            --output-dir {params.output_dir} \
            --report {output.report}
        """

def get_raw_fastqc_outputs(wildcards):

    manifest_path = (
        checkpoints.validation_fastq_manifest
        .get()
        .output.manifest
    )

    manifest = pd.read_csv(
        manifest_path,
        sep="\t",
        dtype=str
    ).fillna("")

    outputs = []

    for run_accession, filename in zip(
        manifest["run_accession"],
        manifest["filename"]
    ):
        basename = filename.removesuffix(".fastq.gz")

        outputs.extend(
            [
                (
                    f"{config['qc']['raw_fastqc_dir']}/"
                    f"{run_accession}/{basename}_fastqc.html"
                ),
                (
                    f"{config['qc']['raw_fastqc_dir']}/"
                    f"{run_accession}/{basename}_fastqc.zip"
                ),
            ]
        )

    return outputs

rule fastqc_raw:
    input:
        fastq=(
            f"{config['fastq']['validation_output_dir']}/"
            "{run_accession}/{basename}.fastq.gz"
        )

    output:
        html=(
            f"{config['qc']['raw_fastqc_dir']}/"
            "{run_accession}/{basename}_fastqc.html"
        ),
        zip=(
            f"{config['qc']['raw_fastqc_dir']}/"
            "{run_accession}/{basename}_fastqc.zip"
        )

    conda:
        "../envs/qc.yaml"

    threads: 1

    params:
        outdir=lambda wildcards: (
            f"{config['qc']['raw_fastqc_dir']}/"
            f"{wildcards.run_accession}"
        )

    shell:
        """
        mkdir -p {params.outdir}

        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            {input.fastq}
        """

rule multiqc_raw:
    input:
        fastqc_outputs=get_raw_fastqc_outputs

    output:
        report=(
            f"{config['qc']['raw_multiqc_dir']}/"
            "multiqc_report.html"
        ),
        data=directory(
            f"{config['qc']['raw_multiqc_dir']}/"
            "multiqc_report_data"
        )

    conda:
        "../envs/qc.yaml"

    params:
        fastqc_dir=config["qc"]["raw_fastqc_dir"],
        outdir=config["qc"]["raw_multiqc_dir"]

    shell:
        """
        mkdir -p {params.outdir}

        multiqc \
            {params.fastqc_dir} \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force
        """

# ============================================================
# FASTP PREPROCESSING
# ============================================================

def _load_validation_manifest_for_preprocessing():
    manifest_path = (
        checkpoints.validation_fastq_manifest
        .get()
        .output.manifest
    )

    return pd.read_csv(
        manifest_path,
        sep="\t",
        dtype=str
    ).fillna("")


def get_fastp_single_input(wildcards):
    """
    Resolve the raw FASTQ corresponding to a SINGLE or UNPAIRED
    file using the structured fastq_role field from the manifest.
    """

    manifest = _load_validation_manifest_for_preprocessing()

    filename = f"{wildcards.basename}.fastq.gz"

    match = manifest[
        (manifest["run_accession"] == wildcards.run_accession)
        & (manifest["fastq_role"] == wildcards.fastq_role)
        & (manifest["filename"] == filename)
    ]

    if len(match) != 1:
        raise ValueError(
            f"Expected exactly one {wildcards.fastq_role} FASTQ for "
            f"{wildcards.run_accession} / {filename}, found {len(match)}."
        )

    filename = match.iloc[0]["filename"]

    return (
        f"{config['fastq']['validation_output_dir']}/"
        f"{wildcards.run_accession}/{filename}"
    )


def get_fastp_paired_input(wildcards, role):
    """
    Resolve R1 or R2 using fastq_role rather than filename conventions.
    """

    manifest = _load_validation_manifest_for_preprocessing()

    match = manifest[
        (manifest["run_accession"] == wildcards.run_accession)
        & (manifest["fastq_role"] == role)
    ]

    if len(match) != 1:
        raise ValueError(
            f"Expected exactly one {role} FASTQ for "
            f"{wildcards.run_accession}, found {len(match)}."
        )

    filename = match.iloc[0]["filename"]

    return (
        f"{config['fastq']['validation_output_dir']}/"
        f"{wildcards.run_accession}/{filename}"
    )


def get_fastp_outputs(wildcards):
    """
    Dynamically build the expected preprocessing outputs from
    the validation FASTQ manifest.
    """

    manifest = _load_validation_manifest_for_preprocessing()

    fastp_dir = config["preprocessing"]["fastp_dir"]
    reports_dir = config["preprocessing"]["reports_dir"]

    outputs = []

    for run_accession, group in manifest.groupby(
        "run_accession",
        sort=False
    ):
        roles = set(group["fastq_role"])

        # Detect inconsistent paired structures explicitly.
        if ("R1" in roles) != ("R2" in roles):
            raise ValueError(
                f"Incomplete paired FASTQ structure for {run_accession}: "
                f"roles={sorted(roles)}"
            )

        # R1 and R2 are processed together.
        if {"R1", "R2"}.issubset(roles):
           outputs.extend([
               (
                   f"{fastp_dir}/{run_accession}/PAIRED/"
                   f"{run_accession}_R1.fastq.gz"
               ),
               (
                   f"{fastp_dir}/{run_accession}/PAIRED/"
                   f"{run_accession}_R2.fastq.gz"
               ),
               (
                   f"{reports_dir}/{run_accession}/PAIRED/"
                   f"{run_accession}.PAIRED.fastp.json"
               ),
               (
                   f"{reports_dir}/{run_accession}/PAIRED/"
                   f"{run_accession}.PAIRED.fastp.html"
               ),
            ])

        # SINGLE and UNPAIRED files are processed independently.
        single_rows = group[
            group["fastq_role"].isin(["SINGLE", "UNPAIRED"])
        ]

        for _, row in single_rows.iterrows():

            basename = row["filename"].removesuffix(".fastq.gz")
            role = row["fastq_role"]

            outputs.extend([
                (
                    f"{fastp_dir}/{run_accession}/{role}/"
                    f"{basename}.fastq.gz"
                ),
                (
                    f"{reports_dir}/{run_accession}/{role}/"
                    f"{basename}.fastp.json"
                ),
                (
                    f"{reports_dir}/{run_accession}/{role}/"
                    f"{basename}.fastp.html"
                ),
            ])

    return outputs


rule fastp_single:
    input:
        fastq=get_fastp_single_input

    output:
        fastq=(
            f"{config['preprocessing']['fastp_dir']}/"
            "{run_accession}/{fastq_role}/{basename}.fastq.gz"
        ),
        json=(
            f"{config['preprocessing']['reports_dir']}/"
            "{run_accession}/{fastq_role}/{basename}.fastp.json"
        ),
        html=(
            f"{config['preprocessing']['reports_dir']}/"
            "{run_accession}/{fastq_role}/{basename}.fastp.html"
        )

    wildcard_constraints:
        fastq_role="SINGLE|UNPAIRED"

    conda:
        "../envs/preprocessing.yaml"

    threads: 4

    log:
        (
              f"{config['preprocessing']['reports_dir']}/"
              "{run_accession}/{fastq_role}/{basename}.fastp.log"
        )

    params:
        fastq_outdir=lambda wildcards: (
            f"{config['preprocessing']['fastp_dir']}/"
            f"{wildcards.run_accession}/{wildcards.fastq_role}"
        ),
        report_outdir=lambda wildcards: (
            f"{config['preprocessing']['reports_dir']}/"
            f"{wildcards.run_accession}/{wildcards.fastq_role}"
        )

    shell:
        """
        mkdir -p {params.fastq_outdir}
        mkdir -p {params.report_outdir}

        fastp \
            --in1 {input.fastq} \
            --out1 {output.fastq} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            -Q \
            -L \
            -G > {log} 2>&1
        """


rule fastp_paired:
    input:
        r1=lambda wildcards: get_fastp_paired_input(
            wildcards, "R1"
        ),
        r2=lambda wildcards: get_fastp_paired_input(
            wildcards, "R2"
        )

    output:
       r1=(
           f"{config['preprocessing']['fastp_dir']}/"
           "{run_accession}/PAIRED/{run_accession}_R1.fastq.gz"
       ),
       r2=(
           f"{config['preprocessing']['fastp_dir']}/"
           "{run_accession}/PAIRED/{run_accession}_R2.fastq.gz"
       ),
       json=(
           f"{config['preprocessing']['reports_dir']}/"
           "{run_accession}/PAIRED/{run_accession}.PAIRED.fastp.json"
       ),
       html=(
           f"{config['preprocessing']['reports_dir']}/"
           "{run_accession}/PAIRED/{run_accession}.PAIRED.fastp.html"
       )

    conda:
        "../envs/preprocessing.yaml"

    threads: 4

    log:
        (
           f"{config['preprocessing']['reports_dir']}/"
           "{run_accession}/PAIRED/"
           "{run_accession}.PAIRED.fastp.log"
        )

    params:
        fastq_outdir=lambda wildcards: (
            f"{config['preprocessing']['fastp_dir']}/"
            f"{wildcards.run_accession}/PAIRED"
        ),
        report_outdir=lambda wildcards: (
            f"{config['preprocessing']['reports_dir']}/"
            f"{wildcards.run_accession}/PAIRED"
        )

    shell:
        """
        mkdir -p {params.fastq_outdir}
        mkdir -p {params.report_outdir}

        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --detect_adapter_for_pe \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            -Q \
            -L \
            -G > {log} 2>&1
        """

# ============================================================
# POST-PREPROCESSING QUALITY CONTROL
# ============================================================

def get_post_fastqc_input(wildcards):
    """
    Resolve a preprocessed FASTQ using the structured fastq_role
    information from the validation manifest.
    """

    manifest = _load_validation_manifest_for_preprocessing()

    run = wildcards.run_accession
    basename = wildcards.basename

    group = manifest[
        manifest["run_accession"] == run
    ]

    matches = []

    for _, row in group.iterrows():

        role = row["fastq_role"]

        if role in ["R1", "R2"]:
            processed_basename = f"{run}_{role}"

            path = (
                f"{config['preprocessing']['fastp_dir']}/"
                f"{run}/PAIRED/{processed_basename}.fastq.gz"
            )

        elif role in ["SINGLE", "UNPAIRED"]:
            processed_basename = (
                row["filename"]
                .removesuffix(".fastq.gz")
            )

            path = (
                f"{config['preprocessing']['fastp_dir']}/"
                f"{run}/{role}/{processed_basename}.fastq.gz"
            )

        else:
            continue

        if processed_basename == basename:
            matches.append(path)

    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one preprocessed FASTQ for "
            f"{run} / {basename}, found {len(matches)}."
        )

    return matches[0]


def get_post_fastqc_outputs(wildcards):
    """
    Dynamically construct FastQC POST outputs for every
    preprocessed FASTQ in the validation manifest.
    """

    manifest = _load_validation_manifest_for_preprocessing()

    post_dir = config["qc"]["post_fastqc_dir"]

    outputs = []

    for run_accession, group in manifest.groupby(
        "run_accession",
        sort=False
    ):

        for _, row in group.iterrows():

            role = row["fastq_role"]

            if role in ["R1", "R2"]:
                basename = f"{run_accession}_{role}"

            elif role in ["SINGLE", "UNPAIRED"]:
                basename = (
                    row["filename"]
                    .removesuffix(".fastq.gz")
                )

            else:
                continue

            outputs.extend([
                (
                    f"{post_dir}/{run_accession}/"
                    f"{basename}_fastqc.html"
                ),
                (
                    f"{post_dir}/{run_accession}/"
                    f"{basename}_fastqc.zip"
                ),
            ])

    return outputs


def get_post_fastqc_zip_outputs(wildcards):
    """
    Return only FastQC POST ZIP files for MultiQC aggregation.
    """

    return [
        path
        for path in get_post_fastqc_outputs(wildcards)
        if path.endswith("_fastqc.zip")
    ]


rule fastqc_post:
    input:
        fastq=get_post_fastqc_input

    output:
        html=(
            f"{config['qc']['post_fastqc_dir']}/"
            "{run_accession}/{basename}_fastqc.html"
        ),
        zip=(
            f"{config['qc']['post_fastqc_dir']}/"
            "{run_accession}/{basename}_fastqc.zip"
        )

    conda:
        "../envs/qc.yaml"

    threads: 1

    log:
        (
            f"{config['qc']['post_fastqc_dir']}/"
            "{run_accession}/{basename}.fastqc.log"
        )

    params:
        outdir=lambda wildcards: (
            f"{config['qc']['post_fastqc_dir']}/"
            f"{wildcards.run_accession}"
        )

    shell:
        """
        mkdir -p {params.outdir}

        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            {input.fastq} \
            > {log} 2>&1
        """


rule multiqc_post:
    input:
        fastqc_outputs=get_post_fastqc_zip_outputs

    output:
        report=(
            f"{config['qc']['post_multiqc_dir']}/"
            "multiqc_report.html"
        ),
        data=directory(
            f"{config['qc']['post_multiqc_dir']}/"
            "multiqc_report_data"
        )

    conda:
        "../envs/qc.yaml"

    log:
        (
            f"{config['qc']['post_multiqc_dir']}/"
            "multiqc.log"
        )

    shell:
        """
        mkdir -p $(dirname {output.report})

        multiqc \
            {input.fastqc_outputs} \
            --outdir $(dirname {output.report}) \
            --filename multiqc_report.html \
            --force \
            > {log} 2>&1
        """

# ============================================================
# Structured preprocessing and PRE/POST QC summary
# ============================================================

rule summarize_preprocessing_qc:
    input:
        raw_fastqc=get_raw_fastqc_outputs,
        fastp_outputs=get_fastp_outputs,
        post_fastqc=get_post_fastqc_outputs,
        script="workflow/scripts/summarize_preprocessing_qc.py"

    output:
        by_fastq=config["preprocessing"]["qc_by_fastq"],
        by_job=config["preprocessing"]["qc_by_job"]

    conda:
        "../envs/reporting.yaml"

    log:
        (
            f"{config['preprocessing']['reports_dir']}/"
            "summarize_preprocessing_qc.log"
        )

    shell:
        r"""
        python {input.script} > {log} 2>&1
        """
