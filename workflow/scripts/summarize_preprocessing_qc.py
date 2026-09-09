#!/usr/bin/env python3

"""
Summarize preprocessing and PRE/POST quality-control results.

Outputs
-------
1. Per-FASTQ QC summary:
   one row per physical FASTQ file.

2. Per-fastp-job summary:
   one row per preprocessing execution
   (single-end/unpaired or paired-end).
"""

from pathlib import Path
import json
import zipfile

import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]

CONFIG_FILE = PROJECT_ROOT / "config/config.yaml"


def load_config(config_file):
    """Load pipeline configuration from YAML."""

    with open(config_file, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)

MANIFEST = (
    PROJECT_ROOT
    / "results/metadata/validation_fastq_manifest.tsv"
)

RAW_DIR = (
    PROJECT_ROOT
    / "data/raw/fastq"
)

FASTP_DIR = (
    PROJECT_ROOT
    / "results/preprocessing/fastp"
)

RAW_FASTQC_DIR = (
    PROJECT_ROOT
    / "results/qc/raw/fastqc"
)

POST_FASTQC_DIR = (
    PROJECT_ROOT
    / "results/qc/post/fastqc"
)


def build_fastq_map(manifest):
    """
    Build one row per physical FASTQ file linking:
    RAW FASTQ -> RAW FastQC -> fastp output -> POST FastQC.
    """

    rows = []

    for _, row in manifest.iterrows():

        run = row["run_accession"]
        role = row["fastq_role"]
        filename = row["filename"]

        raw_basename = filename.removesuffix(".fastq.gz")

        raw_fastq = (
            RAW_DIR
            / run
            / filename
        )

        raw_fastqc_zip = (
            RAW_FASTQC_DIR
            / run
            / f"{raw_basename}_fastqc.zip"
        )

        if role in {"R1", "R2"}:

            processed_basename = f"{run}_{role}"

            processed_fastq = (
                FASTP_DIR
                / run
                / "PAIRED"
                / f"{processed_basename}.fastq.gz"
            )

        elif role in {"SINGLE", "UNPAIRED"}:

            processed_basename = raw_basename

            processed_fastq = (
                FASTP_DIR
                / run
                / role
                / f"{processed_basename}.fastq.gz"
            )

        else:
            raise ValueError(
                f"Unsupported fastq_role={role} "
                f"for {run} / {filename}"
            )

        post_fastqc_zip = (
            POST_FASTQC_DIR
            / run
            / f"{processed_basename}_fastqc.zip"
        )

        rows.append({
            "study_accession": row["study_accession"],
            "run_accession": run,
            "omics": row["omics"],
            "fastq_role": role,
            "filename": filename,
            "raw_fastq": str(raw_fastq),
            "raw_fastqc_zip": str(raw_fastqc_zip),
            "processed_fastq": str(processed_fastq),
            "post_fastqc_zip": str(post_fastqc_zip),
        })

    return pd.DataFrame(rows)

def read_fastqc_data(zip_path):
    """
    Read fastqc_data.txt from a FastQC ZIP archive.
    """

    with zipfile.ZipFile(zip_path) as z:

        candidates = [
            name for name in z.namelist()
            if name.endswith("/fastqc_data.txt")
        ]

        if len(candidates) != 1:
            raise ValueError(
                f"Expected exactly one fastqc_data.txt in "
                f"{zip_path}, found {len(candidates)}."
            )

        return z.read(candidates[0]).decode("utf-8")


def parse_fastqc_summary(zip_path):
    """
    Extract basic statistics and module PASS/WARN/FAIL states
    from one FastQC archive.
    """

    text = read_fastqc_data(zip_path)

    result = {}

    # --------------------------------------------------------
    # Module status
    # --------------------------------------------------------

    for line in text.splitlines():

        if (
            line.startswith(">>")
            and not line.startswith(">>END_MODULE")
        ):

            parts = line.split("\t")

            if len(parts) >= 2:

                module = (
                    parts[0]
                    .replace(">>", "")
                    .strip()
                    .lower()
                    .replace(" ", "_")
                )

                result[f"status_{module}"] = parts[1]

    # --------------------------------------------------------
    # Basic Statistics
    # --------------------------------------------------------

    basic_module = None

    for module in text.split(">>END_MODULE"):

        module = module.strip()

        if ">>Basic Statistics\t" in module:
            basic_module = module
            break

    if basic_module is None:
        raise ValueError(
            f"Basic Statistics module not found in {zip_path}"
        )

    basic = {}

    for line in basic_module.splitlines()[1:]:

        if (
            not line.strip()
            or line.startswith("#")
            or line.startswith(">>")
        ):
            continue

        parts = line.split("\t")

        if len(parts) >= 2:
            basic[parts[0]] = parts[1]

    result["total_sequences"] = int(
        basic["Total Sequences"]
    )

    result["sequence_length"] = (
        basic["Sequence length"]
    )

    result["gc_pct"] = float(
        basic["%GC"]
    )

    # --------------------------------------------------------
    # Duplication
    # --------------------------------------------------------

    deduplicated_pct = None

    for line in text.splitlines():

        if line.startswith(
            "#Total Deduplicated Percentage"
        ):

            deduplicated_pct = float(
                line.split("\t")[1]
            )

            break

    if deduplicated_pct is not None:

        result["deduplicated_pct"] = round(
            deduplicated_pct,
            2
        )

        result["duplicated_pct"] = round(
            100 - deduplicated_pct,
            2
        )

    else:

        result["deduplicated_pct"] = None
        result["duplicated_pct"] = None

    return result

def build_fastp_job_map(manifest):
    """
    Build one row per fastp execution.
    """

    rows = []

    for run_accession, group in manifest.groupby(
        "run_accession",
        sort=False
    ):

        roles = set(group["fastq_role"])
        omics = group["omics"].iloc[0]

        if {"R1", "R2"}.issubset(roles):

            json_path = (
                PROJECT_ROOT
                / "results/preprocessing/fastp/reports"
                / run_accession
                / "PAIRED"
                / f"{run_accession}.PAIRED.fastp.json"
            )

            rows.append({
                "run_accession": run_accession,
                "omics": omics,
                "job_type": "PAIRED",
                "fastq_role": "R1+R2",
                "fastp_json": str(json_path),
            })

        single_rows = group[
            group["fastq_role"].isin(
                ["SINGLE", "UNPAIRED"]
            )
        ]

        for _, row in single_rows.iterrows():

            role = row["fastq_role"]

            basename = (
                row["filename"]
                .removesuffix(".fastq.gz")
            )

            json_path = (
                PROJECT_ROOT
                / "results/preprocessing/fastp/reports"
                / run_accession
                / role
                / f"{basename}.fastp.json"
            )

            rows.append({
                "run_accession": run_accession,
                "omics": row["omics"],
                "job_type": "SINGLE_END",
                "fastq_role": role,
                "fastp_json": str(json_path),
            })

    return pd.DataFrame(rows)


def parse_fastp_json(json_path):
    """
    Extract key preprocessing metrics from one fastp JSON.
    """

    with open(
        json_path,
        "r",
        encoding="utf-8"
    ) as handle:
        data = json.load(handle)

    before = data["summary"]["before_filtering"]
    after = data["summary"]["after_filtering"]

    filtering = data.get(
        "filtering_result",
        {}
    )

    adapter = data.get(
        "adapter_cutting",
        {}
    )

    reads_before = before["total_reads"]
    reads_after = after["total_reads"]

    bases_before = before["total_bases"]
    bases_after = after["total_bases"]

    trimmed_reads = (
        adapter.get("adapter_trimmed_reads")
        or 0
    )

    trimmed_bases = (
        adapter.get("adapter_trimmed_bases")
        or 0
    )

    return {
        "reads_before": reads_before,
        "reads_after": reads_after,
        "reads_removed": (
            reads_before - reads_after
        ),
        "pct_reads_removed": round(
            100
            * (reads_before - reads_after)
            / reads_before,
            4
        ),
        "bases_before": bases_before,
        "bases_after": bases_after,
        "bases_removed": (
            bases_before - bases_after
        ),
        "adapter_trimmed_reads": trimmed_reads,
        "pct_adapter_trimmed_reads": round(
            100 * trimmed_reads / reads_before,
            4
        ),
        "adapter_trimmed_bases": trimmed_bases,
        "pct_adapter_trimmed_bases": round(
            100 * trimmed_bases / bases_before,
            4
        ),
        "low_quality_reads": filtering.get(
            "low_quality_reads",
            0
        ),
        "too_many_N_reads": filtering.get(
            "too_many_N_reads",
            0
        ),
        "too_short_reads": filtering.get(
            "too_short_reads",
            0
        ),
        "too_long_reads": filtering.get(
            "too_long_reads",
            0
        ),
        "adapter_dimer_reads": filtering.get(
            "adapter_dimer_reads",
            0
        ),
    }


def to_project_relative(path_value):
    """
    Convert an absolute project path to a path relative
    to the project root for portable structured outputs.
    """

    path = Path(path_value)

    if not path.is_absolute():
        return path.as_posix()

    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def add_fastqc_transition_characterization(summary):
    """
    Characterize FastQC status transitions between RAW and POST
    without assigning biological or technical interpretation.
    """

    raw_prefix = "raw_status_"
    post_prefix = "post_status_"

    raw_modules = {
        column.removeprefix(raw_prefix)
        for column in summary.columns
        if column.startswith(raw_prefix)
    }

    post_modules = {
        column.removeprefix(post_prefix)
        for column in summary.columns
        if column.startswith(post_prefix)
    }

    modules = sorted(raw_modules & post_modules)

    transition_rows = []

    for _, row in summary.iterrows():

        resolved = []
        reduced_flag = []
        new_flag = []
        persistent_nonpass = []

        for module in modules:

            raw = row[f"raw_status_{module}"]
            post = row[f"post_status_{module}"]

            # No information in either state: no transition.
            if pd.isna(raw) and pd.isna(post):
                continue

            # Do not interpret availability changes as quality changes.
            if pd.isna(raw) or pd.isna(post):
                continue

            raw = str(raw).strip().lower()
            post = str(post).strip().lower()

            # FAIL -> PASS
            if raw == "fail" and post == "pass":
                resolved.append(module)

            # Reduction in FastQC flag severity.
            elif (
                (raw == "fail" and post == "warn")
                or
                (raw == "warn" and post == "pass")
            ):
                reduced_flag.append(module)

            # New or more severe FastQC flag.
            elif (
                (raw == "pass" and post in {"warn", "fail"})
                or
                (raw == "warn" and post == "fail")
            ):
                new_flag.append(module)

            # WARN or FAIL remains unchanged.
            elif raw == post and raw in {"warn", "fail"}:
                persistent_nonpass.append(module)

        transition_rows.append(
            {
                "n_resolved": len(resolved),
                "resolved_modules": ";".join(resolved) or "-",
                "n_reduced_flag": len(reduced_flag),
                "reduced_flag_modules": ";".join(reduced_flag) or "-",
                "n_new_flag": len(new_flag),
                "new_flag_modules": ";".join(new_flag) or "-",
                "n_persistent_nonpass": len(persistent_nonpass),
                "persistent_nonpass_modules": (
                    ";".join(persistent_nonpass) or "-"
                ),
            }
        )

    transition_df = pd.DataFrame(
        transition_rows,
        index=summary.index,
    )

    return pd.concat(
        [
            summary,
            transition_df,
        ],
        axis=1,
    )

def main():

    print("[INFO] Building preprocessing/QC summary...")

    config = load_config(CONFIG_FILE)

    manifest = pd.read_csv(
        MANIFEST,
        sep="\t",
        dtype=str
    ).fillna("")

    fastq_map = build_fastq_map(manifest)

    print("\n--- FASTQ MAP ---\n")
    print(
        fastq_map[
            [
                "run_accession",
                "omics",
                "fastq_role",
                "filename",
                "processed_fastq",
            ]
        ].to_string(index=False)
    )

    print("\n--- FILE EXISTENCE CHECK ---\n")

    checks = {
        "raw_fastq": fastq_map["raw_fastq"].map(
            lambda x: Path(x).exists()
        ),
        "raw_fastqc_zip": fastq_map["raw_fastqc_zip"].map(
            lambda x: Path(x).exists()
        ),
        "processed_fastq": fastq_map["processed_fastq"].map(
            lambda x: Path(x).exists()
        ),
        "post_fastqc_zip": fastq_map["post_fastqc_zip"].map(
            lambda x: Path(x).exists()
        ),
    }

    for name, values in checks.items():
        print(
            f"{name}: "
            f"{values.sum()}/{len(values)} present"
        )

    # --------------------------------------------------------
    # FastQC RAW vs POST summary
    # --------------------------------------------------------

    raw_rows = []
    post_rows = []

    for _, row in fastq_map.iterrows():

        raw = parse_fastqc_summary(
            row["raw_fastqc_zip"]
        )

        post = parse_fastqc_summary(
            row["post_fastqc_zip"]
        )

        raw_rows.append(raw)
        post_rows.append(post)

    raw_df = pd.DataFrame(
        raw_rows
    ).add_prefix("raw_")

    post_df = pd.DataFrame(
        post_rows
    ).add_prefix("post_")

    summary = pd.concat(
        [
            fastq_map.reset_index(drop=True),
            raw_df,
            post_df,
        ],
        axis=1,
    )

    summary = add_fastqc_transition_characterization(
        summary
    )

    print("\n--- FASTQC PRE/POST SUMMARY ---\n")

    columns = [
        "run_accession",
        "fastq_role",
        "raw_total_sequences",
        "post_total_sequences",
        "raw_sequence_length",
        "post_sequence_length",
        "raw_gc_pct",
        "post_gc_pct",
        "raw_duplicated_pct",
        "post_duplicated_pct",
        "raw_status_adapter_content",
        "post_status_adapter_content",
    ]

    print(
        summary[columns].to_string(index=False)
    )

    transition_columns = [
        "run_accession",
        "fastq_role",
        "n_resolved",
        "resolved_modules",
        "n_reduced_flag",
        "reduced_flag_modules",
        "n_new_flag",
        "new_flag_modules",
        "n_persistent_nonpass",
        "persistent_nonpass_modules",
    ]

    print(
        "\n--- FASTQC TRANSITION CHARACTERIZATION ---\n"
    )

    print(
        summary[
            transition_columns
        ].to_string(index=False)
    )

    # --------------------------------------------------------
    # fastp job-level summary
    # --------------------------------------------------------

    fastp_jobs = build_fastp_job_map(
        manifest
    )

    fastp_metrics = []

    for _, row in fastp_jobs.iterrows():
        fastp_metrics.append(
            parse_fastp_json(
                row["fastp_json"]
            )
        )

    fastp_summary = pd.concat(
        [
            fastp_jobs.reset_index(drop=True),
            pd.DataFrame(fastp_metrics),
        ],
        axis=1,
    )

    print("\n--- FASTP JOB SUMMARY ---\n")

    columns = [
        "run_accession",
        "omics",
        "job_type",
        "fastq_role",
        "reads_before",
        "reads_after",
        "pct_reads_removed",
        "pct_adapter_trimmed_reads",
        "pct_adapter_trimmed_bases",
        "low_quality_reads",
        "too_many_N_reads",
        "too_short_reads",
        "adapter_dimer_reads",
    ]

    print(
        fastp_summary[
            columns
        ].to_string(index=False)
    )


    # --------------------------------------------------------
    # Write structured outputs
    # --------------------------------------------------------

    qc_by_fastq_file = (
        PROJECT_ROOT
        / config["preprocessing"]["qc_by_fastq"]
    )

    qc_by_job_file = (
        PROJECT_ROOT
        / config["preprocessing"]["qc_by_job"]
    )

    qc_by_fastq_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    qc_by_job_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    # Keep absolute paths for runtime operations but export
    # project-relative paths for portability.

    summary_export = summary.copy()

    fastq_path_columns = [
        "raw_fastq",
        "raw_fastqc_zip",
        "processed_fastq",
        "post_fastqc_zip",
    ]

    for column in fastq_path_columns:
        summary_export[column] = (
            summary_export[column]
            .map(to_project_relative)
        )

    fastp_summary_export = fastp_summary.copy()

    fastp_summary_export["fastp_json"] = (
        fastp_summary_export["fastp_json"]
        .map(to_project_relative)
    )

    summary_export.to_csv(
        qc_by_fastq_file,
        sep="\t",
        index=False
    )

    fastp_summary_export.to_csv(
        qc_by_job_file,
        sep="\t",
        index=False
    )

    print("\n--- STRUCTURED OUTPUTS ---\n")

    print(
        f"FASTQ-level summary: {qc_by_fastq_file} "
        f"({len(summary)} rows)"
    )

    print(
        f"fastp job-level summary: {qc_by_job_file} "
        f"({len(fastp_summary)} rows)"
    )


if __name__ == "__main__":
    main()
