#!/usr/bin/env python3

"""
Build a fail-closed ChIP-seq alignment input plan from the dynamic
processing cohort and completed preprocessing QC tables.

Current automatic alignment scope:
- ChIP-seq
- SINGLE-end
- ILLUMINA
- one processed FASTQ per run

The processing-run table is authoritative for biological role
(IP/Input) and cohort membership. The preprocessing QC tables are
authoritative for the processed FASTQ path and post-fastp read count.
"""

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


RUN_RE = re.compile(r"[DES]RR[0-9]+")


def sha256(path):
    digest = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_tsv(path):
    path = Path(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(
                f"Missing TSV header: {path}"
            )

        return (
            list(reader),
            list(reader.fieldnames),
        )


def require_columns(
    columns,
    required,
    table_name,
):
    missing = sorted(
        set(required) - set(columns)
    )

    if missing:
        raise ValueError(
            f"{table_name} missing columns: {missing}"
        )


def project_relative(value):
    path = Path(value)

    if (
        path.is_absolute()
        or not path.parts
        or ".." in path.parts
    ):
        raise ValueError(
            f"Expected project-relative path: {value}"
        )

    return path


def unique_by_run(
    rows,
    table_name,
):
    result = {}

    for row in rows:
        run = row["run_accession"]

        if not RUN_RE.fullmatch(run):
            raise ValueError(
                f"Invalid run accession in {table_name}: {run}"
            )

        if run in result:
            raise ValueError(
                f"Duplicate run in {table_name}: {run}"
            )

        result[run] = row

    return result


def positive_int(
    value,
    label,
):
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid integer for {label}: {value}"
        ) from exc

    if result <= 0:
        raise ValueError(
            f"Expected positive integer for {label}: {value}"
        )

    return result


def build_plan(
    project_root,
    processing_runs,
    qc_by_fastq,
    qc_by_job,
):
    project_root = Path(
        project_root
    ).resolve()

    processing_path = Path(
        processing_runs
    )

    fastq_qc_path = Path(
        qc_by_fastq
    )

    job_qc_path = Path(
        qc_by_job
    )

    processing_rows, processing_columns = read_tsv(
        processing_path
    )

    fastq_rows, fastq_columns = read_tsv(
        fastq_qc_path
    )

    job_rows, job_columns = read_tsv(
        job_qc_path
    )

    require_columns(
        processing_columns,
        {
            "run_accession",
            "library_role",
            "study_accession",
            "declared_target",
            "library_layout",
            "instrument_platform",
            "ready_analysis_count",
            "analysis_ids",
            "fastq_file_count",
        },
        "processing_runs",
    )

    require_columns(
        fastq_columns,
        {
            "study_accession",
            "run_accession",
            "omics",
            "fastq_role",
            "processed_fastq",
            "post_total_sequences",
        },
        "preprocessing_qc_by_fastq",
    )

    require_columns(
        job_columns,
        {
            "run_accession",
            "omics",
            "fastq_role",
            "reads_after",
        },
        "preprocessing_qc_by_job",
    )

    if not processing_rows:
        raise ValueError(
            "Dynamic processing cohort is empty"
        )

    processing = unique_by_run(
        processing_rows,
        "processing_runs",
    )

    fastq_qc = unique_by_run(
        fastq_rows,
        "preprocessing_qc_by_fastq",
    )

    job_qc = unique_by_run(
        job_rows,
        "preprocessing_qc_by_job",
    )

    expected_runs = set(
        processing
    )

    if set(fastq_qc) != expected_runs:
        raise ValueError(
            "preprocessing_qc_by_fastq run set does not "
            "exactly match the dynamic processing cohort"
        )

    if set(job_qc) != expected_runs:
        raise ValueError(
            "preprocessing_qc_by_job run set does not "
            "exactly match the dynamic processing cohort"
        )

    plan = []

    for run in sorted(expected_runs):
        process = processing[run]
        fastq = fastq_qc[run]
        job = job_qc[run]

        role = process[
            "library_role"
        ]

        if role not in {
            "ip",
            "input",
        }:
            raise ValueError(
                f"Unsupported biological role for {run}: {role}"
            )

        if process[
            "library_layout"
        ] != "SINGLE":
            raise ValueError(
                f"Automatic alignment currently supports "
                f"SINGLE only: {run}"
            )

        if process[
            "instrument_platform"
        ] != "ILLUMINA":
            raise ValueError(
                f"Automatic alignment currently supports "
                f"ILLUMINA only: {run}"
            )

        if positive_int(
            process["fastq_file_count"],
            f"{run} fastq_file_count",
        ) != 1:
            raise ValueError(
                f"Expected one physical FASTQ for {run}"
            )

        if fastq[
            "fastq_role"
        ] != "SINGLE":
            raise ValueError(
                f"Unexpected preprocessing FASTQ role "
                f"for {run}: {fastq['fastq_role']}"
            )

        if job[
            "fastq_role"
        ] != "SINGLE":
            raise ValueError(
                f"Unexpected preprocessing job role "
                f"for {run}: {job['fastq_role']}"
            )

        if fastq[
            "omics"
        ] != "ChIP-seq":
            raise ValueError(
                f"Unexpected omics in FASTQ QC for {run}: "
                f"{fastq['omics']}"
            )

        if job[
            "omics"
        ] != "ChIP-seq":
            raise ValueError(
                f"Unexpected omics in job QC for {run}: "
                f"{job['omics']}"
            )

        if (
            fastq["study_accession"]
            != process["study_accession"]
        ):
            raise ValueError(
                f"Study mismatch for {run}"
            )

        reads_after = positive_int(
            job["reads_after"],
            f"{run} reads_after",
        )

        post_total = positive_int(
            fastq["post_total_sequences"],
            f"{run} post_total_sequences",
        )

        if reads_after != post_total:
            raise ValueError(
                f"Processed read-count mismatch for {run}: "
                f"{reads_after} != {post_total}"
            )

        relative = project_relative(
            fastq["processed_fastq"]
        )

        source = (
            project_root
            / relative
        )

        if (
            not source.is_file()
            or source.stat().st_size <= 0
        ):
            raise ValueError(
                f"Missing processed FASTQ: {source}"
            )

        ready_analysis_count = positive_int(
            process["ready_analysis_count"],
            f"{run} ready_analysis_count",
        )

        analysis_ids = [
            value
            for value in process[
                "analysis_ids"
            ].split(";")
            if value
        ]

        if (
            len(analysis_ids)
            != ready_analysis_count
        ):
            raise ValueError(
                f"ready_analysis_count/analysis_ids "
                f"mismatch for {run}"
            )

        plan.append({
            "run_accession": run,
            "role": role,
            "study_accession": process[
                "study_accession"
            ],
            "declared_target": process[
                "declared_target"
            ],
            "analysis_ids": analysis_ids,
            "ready_analysis_count": (
                ready_analysis_count
            ),
            "library_layout": "SINGLE",
            "instrument_platform": "ILLUMINA",
            "source_relative": str(relative),
            "destination": (
                f"inputs/{run}.fastq.gz"
            ),
            "reads": reads_after,
            "bytes": source.stat().st_size,
            "sha256": sha256(source),
        })

    role_counts = {
        role: sum(
            item["role"] == role
            for item in plan
        )
        for role in (
            "ip",
            "input",
        )
    }

    if role_counts["ip"] < 1:
        raise ValueError(
            "Dynamic alignment cohort contains no IP runs"
        )

    if role_counts["input"] < 1:
        raise ValueError(
            "Dynamic alignment cohort contains no Input runs"
        )

    return {
        "schema_version": 1,
        "scope": (
            "chipseq_dynamic_single_end_alignment_inputs"
        ),
        "source_root": str(
            project_root
        ),
        "runs": plan,
        "processing_runs_sha256": sha256(
            processing_path
        ),
        "reports_sha256": {
            str(fastq_qc_path): sha256(
                fastq_qc_path
            ),
            str(job_qc_path): sha256(
                job_qc_path
            ),
        },
        "run_count": len(plan),
        "role_counts": role_counts,
        "library_layout": "SINGLE",
        "instrument_platform": "ILLUMINA",
    }


def write_json(
    path,
    data,
):
    path = Path(path)

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_name(
        path.name + ".part"
    )

    temporary.write_text(
        json.dumps(
            data,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    temporary.replace(
        path
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
    )

    parser.add_argument(
        "--project-root",
        required=True,
    )

    parser.add_argument(
        "--processing-runs",
        required=True,
    )

    parser.add_argument(
        "--qc-by-fastq",
        required=True,
    )

    parser.add_argument(
        "--qc-by-job",
        required=True,
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    args = parser.parse_args()

    plan = build_plan(
        project_root=args.project_root,
        processing_runs=args.processing_runs,
        qc_by_fastq=args.qc_by_fastq,
        qc_by_job=args.qc_by_job,
    )

    write_json(
        args.output,
        plan,
    )

    print(
        f"alignment_runs={plan['run_count']}"
    )

    print(
        f"ip_runs={plan['role_counts']['ip']}"
    )

    print(
        f"input_runs={plan['role_counts']['input']}"
    )

    print(
        "DYNAMIC_ALIGNMENT_INPUT_PLAN=PASS"
    )


if __name__ == "__main__":
    main()
