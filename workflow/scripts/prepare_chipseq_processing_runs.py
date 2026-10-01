#!/usr/bin/env python3
"""Derive the unique ChIP-seq processing cohort from ready analyses.

Only analyses with analysis_status == "ready" contribute runs.

IP and Input runs are deduplicated before FASTQ acquisition. Shared Input
controls therefore appear once in the processing cohort while retaining a
count of the ready analyses that reuse them.
"""

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


RUN_RE = re.compile(r"[DES]RR[0-9]+")


def sha256(path):
    h = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        rows = list(reader)
        fields = reader.fieldnames or []

    if not fields:
        raise ValueError(
            f"Missing TSV header: {path}"
        )

    return rows, fields


def index_unique(rows, key, label):
    result = {}

    for row in rows:
        value = str(
            row.get(key, "")
        ).strip()

        if not value:
            raise ValueError(
                f"Missing {key} in {label}"
            )

        if value in result:
            raise ValueError(
                f"Duplicate {key} in {label}: {value}"
            )

        result[value] = row

    return result


def split_semicolon(value):
    return [
        item.strip()
        for item in str(value).split(";")
        if item.strip()
    ]


def fastq_structure(sample, run):
    urls = split_semicolon(
        sample.get("fastq_ftp", "")
    )

    sizes = split_semicolon(
        sample.get("fastq_bytes", "")
    )

    md5s = split_semicolon(
        sample.get("fastq_md5", "")
    )

    if not urls:
        raise ValueError(
            f"No FASTQ files declared for {run}"
        )

    if not (
        len(urls)
        == len(sizes)
        == len(md5s)
    ):
        raise ValueError(
            f"Inconsistent FASTQ metadata for {run}: "
            f"urls={len(urls)}, sizes={len(sizes)}, md5s={len(md5s)}"
        )

    total_bytes = 0

    for size in sizes:
        value = int(size)

        if value <= 0:
            raise ValueError(
                f"Invalid FASTQ size for {run}: {value}"
            )

        total_bytes += value

    for md5 in md5s:
        if not re.fullmatch(
            r"[a-fA-F0-9]{32}",
            md5,
        ):
            raise ValueError(
                f"Invalid FASTQ MD5 for {run}: {md5}"
            )

    return len(urls), total_bytes


def build_processing_rows(
    plan_rows,
    sample_rows,
):
    samples = index_unique(
        sample_rows,
        "run_accession",
        "sample catalogue",
    )

    analysis_ids = set()
    ready_ips = set()
    records = {}

    ready_count = 0

    for row in plan_rows:
        analysis_id = str(
            row.get(
                "analysis_id",
                "",
            )
        ).strip()

        if not analysis_id:
            raise ValueError(
                "Analysis plan contains an empty analysis_id"
            )

        if analysis_id in analysis_ids:
            raise ValueError(
                f"Duplicate analysis_id: {analysis_id}"
            )

        analysis_ids.add(
            analysis_id
        )

        status = str(
            row.get(
                "analysis_status",
                "",
            )
        ).strip()

        if status != "ready":
            continue

        ready_count += 1

        ip = str(
            row.get(
                "ip_run_accession",
                "",
            )
        ).strip()

        control = str(
            row.get(
                "control_run_accession",
                "",
            )
        ).strip()

        study = str(
            row.get(
                "study_accession",
                "",
            )
        ).strip()

        target = str(
            row.get(
                "declared_target",
                "",
            )
        ).strip()

        layout = str(
            row.get(
                "library_layout",
                "",
            )
        ).strip()

        platform = str(
            row.get(
                "instrument_platform",
                "",
            )
        ).strip()

        if not ip or not control:
            raise ValueError(
                f"Ready analysis lacks IP/Input: {analysis_id}"
            )

        if ip == control:
            raise ValueError(
                f"IP and Input must be distinct: {analysis_id}"
            )

        if not RUN_RE.fullmatch(ip):
            raise ValueError(
                f"Unexpected IP run accession: {ip}"
            )

        if not RUN_RE.fullmatch(control):
            raise ValueError(
                f"Unexpected Input run accession: {control}"
            )

        if not study:
            raise ValueError(
                f"Missing study for ready analysis: {analysis_id}"
            )

        if not target:
            raise ValueError(
                f"Missing histone target for ready analysis: {analysis_id}"
            )

        if layout != "SINGLE":
            raise ValueError(
                f"Ready analysis has unsupported layout: "
                f"{analysis_id}/{layout}"
            )

        if platform != "ILLUMINA":
            raise ValueError(
                f"Ready analysis has unsupported platform: "
                f"{analysis_id}/{platform}"
            )

        if ip in ready_ips:
            raise ValueError(
                f"IP run occurs in more than one ready analysis: {ip}"
            )

        ready_ips.add(
            ip
        )

        for run, role in (
            (ip, "ip"),
            (control, "input"),
        ):
            if run not in samples:
                raise ValueError(
                    f"Ready run missing from sample catalogue: {run}"
                )

            sample = samples[
                run
            ]

            sample_study = str(
                sample.get(
                    "study_accession",
                    "",
                )
            ).strip()

            sample_omics = str(
                sample.get(
                    "omics",
                    "",
                )
            ).strip()

            sample_layout = str(
                sample.get(
                    "library_layout",
                    "",
                )
            ).strip()

            sample_platform = str(
                sample.get(
                    "instrument_platform",
                    "",
                )
            ).strip()

            if sample_study != study:
                raise ValueError(
                    f"Study mismatch for {run}: "
                    f"plan={study}, samples={sample_study}"
                )

            if sample_omics != "ChIP-seq":
                raise ValueError(
                    f"Unexpected omics for {run}: {sample_omics}"
                )

            if sample_layout != "SINGLE":
                raise ValueError(
                    f"Unsupported sample layout for {run}: "
                    f"{sample_layout}"
                )

            if sample_platform != "ILLUMINA":
                raise ValueError(
                    f"Unsupported sample platform for {run}: "
                    f"{sample_platform}"
                )

            fastq_files, fastq_bytes = fastq_structure(
                sample,
                run,
            )

            if run not in records:
                records[
                    run
                ] = {
                    "run_accession": run,
                    "library_role": role,
                    "study_accession": study,
                    "declared_target": (
                        target
                        if role == "ip"
                        else ""
                    ),
                    "library_layout": sample_layout,
                    "instrument_platform": sample_platform,
                    "ready_analysis_count": 0,
                    "analysis_ids": [],
                    "fastq_file_count": fastq_files,
                    "fastq_bytes": fastq_bytes,
                }

            item = records[
                run
            ]

            if item[
                "library_role"
            ] != role:
                raise ValueError(
                    f"Run assigned conflicting roles: {run}"
                )

            if item[
                "study_accession"
            ] != study:
                raise ValueError(
                    f"Run assigned conflicting studies: {run}"
                )

            if item[
                "fastq_file_count"
            ] != fastq_files:
                raise ValueError(
                    f"FASTQ file count changed within plan: {run}"
                )

            if item[
                "fastq_bytes"
            ] != fastq_bytes:
                raise ValueError(
                    f"FASTQ size changed within plan: {run}"
                )

            item[
                "ready_analysis_count"
            ] += 1

            item[
                "analysis_ids"
            ].append(
                analysis_id
            )

    if ready_count == 0:
        raise ValueError(
            "Analysis plan contains no ready analyses"
        )

    role_sets = {
        role: {
            run
            for run, item in records.items()
            if item[
                "library_role"
            ] == role
        }
        for role in (
            "ip",
            "input",
        )
    }

    if role_sets[
        "ip"
    ] & role_sets[
        "input"
    ]:
        raise ValueError(
            "A run cannot be both IP and Input"
        )

    output = []

    for run in sorted(
        records
    ):
        item = dict(
            records[
                run
            ]
        )

        item[
            "analysis_ids"
        ] = ";".join(
            sorted(
                item[
                    "analysis_ids"
                ]
            )
        )

        output.append(
            item
        )

    return output, ready_count


def write_tsv(path, rows):
    fields = [
        "run_accession",
        "library_role",
        "study_accession",
        "declared_target",
        "library_layout",
        "instrument_platform",
        "ready_analysis_count",
        "analysis_ids",
        "fastq_file_count",
        "fastq_bytes",
    ]

    path = Path(path)

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_name(
        path.name + ".tmp"
    )

    with temporary.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(
            rows
        )

    temporary.replace(
        path
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--analysis-plan",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--samples",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    plan_rows, _ = read_tsv(
        args.analysis_plan
    )

    sample_rows, _ = read_tsv(
        args.samples
    )

    rows, ready_count = build_processing_rows(
        plan_rows,
        sample_rows,
    )

    write_tsv(
        args.output,
        rows,
    )

    role_counts = Counter(
        row[
            "library_role"
        ]
        for row in rows
    )

    layout_counts = Counter(
        row[
            "library_layout"
        ]
        for row in rows
    )

    shared_inputs = {
        row[
            "run_accession"
        ]: int(
            row[
                "ready_analysis_count"
            ]
        )
        for row in rows
        if (
            row[
                "library_role"
            ] == "input"
            and int(
                row[
                    "ready_analysis_count"
                ]
            ) > 1
        )
    }

    summary = {
        "schema_version": 1,
        "scope": "chipseq_ready_processing_run_set",
        "ready_analysis_count": ready_count,
        "unique_processing_runs": len(
            rows
        ),
        "role_counts": dict(
            sorted(
                role_counts.items()
            )
        ),
        "layout_counts": dict(
            sorted(
                layout_counts.items()
            )
        ),
        "shared_input_reuse": dict(
            sorted(
                shared_inputs.items()
            )
        ),
        "total_fastq_files": sum(
            int(
                row[
                    "fastq_file_count"
                ]
            )
            for row in rows
        ),
        "total_fastq_bytes": sum(
            int(
                row[
                    "fastq_bytes"
                ]
            )
            for row in rows
        ),
        "input_sha256": {
            "analysis_plan": sha256(
                args.analysis_plan
            ),
            "samples": sha256(
                args.samples
            ),
        },
    }

    args.summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = args.summary.with_name(
        args.summary.name + ".tmp"
    )

    temporary.write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    temporary.replace(
        args.summary
    )

    print(
        "[OK] ChIP-seq processing run set built"
    )

    print(
        f"ready_analyses={ready_count}"
    )

    print(
        f"unique_processing_runs={len(rows)}"
    )

    for role, count in sorted(
        role_counts.items()
    ):
        print(
            f"{role}: {count}"
        )

    for run, count in sorted(
        shared_inputs.items()
    ):
        print(
            f"shared_input {run}: {count} analyses"
        )


if __name__ == "__main__":
    main()
