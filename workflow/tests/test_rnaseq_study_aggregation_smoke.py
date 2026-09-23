#!/usr/bin/env python3

import csv
import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "workflow/scripts/aggregate_rnaseq_by_study.py"
)


def write_gz_matrix(path, header, rows):
    with gzip.open(
        path,
        "wt",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(header)
        writer.writerows(rows)


def read_gz_matrix(path):
    with gzip.open(
        path,
        "rt",
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(csv.reader(handle, delimiter="\t"))


with tempfile.TemporaryDirectory() as td:
    td = Path(td)

    metadata = td / "metadata.tsv"
    logcpm = td / "logcpm.tsv.gz"
    tpm = td / "tpm.tsv.gz"
    out = td / "out"

    with open(
        metadata,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(
            [
                "experiment_accession",
                "study_accession",
                "sample_accession",
            ]
        )
        writer.writerow(["E1", "S1", "A1"])
        writer.writerow(["E2", "S1", "A2"])
        writer.writerow(["E3", "S2", "A3"])
        writer.writerow(["E4", "S3", "A4"])

    header = ["Geneid", "E1", "E2", "E3", "E4"]

    write_gz_matrix(
        logcpm,
        header,
        [
            ["G1", 1, 3, 5, 7],
            ["G2", 2, 6, 8, 10],
        ],
    )

    write_gz_matrix(
        tpm,
        header,
        [
            ["G1", 100, 300, 500, 700],
            ["G2", 900, 700, 500, 300],
        ],
    )

    cmd = [
        sys.executable,
        str(SCRIPT),
        "--logcpm",
        str(logcpm),
        "--tpm",
        str(tpm),
        "--metadata",
        str(metadata),
        "--outdir",
        str(out),
    ]

    result = subprocess.run(
        cmd,
        cwd=ROOT,
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(
            "Smoke aggregation unexpectedly failed"
        )

    expected = {
        "rnaseq_tmm_logcpm_by_study_median_3.tsv.gz",
        "rnaseq_tpm_by_study_median_3.tsv.gz",
        "study_metadata.tsv",
        "study_aggregation_qc_summary.json",
        "provenance.json",
        "SHA256SUMS.txt",
    }

    observed = {
        p.name for p in out.iterdir()
        if p.is_file()
    }

    assert observed == expected

    for name in expected:
        assert (out / name).stat().st_size > 0

    log_rows = read_gz_matrix(
        out
        / "rnaseq_tmm_logcpm_by_study_median_3.tsv.gz"
    )

    tpm_rows = read_gz_matrix(
        out
        / "rnaseq_tpm_by_study_median_3.tsv.gz"
    )

    assert log_rows[0] == [
        "Geneid", "S1", "S2", "S3"
    ]

    assert tpm_rows[0] == [
        "Geneid", "S1", "S2", "S3"
    ]

    assert [float(x) for x in log_rows[1][1:]] == [
        2.0, 5.0, 7.0
    ]

    assert [float(x) for x in log_rows[2][1:]] == [
        4.0, 8.0, 10.0
    ]

    assert [float(x) for x in tpm_rows[1][1:]] == [
        200.0, 500.0, 700.0
    ]

    assert [float(x) for x in tpm_rows[2][1:]] == [
        800.0, 500.0, 300.0
    ]

    subprocess.run(
        [
            "gzip",
            "-t",
            str(
                out
                / "rnaseq_tmm_logcpm_by_study_median_3.tsv.gz"
            ),
        ],
        check=True,
    )

    subprocess.run(
        [
            "gzip",
            "-t",
            str(
                out
                / "rnaseq_tpm_by_study_median_3.tsv.gz"
            ),
        ],
        check=True,
    )

    subprocess.run(
        ["sha256sum", "--check", "SHA256SUMS.txt"],
        cwd=out,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )

    summary = json.loads(
        (
            out
            / "study_aggregation_qc_summary.json"
        ).read_text()
    )

    assert summary["input_dimensions"] == {
        "genes": 2,
        "experiments": 4,
        "studies": 3,
    }

    assert (
        summary["aggregation_statistic"]
        == "median"
    )

    assert (
        summary["post_aggregation_renormalization"]
        is False
    )

    print(
        "SMOKE PASS: 4 experiments -> "
        "3 studies using within-study gene-wise median"
    )
    print(
        "SMOKE PASS: singleton studies preserved exactly"
    )
    print(
        "SMOKE PASS: no post-aggregation TPM "
        "renormalization"
    )
    print(
        "SMOKE PASS: exact 6 outputs, gzip integrity "
        "and SHA256SUMS verified"
    )

    bad_tpm = td / "bad_tpm.tsv.gz"
    bad_out = td / "bad_out"

    write_gz_matrix(
        bad_tpm,
        header,
        [
            ["G1", -1, 300, 500, 700],
            ["G2", 900, 700, 500, 300],
        ],
    )

    bad_cmd = [
        sys.executable,
        str(SCRIPT),
        "--logcpm",
        str(logcpm),
        "--tpm",
        str(bad_tpm),
        "--metadata",
        str(metadata),
        "--outdir",
        str(bad_out),
    ]

    bad = subprocess.run(
        bad_cmd,
        cwd=ROOT,
        text=True,
        capture_output=True,
    )

    assert bad.returncode != 0
    assert not bad_out.exists()

    print(
        "SMOKE PASS: invalid negative TPM rejected "
        "with no final output publication"
    )
