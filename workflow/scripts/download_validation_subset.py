#!/usr/bin/env python3

"""
Download and validate the FASTQ files included in the local
validation subset.

Each FASTQ is downloaded using the logic implemented in
download_fastq.py. After processing all files, an audit table is
generated containing the expected and observed file size and MD5
checksum.
"""

import argparse
import hashlib
from pathlib import Path

import pandas as pd

from download_fastq import download_and_validate


def calculate_md5(path):
    """
    Calculate the MD5 checksum of a file.
    """
    md5 = hashlib.md5()

    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(chunk)

    return md5.hexdigest()


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Download and validate all FASTQ files "
            "included in the validation manifest."
        )
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="Validation FASTQ manifest TSV.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where raw FASTQ files are stored.",
    )

    parser.add_argument(
        "--report",
        required=True,
        help="Output TSV containing download validation results.",
    )

    args = parser.parse_args()

    manifest = pd.read_csv(
        args.manifest,
        sep="\t",
        dtype=str,
    ).fillna("")

    required_columns = {
        "run_accession",
        "filename",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
        "fastq_role",
    }

    missing_columns = required_columns - set(manifest.columns)

    if missing_columns:
        raise ValueError(
            "Missing required manifest columns: "
            f"{sorted(missing_columns)}"
        )

    output_dir = Path(args.output_dir)
    report_path = Path(args.report)

    report_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    results = []

    print("\n--- VALIDATION FASTQ DOWNLOAD ---")
    print("Archivos incluidos:", len(manifest))

    for index, row in manifest.iterrows():

        run = row["run_accession"]
        filename = row["filename"]

        output = (
            output_dir
            / run
            / filename
        )

        expected_bytes = int(
            row["fastq_bytes"]
        )

        expected_md5 = (
            row["fastq_md5"]
            .lower()
        )

        print(
            f"\n[{index + 1}/{len(manifest)}] "
            f"{run} / {filename}"
        )

        download_and_validate(
            url=row["fastq_ftp"],
            output=output,
            expected_bytes=expected_bytes,
            expected_md5=expected_md5,
        )

        observed_bytes = output.stat().st_size
        observed_md5 = calculate_md5(output)

        size_ok = observed_bytes == expected_bytes
        md5_ok = observed_md5.lower() == expected_md5

        results.append(
            {
                "run_accession": run,
                "filename": filename,
                "fastq_role": row["fastq_role"],
                "output_path": str(output),
                "expected_bytes": expected_bytes,
                "observed_bytes": observed_bytes,
                "size_ok": size_ok,
                "expected_md5": expected_md5,
                "observed_md5": observed_md5,
                "md5_ok": md5_ok,
            }
        )

    audit = pd.DataFrame(results)

    audit.to_csv(
        report_path,
        sep="\t",
        index=False,
    )

    n_files = len(audit)
    n_size_ok = int(audit["size_ok"].sum())
    n_md5_ok = int(audit["md5_ok"].sum())

    part_files = list(
        output_dir.rglob("*.part")
    )

    print("\n--- DOWNLOAD VALIDATION SUMMARY ---")
    print("FASTQ esperados:", n_files)
    print("Tamaño correcto:", n_size_ok)
    print("MD5 correcto:", n_md5_ok)
    print("Archivos .part:", len(part_files))

    if (
        n_size_ok != n_files
        or n_md5_ok != n_files
        or part_files
    ):
        raise RuntimeError(
            "FASTQ validation failed."
        )

    print(
        f"\n[OK] Informe escrito en: {report_path}"
    )


if __name__ == "__main__":
    main()
