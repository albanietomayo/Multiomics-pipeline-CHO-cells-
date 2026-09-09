#!/usr/bin/env python3

"""
Build a FASTQ manifest for the local validation subset.

The script selects from the complete FASTQ manifest all physical
FASTQ files belonging to the runs listed in config/validation_runs.tsv.
"""

import argparse
from pathlib import Path

import pandas as pd


def build_validation_manifest(manifest, validation_runs):
    """
    Select all FASTQ files associated with the validation runs.
    """

    if "run_accession" not in manifest.columns:
        raise ValueError(
            "The FASTQ manifest does not contain 'run_accession'."
        )

    if "run_accession" not in validation_runs.columns:
        raise ValueError(
            "The validation table does not contain 'run_accession'."
        )

    selected_runs = set(
        validation_runs["run_accession"]
        .dropna()
        .astype(str)
    )

    if not selected_runs:
        raise ValueError(
            "No validation runs were found."
        )

    available_runs = set(
        manifest["run_accession"]
        .dropna()
        .astype(str)
    )

    missing_runs = selected_runs - available_runs

    if missing_runs:
        raise ValueError(
            "Validation runs missing from FASTQ manifest: "
            f"{sorted(missing_runs)}"
        )

    subset = manifest[
        manifest["run_accession"].isin(selected_runs)
    ].copy()

    subset["fastq_bytes"] = pd.to_numeric(
        subset["fastq_bytes"],
        errors="raise",
    )

    subset = subset.sort_values(
        [
            "omics",
            "run_accession",
            "file_index",
        ]
    )

    return subset


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate the FASTQ manifest for "
            "the local validation subset."
        )
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="Complete FASTQ manifest TSV.",
    )

    parser.add_argument(
        "--validation-runs",
        required=True,
        help="TSV containing validation run accessions.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output validation FASTQ manifest.",
    )

    args = parser.parse_args()

    manifest = pd.read_csv(
        args.manifest,
        sep="\t",
        dtype=str,
    ).fillna("")

    validation_runs = pd.read_csv(
        args.validation_runs,
        sep="\t",
        dtype=str,
    ).fillna("")

    subset = build_validation_manifest(
        manifest,
        validation_runs,
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    subset.to_csv(
        output_path,
        sep="\t",
        index=False,
    )

    total_bytes = subset["fastq_bytes"].sum()

    print("\n--- VALIDATION FASTQ MANIFEST ---")
    print(
        "Corridas:",
        subset["run_accession"].nunique(),
    )
    print(
        "Archivos FASTQ:",
        len(subset),
    )
    print(
        "Volumen:",
        f"{total_bytes / 1024**3:.3f} GiB",
    )

    print("\nRoles FASTQ:")
    print(
        subset["fastq_role"].value_counts()
    )

    print(
        f"\nManifest escrito en: {output_path}"
    )


if __name__ == "__main__":
    main()
