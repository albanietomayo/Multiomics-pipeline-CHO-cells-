#!/usr/bin/env python3

"""
Build a normalized FASTQ manifest from config/samples.tsv.

The input table contains one row per sequencing run, whereas the
output manifest contains one row per physical FASTQ file.

The original ENA library_layout value is preserved. An additional
processing_layout field is derived from the actual FASTQ structure
available in ENA.
"""

import argparse
from pathlib import Path, PurePosixPath

import pandas as pd


def split_ena_field(value):
    """
    Split semicolon-separated ENA fields.

    ENA stores fastq_ftp, fastq_md5 and fastq_bytes as semicolon-
    separated lists when more than one FASTQ belongs to the same run.
    """
    if pd.isna(value):
        return []

    return [
        item.strip()
        for item in str(value).split(";")
        if item.strip()
    ]


def classify_fastq_structure(run_accession, library_layout, filenames):
    """
    Determine the effective FASTQ structure of a sequencing run.

    The original library_layout is never modified.

    Returns
    -------
    processing_layout : str
        SINGLE, PAIRED or UNRESOLVED.

    fastq_structure : str
        SINGLE
        UNPAIRED_FROM_PAIRED
        PAIRED
        PAIRED_PLUS_UNPAIRED
        UNRESOLVED
    """

    single_file = f"{run_accession}.fastq.gz"
    r1_file = f"{run_accession}_1.fastq.gz"
    r2_file = f"{run_accession}_2.fastq.gz"

    observed = set(filenames)

    pacbio_subreads_file = f"{run_accession}_subreads.fastq.gz"

    # PacBio Iso-Seq / SMRT subreads
    if observed == {pacbio_subreads_file}:
        return "SINGLE", "PACBIO_SUBREADS"

    # One FASTQ only
    if observed == {single_file}:

        if library_layout == "PAIRED":
            return "SINGLE", "UNPAIRED_FROM_PAIRED"

        if library_layout == "SINGLE":
            return "SINGLE", "SINGLE"

    # Two FASTQ: R1 + R2
    if observed == {r1_file, r2_file}:
        return "PAIRED", "PAIRED"

    # Three FASTQ: R1 + R2 + unpaired reads
    if observed == {single_file, r1_file, r2_file}:
        return "PAIRED", "PAIRED_PLUS_UNPAIRED"

    return "UNRESOLVED", "UNRESOLVED"


def assign_fastq_role(
    run_accession,
    processing_layout,
    fastq_structure,
    filename,
):
    """
    Assign the biological/technical role of an individual FASTQ file.
    """

    single_file = f"{run_accession}.fastq.gz"
    r1_file = f"{run_accession}_1.fastq.gz"
    r2_file = f"{run_accession}_2.fastq.gz"
    pacbio_subreads_file = f"{run_accession}_subreads.fastq.gz"

    if filename == r1_file:
        return "R1"

    if filename == r2_file:
        return "R2"

    if filename == single_file:

        if fastq_structure == "SINGLE":
            return "SINGLE"

        if fastq_structure in {
            "UNPAIRED_FROM_PAIRED",
            "PAIRED_PLUS_UNPAIRED",
        }:
            return "UNPAIRED"
    if (
           fastq_structure == "PACBIO_SUBREADS"
           and filename == pacbio_subreads_file
        ):
           return "SUBREADS"

    return "UNRESOLVED"


def build_fastq_manifest(samples):
    """
    Transform a one-row-per-run samples table into a
    one-row-per-FASTQ manifest.
    """

    manifest_rows = []

    for _, row in samples.iterrows():

        run_accession = row["run_accession"]
        library_layout = row["library_layout"]

        fastq_urls = split_ena_field(row["fastq_ftp"])
        fastq_md5s = split_ena_field(row["fastq_md5"])
        fastq_sizes = split_ena_field(row["fastq_bytes"])

        # The three ENA fields must describe the same number of files.
        if not (
            len(fastq_urls)
            == len(fastq_md5s)
            == len(fastq_sizes)
        ):
            raise ValueError(
                f"{run_accession}: inconsistent number of entries in "
                "fastq_ftp, fastq_md5 and fastq_bytes"
            )

        filenames = [
            PurePosixPath(url).name
            for url in fastq_urls
        ]

        processing_layout, fastq_structure = (
            classify_fastq_structure(
                run_accession,
                library_layout,
                filenames,
            )
        )

        for file_index, (url, md5, size, filename) in enumerate(
            zip(
                fastq_urls,
                fastq_md5s,
                fastq_sizes,
                filenames,
            ),
            start=1,
        ):

            fastq_role = assign_fastq_role(
                run_accession,
                processing_layout,
                fastq_structure,
                filename,
            )

            manifest_rows.append(
                {
                    "study_accession": row.get(
                        "study_accession", ""
                    ),
                    "sample_accession": row.get(
                        "sample_accession", ""
                    ),
                    "experiment_accession": row.get(
                        "experiment_accession", ""
                    ),
                    "run_accession": run_accession,
                    "omics": row.get("omics", ""),
                    "library_strategy": row.get(
                        "library_strategy", ""
                    ),
                    "library_layout": library_layout,
                    "processing_layout": processing_layout,
                    "fastq_structure": fastq_structure,
                    "fastq_role": fastq_role,
                    "file_index": file_index,
                    "filename": filename,
                    "fastq_ftp": url,
                    "fastq_md5": md5,
                    "fastq_bytes": int(size),
                }
            )

    return pd.DataFrame(manifest_rows)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate one-row-per-FASTQ manifest "
            "from config/samples.tsv."
        )
    )

    parser.add_argument(
        "--samples",
        required=True,
        help="Input samples.tsv file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTQ manifest TSV.",
    )

    args = parser.parse_args()

    samples = pd.read_csv(
        args.samples,
        sep="\t",
        dtype=str,
    ).fillna("")

    manifest = build_fastq_manifest(samples)

    output_path = Path(args.output)
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    manifest.to_csv(
        output_path,
        sep="\t",
        index=False,
    )

    # One row per run for summary statistics.
    run_summary = manifest[
        [
            "run_accession",
            "library_layout",
            "processing_layout",
            "fastq_structure",
        ]
    ].drop_duplicates()

    print("\n--- FASTQ MANIFEST ---")
    print(
        "Corridas:",
        run_summary["run_accession"].nunique(),
    )
    print(
        "Archivos FASTQ:",
        len(manifest),
    )

    print("\nLibrary layout original:")
    print(
        run_summary[
            "library_layout"
        ].value_counts()
    )

    print("\nProcessing layout:")
    print(
        run_summary[
            "processing_layout"
        ].value_counts()
    )

    print("\nEstructura FASTQ:")
    print(
        run_summary[
            "fastq_structure"
        ].value_counts()
    )

    unresolved = run_summary[
        run_summary["processing_layout"]
        == "UNRESOLVED"
    ]

    print(
        "\nCorridas no resueltas:",
        len(unresolved),
    )

    if not unresolved.empty:
        print(
            unresolved.to_string(
                index=False
            )
        )

    print(
        f"\nManifest escrito en: {output_path}"
    )


if __name__ == "__main__":
    main()
