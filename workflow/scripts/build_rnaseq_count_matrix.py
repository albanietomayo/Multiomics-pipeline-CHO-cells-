from pathlib import Path
import csv
import sys


ANNOTATION_COLUMNS = [
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length",
]


def read_run_counts(counts_file):
    """
    Read one run-level gene-count table.
    """

    rows = {}

    with open(
        counts_file,
        "r",
        encoding="utf-8",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(
                f"Missing header in {counts_file}."
            )

        if reader.fieldnames[:6] != ANNOTATION_COLUMNS:
            raise ValueError(
                f"Unexpected annotation columns in {counts_file}: "
                f"{reader.fieldnames[:6]}"
            )

        if len(reader.fieldnames) != 7:
            raise ValueError(
                f"Expected exactly one count column in {counts_file}, "
                f"found {len(reader.fieldnames) - 6}."
            )

        run_accession = reader.fieldnames[6]

        for row in reader:
            gene_id = row["Geneid"]

            if gene_id in rows:
                raise ValueError(
                    f"Duplicate gene {gene_id} in {counts_file}."
                )

            rows[gene_id] = {
                "annotation": [
                    row[column]
                    for column in ANNOTATION_COLUMNS
                ],
                "count": int(row[run_accession]),
            }

    return run_accession, rows


def main():
    if len(sys.argv) < 3:
        raise ValueError(
            "Usage: build_rnaseq_count_matrix.py "
            "<output.tsv> <run_counts_1.tsv> "
            "[run_counts_2.tsv ...]"
        )

    output_file = Path(sys.argv[1])
    input_files = [
        Path(path)
        for path in sys.argv[2:]
    ]

    run_data = []

    for counts_file in input_files:
        run_accession, rows = read_run_counts(
            counts_file
        )

        run_data.append(
            (
                run_accession,
                rows,
            )
        )

    reference_genes = set(
        run_data[0][1].keys()
    )

    for run_accession, rows in run_data[1:]:
        genes = set(rows.keys())

        if genes != reference_genes:
            raise ValueError(
                f"Gene set differs for run {run_accession}."
            )

    first_rows = run_data[0][1]

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            ["Geneid"]
            + [
                run_accession
                for run_accession, _ in run_data
            ]
        )

        for gene_id in first_rows:
            counts = []

            for run_accession, rows in run_data:
                if (
                    rows[gene_id]["annotation"]
                    != first_rows[gene_id]["annotation"]
                ):
                    raise ValueError(
                        f"Inconsistent annotation for gene "
                        f"{gene_id} in run {run_accession}."
                    )

                counts.append(
                    rows[gene_id]["count"]
                )

            writer.writerow(
                [gene_id] + counts
            )


if __name__ == "__main__":
    main()
