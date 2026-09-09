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


def read_featurecounts(counts_file):
    """
    Read one featureCounts output and return annotation plus gene counts.
    """

    rows = {}

    with open(
        counts_file,
        "r",
        encoding="utf-8",
    ) as handle:

        reader = csv.reader(
            (
                line
                for line in handle
                if not line.startswith("#")
            ),
            delimiter="\t",
        )

        header = next(reader)

        if len(header) != 7:
            raise ValueError(
                f"Expected 7 columns in {counts_file}, "
                f"found {len(header)}."
            )

        if header[:6] != ANNOTATION_COLUMNS:
            raise ValueError(
                f"Unexpected featureCounts annotation columns "
                f"in {counts_file}: {header[:6]}"
            )

        for row in reader:
            gene_id = row[0]

            rows[gene_id] = {
                "annotation": row[:6],
                "count": int(row[6]),
            }

    return rows


def main():
    if len(sys.argv) < 4:
        raise ValueError(
            "Usage: combine_featurecounts_units.py "
            "<run_accession> <output.tsv> "
            "<featurecounts_1.txt> [featurecounts_2.txt ...]"
        )

    run_accession = sys.argv[1]
    output_file = Path(sys.argv[2])
    input_files = [
        Path(path)
        for path in sys.argv[3:]
    ]

    combined = {}

    for counts_file in input_files:
        rows = read_featurecounts(counts_file)

        for gene_id, record in rows.items():
            if gene_id not in combined:
                combined[gene_id] = {
                    "annotation": record["annotation"],
                    "count": 0,
                }

            elif (
                combined[gene_id]["annotation"]
                != record["annotation"]
            ):
                raise ValueError(
                    f"Inconsistent annotation for gene {gene_id} "
                    f"across featureCounts files."
                )

            combined[gene_id]["count"] += record["count"]

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
            ANNOTATION_COLUMNS
            + [run_accession]
        )

        for record in combined.values():
            writer.writerow(
                record["annotation"]
                + [record["count"]]
            )


if __name__ == "__main__":
    main()
