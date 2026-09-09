from pathlib import Path
import csv
import sys


def read_tsv(tsv_file):
    """
    Read the single data row from one STAR metrics TSV.
    """

    with open(tsv_file, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly one data row in {tsv_file}, "
            f"found {len(rows)}."
        )

    return rows[0]


def main():
    if len(sys.argv) < 3:
        raise ValueError(
            "Usage: combine_star_metrics.py "
            "<output.tsv> <star_metrics_1.tsv> [star_metrics_2.tsv ...]"
        )

    output_file = Path(sys.argv[1])
    input_files = [Path(path) for path in sys.argv[2:]]

    rows = [read_tsv(path) for path in input_files]

    fieldnames = []

    for row in rows:
        for field in row:
            if field not in fieldnames:
                fieldnames.append(field)

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
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(row)


if __name__ == "__main__":
    main()
