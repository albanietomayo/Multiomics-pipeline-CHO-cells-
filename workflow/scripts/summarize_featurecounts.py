from pathlib import Path
import csv
import re
import sys


def normalize_status(status):
    """
    Convert featureCounts status labels to stable snake_case names.
    """

    status = status.strip()

    status = re.sub(
        r"[^A-Za-z0-9]+",
        "_",
        status,
    )

    return status.strip("_").lower()


def parse_featurecounts_summary(summary_file):
    """
    Parse a featureCounts .summary file.
    """

    counts = {}

    with open(
        summary_file,
        "r",
        encoding="utf-8",
    ) as handle:

        reader = csv.reader(
            handle,
            delimiter="\t",
        )

        header = next(reader)

        if len(header) != 2:
            raise ValueError(
                f"Expected two columns in {summary_file}, "
                f"found {len(header)}."
            )

        for row in reader:
            if len(row) != 2:
                raise ValueError(
                    f"Unexpected row in {summary_file}: {row}"
                )

            status = normalize_status(row[0])
            count = int(row[1])

            counts[status] = count

    if "assigned" not in counts:
        raise ValueError(
            f"'Assigned' category not found in {summary_file}."
        )

    total = sum(counts.values())

    if total == 0:
        raise ValueError(
            f"No counting units found in {summary_file}."
        )

    assigned_pct = (
        counts["assigned"] / total * 100
    )

    return counts, total, assigned_pct


def main():
    if len(sys.argv) != 5:
        raise ValueError(
            "Usage: summarize_featurecounts.py "
            "<featurecounts.txt.summary> "
            "<run_accession> <alignment_unit> "
            "<output.tsv>"
        )

    summary_file = Path(sys.argv[1])
    run_accession = sys.argv[2]
    alignment_unit = sys.argv[3]
    output_file = Path(sys.argv[4])

    counts, total, assigned_pct = (
        parse_featurecounts_summary(summary_file)
    )

    counting_unit = (
        "fragments"
        if alignment_unit == "PAIRED"
        else "reads"
    )

    row = {
        "run_accession": run_accession,
        "alignment_unit": alignment_unit,
        "counting_unit": counting_unit,
        "total_counting_units": total,
        "assigned": counts["assigned"],
        "assigned_pct": f"{assigned_pct:.4f}",
        **{
            key: value
            for key, value in counts.items()
            if key != "assigned"
        },
    }

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
            fieldnames=row.keys(),
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
