from pathlib import Path
import csv
import re
import sys


def normalize_metric_name(metric):
    """
    Convert STAR metric labels into stable snake_case column names.
    """

    metric = metric.strip()

    if metric.startswith("% of "):
        metric = metric[len("% of "):]
        suffix = "_pct"
    elif metric.endswith(" %"):
        metric = metric[:-2]
        suffix = "_pct"
    else:
        suffix = ""

    metric = metric.replace(":", " ")
    metric = metric.replace(",", " ")
    metric = metric.replace("/", " per ")

    metric = re.sub(
        r"[^A-Za-z0-9]+",
        "_",
        metric,
    )

    metric = metric.strip("_").lower()

    return metric + suffix


def normalize_metric_value(value):
    """
    Remove presentation symbols while preserving the metric value.
    """

    value = value.strip()

    if value.endswith("%"):
        value = value[:-1].strip()

    return value


def parse_star_log(log_file):
    """
    Parse metrics from a STAR Log.final.out file.
    """

    metrics = {}

    with open(log_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if "|" not in line:
                continue

            key, value = line.split("|", 1)

            key = normalize_metric_name(key)
            value = normalize_metric_value(value)

            if key:
                metrics[key] = value

    return metrics


def main():
    if len(sys.argv) != 5:
        raise ValueError(
            "Usage: summarize_star_log.py "
            "<Log.final.out> <run_accession> "
            "<alignment_unit> <output.tsv>"
        )

    log_file = Path(sys.argv[1])
    run_accession = sys.argv[2]
    alignment_unit = sys.argv[3]
    output_file = Path(sys.argv[4])

    if not log_file.exists():
        raise FileNotFoundError(
            f"STAR log not found: {log_file}"
        )

    metrics = parse_star_log(log_file)

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    row = {
        "run_accession": run_accession,
        "alignment_unit": alignment_unit,
        **metrics,
    }

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
