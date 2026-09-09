from pathlib import Path
import csv
import re
import sys


def parse_rseqc_report(report_file):
    """
    Parse the output produced by RSeQC infer_experiment.py.

    PairEnd and SingleEnd outputs contain orientation fractions.
    RSeQC may also report 'Unknown data type: Mixture'. This is
    retained explicitly instead of being treated as a parser error.
    """

    layout = None
    rseqc_status = "ok"
    failed_fraction = None
    patterns = []

    with open(report_file, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if line == "This is PairEnd Data":
                layout = "PAIRED"

            elif line == "This is SingleEnd Data":
                layout = "SINGLE"

            elif line == "Unknown data type: Mixture":
                layout = "MIXTURE"
                rseqc_status = "mixture"

            elif line.startswith(
                "Fraction of reads failed to determine:"
            ):
                failed_fraction = line.split(":", 1)[1].strip()

            elif line.startswith(
                "Fraction of reads explained by"
            ):
                match = re.match(
                    r'Fraction of reads explained by "(.*?)":\s*(\S+)',
                    line,
                )

                if match:
                    patterns.append(
                        (
                            match.group(1),
                            match.group(2),
                        )
                    )

    if layout is None:
        raise ValueError(
            f"Could not determine RSeQC library layout from {report_file}."
        )

    if layout == "MIXTURE":
        return {
            "rseqc_layout": layout,
            "rseqc_status": rseqc_status,
            "failed_to_determine_fraction": "",
            "pattern_1": "",
            "pattern_1_fraction": "",
            "pattern_2": "",
            "pattern_2_fraction": "",
        }

    if failed_fraction is None:
        raise ValueError(
            f"Could not find failed-to-determine fraction in {report_file}."
        )

    if len(patterns) != 2:
        raise ValueError(
            f"Expected exactly two strand-orientation patterns in "
            f"{report_file}, found {len(patterns)}."
        )

    return {
        "rseqc_layout": layout,
        "rseqc_status": rseqc_status,
        "failed_to_determine_fraction": failed_fraction,
        "pattern_1": patterns[0][0],
        "pattern_1_fraction": patterns[0][1],
        "pattern_2": patterns[1][0],
        "pattern_2_fraction": patterns[1][1],
    }


def main():
    if len(sys.argv) != 5:
        raise ValueError(
            "Usage: summarize_rseqc_strandedness.py "
            "<infer_experiment.txt> <run_accession> "
            "<alignment_unit> <output.tsv>"
        )

    report_file = Path(sys.argv[1])
    run_accession = sys.argv[2]
    alignment_unit = sys.argv[3]
    output_file = Path(sys.argv[4])

    if not report_file.exists():
        raise FileNotFoundError(
            f"RSeQC strandedness report not found: {report_file}"
        )

    metrics = parse_rseqc_report(report_file)

    row = {
        "run_accession": run_accession,
        "alignment_unit": alignment_unit,
        **metrics,
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
