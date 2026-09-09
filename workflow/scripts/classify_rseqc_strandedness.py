from pathlib import Path
import csv
import sys


PAIRED_FORWARD = "1++,1--,2+-,2-+"
PAIRED_REVERSE = "1+-,1-+,2++,2--"

SINGLE_FORWARD = "++,--"
SINGLE_REVERSE = "+-,-+"


def get_orientation_fractions(row):
    """
    Convert RSeQC orientation patterns into forward and reverse fractions.
    """

    layout = row["rseqc_layout"]

    if layout == "MIXTURE":
        return None, None

    pattern_1 = row["pattern_1"]
    pattern_2 = row["pattern_2"]

    fraction_1 = float(row["pattern_1_fraction"])
    fraction_2 = float(row["pattern_2_fraction"])

    if layout == "PAIRED":
        expected = {
            PAIRED_FORWARD: "forward",
            PAIRED_REVERSE: "reverse",
        }

    elif layout == "SINGLE":
        expected = {
            SINGLE_FORWARD: "forward",
            SINGLE_REVERSE: "reverse",
        }

    else:
        raise ValueError(
            f"Unsupported RSeQC layout: {layout}"
        )

    fractions = {}

    for pattern, fraction in [
        (pattern_1, fraction_1),
        (pattern_2, fraction_2),
    ]:
        if pattern not in expected:
            raise ValueError(
                f"Unexpected RSeQC pattern for {layout}: {pattern}"
            )

        fractions[expected[pattern]] = fraction

    if set(fractions) != {"forward", "reverse"}:
        raise ValueError(
            "Could not resolve both forward and reverse fractions."
        )

    return fractions["forward"], fractions["reverse"]


def classify_strandedness(
    forward_fraction,
    reverse_fraction,
    stranded_threshold,
    unstranded_threshold,
):
    """
    Classify RNA-seq strandedness using explicit thresholds.
    """

    if forward_fraction >= stranded_threshold:
        return "forward"

    if reverse_fraction >= stranded_threshold:
        return "reverse"

    if abs(forward_fraction - reverse_fraction) < unstranded_threshold:
        return "unstranded"

    return "undetermined"


def main():
    if len(sys.argv) != 5:
        raise ValueError(
            "Usage: classify_rseqc_strandedness.py "
            "<input.tsv> <stranded_threshold> "
            "<unstranded_threshold> <output.tsv>"
        )

    input_file = Path(sys.argv[1])
    stranded_threshold = float(sys.argv[2])
    unstranded_threshold = float(sys.argv[3])
    output_file = Path(sys.argv[4])

    rows = []

    with open(input_file, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            if row["rseqc_layout"] == "MIXTURE":
                forward_fraction = ""
                reverse_fraction = ""
                classification = "undetermined"

            else:
                (
                    forward_fraction,
                    reverse_fraction,
                ) = get_orientation_fractions(row)

                classification = classify_strandedness(
                    forward_fraction,
                    reverse_fraction,
                    stranded_threshold,
                    unstranded_threshold,
                )

            row["forward_fraction"] = forward_fraction
            row["reverse_fraction"] = reverse_fraction
            row["strandedness"] = classification

            rows.append(row)

    if not rows:
        raise ValueError(
            f"No strandedness records found in {input_file}."
        )

    fieldnames = list(rows[0].keys())

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
        )

        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
