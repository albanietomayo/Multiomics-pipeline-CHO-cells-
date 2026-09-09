#!/usr/bin/env python3

from pathlib import Path
import csv
import subprocess
import sys


STRANDEDNESS_TO_FEATURECOUNTS = {
    "unstranded": "0",
    "forward": "1",
    "reverse": "2",
}


def get_counting_configuration(
    classification_file,
    run_accession,
    alignment_unit,
):
    """
    Retrieve library and counting strandedness for one RNA-seq
    alignment unit.
    """

    matches = []

    with open(
        classification_file,
        "r",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            if (
                row["run_accession"] == run_accession
                and row["alignment_unit"] == alignment_unit
            ):
                matches.append(row)

    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one resolved strandedness record for "
            f"{run_accession}/{alignment_unit}, found {len(matches)}."
        )

    row = matches[0]

    if "counting_strandedness" not in row:
        raise ValueError(
            f"counting_strandedness is missing for "
            f"{run_accession}/{alignment_unit}. "
            "The strandedness resolver output must be regenerated."
        )

    library_strandedness = row.get(
        "library_strandedness",
        "",
    )

    counting_strandedness = row[
        "counting_strandedness"
    ]

    resolution_status = row.get(
        "resolution_status",
        "",
    )

    mate_provenance = row.get(
        "mate_provenance",
        "",
    )

    if counting_strandedness == "requires_mate_split":
        raise ValueError(
            f"{run_accession}/{alignment_unit} contains mixed R1/R2 "
            "orphan reads in a stranded library. The UNPAIRED FASTQ "
            "must be split by mate before automatic featureCounts "
            "quantification."
        )

    if counting_strandedness == "undetermined":
        raise ValueError(
            f"Counting strandedness is undetermined for "
            f"{run_accession}/{alignment_unit}. "
            f"Library resolution status: {resolution_status}. "
            "Automatic featureCounts quantification was not performed."
        )

    if counting_strandedness not in STRANDEDNESS_TO_FEATURECOUNTS:
        raise ValueError(
            f"Unsupported counting strandedness for "
            f"{run_accession}/{alignment_unit}: "
            f"{counting_strandedness}"
        )

    return {
        "library_strandedness": library_strandedness,
        "counting_strandedness": counting_strandedness,
        "resolution_status": resolution_status,
        "mate_provenance": mate_provenance,
    }


def main():
    if len(sys.argv) != 9:
        raise ValueError(
            "Usage: run_featurecounts.py "
            "<resolved_strandedness.tsv> <run_accession> "
            "<alignment_unit> <annotation.gtf> <input.bam> "
            "<threads> <output.counts.txt> "
            "<output.strandedness.txt>"
        )

    classification_file = Path(sys.argv[1])
    run_accession = sys.argv[2]
    alignment_unit = sys.argv[3]
    annotation_gtf = Path(sys.argv[4])
    input_bam = Path(sys.argv[5])
    threads = sys.argv[6]
    output_counts = Path(sys.argv[7])
    output_strandedness = Path(sys.argv[8])

    configuration = get_counting_configuration(
        classification_file,
        run_accession,
        alignment_unit,
    )

    counting_strandedness = configuration[
        "counting_strandedness"
    ]

    featurecounts_strand = (
        STRANDEDNESS_TO_FEATURECOUNTS[
            counting_strandedness
        ]
    )

    output_counts.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    command = [
        "featureCounts",
        "-T",
        threads,
        "-a",
        str(annotation_gtf),
        "-t",
        "exon",
        "-g",
        "gene_id",
        "-s",
        featurecounts_strand,
    ]

    if alignment_unit == "PAIRED":
        command.extend(
            [
                "-p",
                "--countReadPairs",
            ]
        )

    elif alignment_unit not in {
        "SINGLE",
        "UNPAIRED",
        "UNPAIRED_R1",
        "UNPAIRED_R2",
    }:
        raise ValueError(
            f"Unsupported alignment unit: {alignment_unit}"
        )

    command.extend(
        [
            "-o",
            str(output_counts),
            str(input_bam),
        ]
    )

    subprocess.run(
        command,
        check=True,
    )

    with open(
        output_strandedness,
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            f"strandedness\t{counting_strandedness}\n"
            f"library_strandedness\t"
            f"{configuration['library_strandedness']}\n"
            f"counting_strandedness\t"
            f"{counting_strandedness}\n"
            f"featureCounts_s\t{featurecounts_strand}\n"
            f"mate_provenance\t"
            f"{configuration['mate_provenance']}\n"
            f"resolution_status\t"
            f"{configuration['resolution_status']}\n"
        )


if __name__ == "__main__":
    main()
