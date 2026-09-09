#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


VALID_STRANDEDNESS = {
    "forward",
    "reverse",
    "unstranded",
    "undetermined",
}

VALID_MATE_PROVENANCE = {
    "R1_ONLY",
    "R2_ONLY",
    "MIXED_R1_R2",
    "UNKNOWN",
    "EMPTY",
}


def invert_strandedness(strandedness):
    """
    Convert R1-oriented strandedness to the corresponding R2 orientation,
    or vice versa.
    """

    if strandedness == "forward":
        return "reverse"

    if strandedness == "reverse":
        return "forward"

    if strandedness in {"unstranded", "undetermined"}:
        return strandedness

    raise ValueError(
        f"Unsupported strandedness for inversion: {strandedness}"
    )


def read_single_classification(path):
    with open(path, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly one strandedness record in {path}, "
            f"found {len(rows)}."
        )

    row = rows[0]

    strandedness = row.get("strandedness", "")

    if strandedness not in VALID_STRANDEDNESS:
        raise ValueError(
            f"Unsupported strandedness classification in {path}: "
            f"{strandedness}"
        )

    return row


def read_mate_provenance(path, run_accession):
    with open(path, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly one mate-provenance record in {path}, "
            f"found {len(rows)}."
        )

    row = rows[0]

    if row.get("run_accession") != run_accession:
        raise ValueError(
            f"Mate-provenance run mismatch: expected {run_accession}, "
            f"found {row.get('run_accession')}."
        )

    provenance = row.get("mate_provenance", "")

    if provenance not in VALID_MATE_PROVENANCE:
        raise ValueError(
            f"Unsupported mate provenance in {path}: {provenance}"
        )

    return row


def derive_library_evidence(
    alignment_unit,
    unit_strandedness,
    mate_provenance,
):
    """
    Convert unit-level RSeQC strandedness into evidence about the
    original library strandedness.
    """

    if alignment_unit in {"PAIRED", "SINGLE"}:
        return unit_strandedness, "direct_unit_evidence"

    if alignment_unit not in {
        "UNPAIRED",
        "UNPAIRED_R1",
        "UNPAIRED_R2",
    }:
        raise ValueError(
            f"Unsupported alignment unit: {alignment_unit}"
        )

    if unit_strandedness == "undetermined":
        return "undetermined", "unit_not_informative"

    if mate_provenance == "R1_ONLY":
        return unit_strandedness, "r1_same_orientation"

    if mate_provenance == "R2_ONLY":
        return (
            invert_strandedness(unit_strandedness),
            "r2_orientation_inverted",
        )

    if mate_provenance == "MIXED_R1_R2":
        return (
            "undetermined",
            "mixed_mates_not_used_for_library_inference",
        )

    return (
        "undetermined",
        "mate_provenance_not_resolved",
    )


def resolve_library_strandedness(rows):
    informative = [
        row["library_evidence_strandedness"]
        for row in rows
        if row["library_evidence_strandedness"] != "undetermined"
    ]

    unique = sorted(set(informative))

    if not unique:
        return "undetermined", "no_informative_units"

    if len(unique) == 1:
        resolved = unique[0]

        if len(informative) == len(rows):
            status = "all_informative_units_consensus"
        else:
            status = "informative_consensus_with_undetermined_units"

        return resolved, status

    return "undetermined", "conflicting_library_evidence"


def determine_counting_strandedness(
    alignment_unit,
    library_strandedness,
    mate_provenance,
):
    """
    Determine the strand mode that featureCounts should use for one
    alignment unit.
    """

    if library_strandedness == "undetermined":
        return "undetermined"

    if alignment_unit in {"PAIRED", "SINGLE"}:
        return library_strandedness

    if alignment_unit not in {
        "UNPAIRED",
        "UNPAIRED_R1",
        "UNPAIRED_R2",
    }:
        raise ValueError(
            f"Unsupported alignment unit: {alignment_unit}"
        )

    if mate_provenance == "R1_ONLY":
        return library_strandedness

    if mate_provenance == "R2_ONLY":
        return invert_strandedness(library_strandedness)

    if mate_provenance == "MIXED_R1_R2":
        if library_strandedness == "unstranded":
            return "unstranded"

        return "requires_mate_split"

    if library_strandedness == "unstranded":
        return "unstranded"

    return "undetermined"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Resolve RNA-seq library strandedness and derive "
            "featureCounts counting strandedness for each alignment unit."
        )
    )

    parser.add_argument("run_accession")
    parser.add_argument("output_file", type=Path)

    parser.add_argument(
        "--mate-provenance",
        type=Path,
        default=None,
    )

    parser.add_argument(
        "unit_classifications",
        nargs="+",
        type=Path,
    )

    args = parser.parse_args()

    rows = [
        read_single_classification(path)
        for path in args.unit_classifications
    ]

    for row in rows:
        if row.get("run_accession") != args.run_accession:
            raise ValueError(
                f"Expected run_accession {args.run_accession}, "
                f"found {row.get('run_accession')}."
            )

    alignment_units = [
        row["alignment_unit"]
        for row in rows
    ]

    if len(alignment_units) != len(set(alignment_units)):
        raise ValueError(
            f"Duplicate alignment units for {args.run_accession}: "
            f"{alignment_units}"
        )

    has_unpaired = any(
        unit in {
            "UNPAIRED",
            "UNPAIRED_R1",
            "UNPAIRED_R2",
        }
        for unit in alignment_units
    )

    mate_row = None

    if has_unpaired:
        if args.mate_provenance is None:
            raise ValueError(
                f"Run {args.run_accession} contains an UNPAIRED unit "
                "but no mate-provenance file was provided."
            )

        mate_row = read_mate_provenance(
            args.mate_provenance,
            args.run_accession,
        )

    enriched_rows = []

    for row in rows:
        alignment_unit = row["alignment_unit"]
        unit_strandedness = row["strandedness"]

        if alignment_unit == "UNPAIRED":
            mate_provenance = mate_row["mate_provenance"]

        elif alignment_unit == "UNPAIRED_R1":
            mate_provenance = "R1_ONLY"

        elif alignment_unit == "UNPAIRED_R2":
            mate_provenance = "R2_ONLY"

        else:
            mate_provenance = "NOT_APPLICABLE"

        evidence, evidence_status = derive_library_evidence(
            alignment_unit,
            unit_strandedness,
            mate_provenance,
        )

        enriched_rows.append(
            {
                "source": row,
                "alignment_unit": alignment_unit,
                "unit_strandedness": unit_strandedness,
                "mate_provenance": mate_provenance,
                "library_evidence_strandedness": evidence,
                "library_evidence_status": evidence_status,
            }
        )

    (
        library_strandedness,
        resolution_status,
    ) = resolve_library_strandedness(enriched_rows)

    n_units = len(enriched_rows)

    n_informative_units = sum(
        row["library_evidence_strandedness"] != "undetermined"
        for row in enriched_rows
    )

    output_rows = []

    for row in enriched_rows:
        alignment_unit = row["alignment_unit"]
        mate_provenance = row["mate_provenance"]

        counting_strandedness = determine_counting_strandedness(
            alignment_unit,
            library_strandedness,
            mate_provenance,
        )

        if (
            alignment_unit in {
                "UNPAIRED",
                "UNPAIRED_R1",
                "UNPAIRED_R2",
            }
            and mate_row is not None
        ):
            mate_records = mate_row.get("records", "")
            mate_r1_fraction = mate_row.get("r1_fraction", "")
            mate_r2_fraction = mate_row.get("r2_fraction", "")
        else:
            mate_records = ""
            mate_r1_fraction = ""
            mate_r2_fraction = ""

        output_rows.append(
            {
                "run_accession": args.run_accession,
                "alignment_unit": alignment_unit,
                "unit_strandedness": row["unit_strandedness"],
                "mate_provenance": mate_provenance,
                "mate_records": mate_records,
                "mate_r1_fraction": mate_r1_fraction,
                "mate_r2_fraction": mate_r2_fraction,
                "library_evidence_strandedness":
                    row["library_evidence_strandedness"],
                "library_evidence_status":
                    row["library_evidence_status"],
                "resolved_strandedness": library_strandedness,
                "library_strandedness": library_strandedness,
                "counting_strandedness": counting_strandedness,
                "resolution_status": resolution_status,
                "n_units": n_units,
                "n_informative_units": n_informative_units,
            }
        )

    args.output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = list(output_rows[0].keys())

    with open(
        args.output_file,
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
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
