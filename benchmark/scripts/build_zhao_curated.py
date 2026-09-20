#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path


TARGET_ASSEMBLY = "GCF_003668045.3"

EXPECTED_ZHAO_OBS_ROLES = {
    "positive": 7,
    "negative": 5,
}

EXPECTED_ZHAO_CLASSES = {
    "observed_stable": 7,
    "observed_unstable": 5,
}

EXPECTED_ZHAO_LOCUS_ROLES = {
    "C12orf35": "positive",
    "HPRT": "support_only",
    "GRIK1": "support_only",
}


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

        if reader.fieldnames is None:
            raise SystemExit(
                f"ERROR: missing TSV header: {path}"
            )

        return reader.fieldnames, rows


def write_tsv(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(
        path,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


def require_unique(rows, field, label):
    counts = Counter(
        row[field]
        for row in rows
    )

    duplicates = sorted(
        value
        for value, count in counts.items()
        if value and count > 1
    )

    if duplicates:
        raise SystemExit(
            f"ERROR: duplicate {label}: "
            + ", ".join(duplicates)
        )


def integer(value, label):
    try:
        return int(value)
    except Exception:
        raise SystemExit(
            f"ERROR: invalid integer "
            f"{label}={value!r}"
        )


def normalized_dict(row, fields):
    return {
        field: row.get(field, "")
        for field in fields
    }


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--observations",
        required=True,
    )
    parser.add_argument(
        "--loci",
        required=True,
    )
    parser.add_argument(
        "--gold",
        required=True,
    )
    parser.add_argument(
        "--zhao-observations",
        required=True,
    )
    parser.add_argument(
        "--zhao-loci",
        required=True,
    )
    parser.add_argument(
        "--output-dir",
        required=True,
    )

    args = parser.parse_args()

    obs_fields, input_obs = read_tsv(
        args.observations
    )
    loci_fields, input_loci = read_tsv(
        args.loci
    )
    gold_fields, input_gold = read_tsv(
        args.gold
    )

    zobs_fields, zhao_obs = read_tsv(
        args.zhao_observations
    )
    zloci_fields, zhao_loci = read_tsv(
        args.zhao_loci
    )

    outdir = Path(args.output_dir)

    # ========================================================
    # Schema
    # ========================================================

    if zobs_fields != obs_fields:
        raise SystemExit(
            "ERROR: Zhao observation schema differs "
            "from canonical observation schema"
        )

    if zloci_fields != loci_fields:
        raise SystemExit(
            "ERROR: Zhao locus schema differs "
            "from canonical locus schema"
        )

    if gold_fields != loci_fields:
        raise SystemExit(
            "ERROR: binary Gold schema differs "
            "from canonical locus schema"
        )

    # ========================================================
    # Zhao input invariants
    # ========================================================

    if len(zhao_obs) != 12:
        raise SystemExit(
            f"ERROR: expected 12 Zhao observations; "
            f"observed {len(zhao_obs)}"
        )

    if len(zhao_loci) != 3:
        raise SystemExit(
            f"ERROR: expected 3 Zhao loci; "
            f"observed {len(zhao_loci)}"
        )

    require_unique(
        zhao_obs,
        "observation_id",
        "Zhao observation_id",
    )

    require_unique(
        zhao_loci,
        "canonical_locus_id",
        "Zhao canonical_locus_id",
    )

    obs_roles = Counter(
        x["benchmark_role"]
        for x in zhao_obs
    )

    if dict(obs_roles) != EXPECTED_ZHAO_OBS_ROLES:
        raise SystemExit(
            "ERROR: unexpected Zhao observation roles: "
            f"{dict(obs_roles)}"
        )

    obs_classes = Counter(
        x["evidence_class"]
        for x in zhao_obs
    )

    if dict(obs_classes) != EXPECTED_ZHAO_CLASSES:
        raise SystemExit(
            "ERROR: unexpected Zhao evidence classes: "
            f"{dict(obs_classes)}"
        )

    locus_roles = {
        x["locus_name"]: x["benchmark_role"]
        for x in zhao_loci
    }

    if locus_roles != EXPECTED_ZHAO_LOCUS_ROLES:
        raise SystemExit(
            "ERROR: unexpected Zhao locus roles: "
            f"{locus_roles}"
        )

    for row in zhao_obs:
        oid = row["observation_id"]

        if row["study_id"] != "ZHAO2018":
            raise SystemExit(
                f"ERROR: {oid}: study_id is not ZHAO2018"
            )

        if row["target_assembly"] != TARGET_ASSEMBLY:
            raise SystemExit(
                f"ERROR: {oid}: unexpected target assembly"
            )

        if row["mapping_status"] != "mapped_unique":
            raise SystemExit(
                f"ERROR: {oid}: mapping is not mapped_unique"
            )

        if row["mapping_confidence"] != "high":
            raise SystemExit(
                f"ERROR: {oid}: mapping confidence is not high"
            )

        start = integer(
            row["target_start"],
            f"{oid}.target_start",
        )
        end = integer(
            row["target_end"],
            f"{oid}.target_end",
        )

        if end != start + 1:
            raise SystemExit(
                f"ERROR: {oid}: expected adjacent bases "
                "flanking an interbase Cas9 cut"
            )

    # ========================================================
    # Identify Zhao IDs from authoritative curated inputs.
    # This makes the builder idempotent:
    # - pre-Zhao input: no matching IDs are stripped
    # - post-Zhao input: existing Zhao rows are stripped
    #   and reconstructed from the curated Zhao inputs
    # ========================================================

    zhao_obs_ids = {
        x["observation_id"]
        for x in zhao_obs
    }

    zhao_locus_ids = {
        x["canonical_locus_id"]
        for x in zhao_loci
    }

    preview_obs_by_id = {
        x["observation_id"]: x
        for x in zhao_obs
    }

    preview_loci_by_id = {
        x["canonical_locus_id"]: x
        for x in zhao_loci
    }

    unknown_zhao_obs = [
        x["observation_id"]
        for x in input_obs
        if x.get("study_id") == "ZHAO2018"
        and x["observation_id"] not in zhao_obs_ids
    ]

    if unknown_zhao_obs:
        raise SystemExit(
            "ERROR: canonical table contains unknown Zhao "
            "observations: "
            + ", ".join(sorted(unknown_zhao_obs))
        )

    unknown_zhao_loci = [
        x["canonical_locus_id"]
        for x in input_loci
        if "ZHAO2018" in x.get(
            "supporting_studies",
            ""
        ).split(";")
        and x["canonical_locus_id"]
        not in zhao_locus_ids
    ]

    if unknown_zhao_loci:
        raise SystemExit(
            "ERROR: canonical table contains unknown Zhao loci: "
            + ", ".join(sorted(unknown_zhao_loci))
        )

    existing_zhao_obs = [
        x
        for x in input_obs
        if x["observation_id"] in zhao_obs_ids
    ]

    existing_zhao_loci = [
        x
        for x in input_loci
        if x["canonical_locus_id"] in zhao_locus_ids
    ]

    existing_zhao_gold = [
        x
        for x in input_gold
        if x["canonical_locus_id"] in zhao_locus_ids
    ]

    state = (
        len(existing_zhao_obs),
        len(existing_zhao_loci),
        len(existing_zhao_gold),
    )

    if state not in {
        (0, 0, 0),
        (12, 3, 1),
    }:
        raise SystemExit(
            "ERROR: partial/inconsistent Zhao integration "
            f"detected: obs/loci/gold={state}"
        )

    # If Zhao is already present, require exact agreement
    # with the authoritative curated Zhao inputs.
    if state == (12, 3, 1):

        for row in existing_zhao_obs:
            expected = preview_obs_by_id[
                row["observation_id"]
            ]

            if normalized_dict(
                row,
                obs_fields,
            ) != normalized_dict(
                expected,
                obs_fields,
            ):
                raise SystemExit(
                    "ERROR: existing canonical Zhao observation "
                    "differs from curated Zhao input: "
                    f"{row['observation_id']}"
                )

        for row in existing_zhao_loci:
            expected = preview_loci_by_id[
                row["canonical_locus_id"]
            ]

            if normalized_dict(
                row,
                loci_fields,
            ) != normalized_dict(
                expected,
                loci_fields,
            ):
                raise SystemExit(
                    "ERROR: existing canonical Zhao locus "
                    "differs from curated Zhao input: "
                    f"{row['canonical_locus_id']}"
                )

        expected_gold_ids = {
            row["canonical_locus_id"]
            for row in zhao_loci
            if row["benchmark_role"]
            in {"positive", "negative"}
        }

        observed_gold_ids = {
            row["canonical_locus_id"]
            for row in existing_zhao_gold
        }

        if observed_gold_ids != expected_gold_ids:
            raise SystemExit(
                "ERROR: existing Zhao binary-Gold membership "
                "does not match curated Zhao locus roles"
            )

    # ========================================================
    # Strip Zhao if already present.
    # This is the idempotence step.
    # ========================================================

    base_obs = [
        x
        for x in input_obs
        if x["observation_id"] not in zhao_obs_ids
    ]

    base_loci = [
        x
        for x in input_loci
        if x["canonical_locus_id"] not in zhao_locus_ids
    ]

    # Gold is rebuilt later from loci, so no base_gold
    # manipulation is needed beyond provenance reporting.
    base_gold = [
        x
        for x in input_gold
        if x["canonical_locus_id"] not in zhao_locus_ids
    ]

    # ========================================================
    # Collision checks against non-Zhao benchmark content
    # ========================================================

    require_unique(
        base_obs,
        "observation_id",
        "base observation_id",
    )

    require_unique(
        base_loci,
        "canonical_locus_id",
        "base canonical_locus_id",
    )

    existing_obs_ids = {
        x["observation_id"]
        for x in base_obs
    }

    if existing_obs_ids & zhao_obs_ids:
        raise SystemExit(
            "ERROR: unexpected Zhao observation-ID collision"
        )

    existing_locus_ids = {
        x["canonical_locus_id"]
        for x in base_loci
    }

    if existing_locus_ids & zhao_locus_ids:
        raise SystemExit(
            "ERROR: unexpected Zhao locus-ID collision"
        )

    # Physical overlap with existing non-Zhao canonical loci.
    for new in zhao_loci:

        ns = integer(
            new["target_start"],
            new["canonical_locus_id"],
        )
        ne = integer(
            new["target_end"],
            new["canonical_locus_id"],
        )

        for old in base_loci:

            if (
                new["target_assembly"]
                != old["target_assembly"]
            ):
                continue

            if (
                new["target_seqname"]
                != old["target_seqname"]
            ):
                continue

            os = integer(
                old["target_start"],
                old["canonical_locus_id"],
            )
            oe = integer(
                old["target_end"],
                old["canonical_locus_id"],
            )

            if ns <= oe and os <= ne:
                raise SystemExit(
                    "ERROR: Zhao locus overlaps an existing "
                    "non-Zhao canonical locus: "
                    f"{new['canonical_locus_id']} vs "
                    f"{old['canonical_locus_id']}"
                )

    # ========================================================
    # Reconstruct observations
    # ========================================================

    staged_obs = (
        base_obs
        + zhao_obs
    )

    require_unique(
        staged_obs,
        "observation_id",
        "final observation_id",
    )

    # ========================================================
    # Reconstruct canonical loci
    # ========================================================

    staged_loci = (
        base_loci
        + zhao_loci
    )

    require_unique(
        staged_loci,
        "canonical_locus_id",
        "final canonical_locus_id",
    )

    staged_loci.sort(
        key=lambda row: (
            row["target_seqname"],
            integer(
                row["target_start"],
                row["canonical_locus_id"],
            ),
            integer(
                row["target_end"],
                row["canonical_locus_id"],
            ),
            row["canonical_locus_id"],
        )
    )

    # ========================================================
    # Rebuild binary Gold from canonical loci
    # ========================================================

    staged_gold = [
        row
        for row in staged_loci
        if row["benchmark_role"]
        in {"positive", "negative"}
    ]

    if any(
        row["benchmark_role"] == "support_only"
        for row in staged_gold
    ):
        raise SystemExit(
            "ERROR: support_only locus entered binary Gold"
        )

    zhao_gold = [
        row
        for row in staged_gold
        if row["canonical_locus_id"]
        in zhao_locus_ids
    ]

    if len(zhao_gold) != 1:
        raise SystemExit(
            "ERROR: Zhao must contribute exactly one "
            "binary-Gold locus"
        )

    if (
        zhao_gold[0]["canonical_locus_id"]
        != "CHO_BM_ZHAO_001"
        or zhao_gold[0]["locus_name"]
        != "C12orf35"
    ):
        raise SystemExit(
            "ERROR: unexpected Zhao Gold locus"
        )

    # ========================================================
    # Write outputs
    # ========================================================

    write_tsv(
        outdir / "benchmark_observations.tsv",
        obs_fields,
        staged_obs,
    )

    write_tsv(
        outdir / "benchmark_loci.tsv",
        loci_fields,
        staged_loci,
    )

    write_tsv(
        outdir / "benchmark_gold_binary.tsv",
        gold_fields,
        staged_gold,
    )

    print("PASS: Zhao curated build")
    print(
        "input_state="
        f"{len(input_obs)}/"
        f"{len(input_loci)}/"
        f"{len(input_gold)}"
    )
    print(
        "zhao_existing_state="
        f"{state[0]}/{state[1]}/{state[2]}"
    )
    print(
        "non_zhao_base="
        f"{len(base_obs)}/"
        f"{len(base_loci)}/"
        f"{len(base_gold)}"
    )
    print(
        "output_state="
        f"{len(staged_obs)}/"
        f"{len(staged_loci)}/"
        f"{len(staged_gold)}"
    )
    print(
        "Zhao_binary_gold="
        "CHO_BM_ZHAO_001:C12orf35"
    )


if __name__ == "__main__":
    main()
