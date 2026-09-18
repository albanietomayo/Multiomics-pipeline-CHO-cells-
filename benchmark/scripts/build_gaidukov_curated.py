#!/usr/bin/env python3

import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MAPPED = (
    ROOT / "intermediate" / "gaidukov2018" /
    "liftover" / "gaidukov_loci_mapped.tsv"
)

ASM5_DIAG = (
    ROOT / "reports" / "gaidukov2018" /
    "gaidukov_unresolved_asm5_diagnostics.tsv"
)

OBS_FILE = ROOT / "curated" / "benchmark_observations.tsv"
LOCI_FILE = ROOT / "curated" / "benchmark_loci.tsv"

GOLD_FILE = ROOT / "curated" / "benchmark_gold_binary.tsv"
BED_FILE = ROOT / "curated" / "benchmark_loci.bed"
GOLD_BED = ROOT / "curated" / "benchmark_gold_binary.bed"

UNRESOLVED_FILE = (
    ROOT / "reports" / "gaidukov2018" /
    "unresolved_loci.tsv"
)

SUMMARY_FILE = (
    ROOT / "reports" / "gaidukov2018" /
    "curation_summary.tsv"
)

TARGET = "GCF_003668045.3"


def read_tsv(path):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


mapped_fields, rows = read_tsv(MAPPED)

if len(rows) != 21:
    raise SystemExit(
        f"ERROR: expected 21 Gaidukov loci, found {len(rows)}"
    )

status_counts = Counter(
    x["mapping_status"]
    for x in rows
)

if status_counts != Counter(
    {
        "mapped_unique": 17,
        "unresolved": 4,
    }
):
    raise SystemExit(
        f"ERROR: unexpected mapping counts: {status_counts}"
    )

unresolved_expected = {
    "GAIDUKOV2018_L004",
    "GAIDUKOV2018_L006",
    "GAIDUKOV2018_L012",
    "GAIDUKOV2018_L020",
}

unresolved_actual = {
    x["source_locus_id"]
    for x in rows
    if x["mapping_status"] == "unresolved"
}

if unresolved_actual != unresolved_expected:
    raise SystemExit(
        f"ERROR: unresolved loci differ: "
        f"{sorted(unresolved_actual)}"
    )

# ------------------------------------------------------------
# Existing benchmark schema
# ------------------------------------------------------------

obs_fields, existing_obs = read_tsv(OBS_FILE)
locus_fields, existing_loci = read_tsv(LOCI_FILE)

required_obs = {
    "observation_id",
    "study_id",
    "source_locus_id",
    "locus_name",
    "cell_line",
    "integration_method",
    "transgene_or_product",
    "phenotype",
    "evidence_class",
    "evidence_tier",
    "benchmark_role",
    "expression_evidence",
    "productivity_evidence",
    "stability_evidence",
    "stability_duration",
    "source_table_or_figure",
    "source_assembly",
    "source_coordinate_system",
    "source_seqname",
    "source_start",
    "source_end",
    "source_strand",
    "target_assembly",
    "target_seqname",
    "target_start",
    "target_end",
    "mapping_method",
    "mapping_status",
    "mapping_confidence",
    "notes",
}

required_locus = {
    "canonical_locus_id",
    "locus_name",
    "benchmark_role",
    "best_evidence_tier",
    "n_observations",
    "target_assembly",
    "target_seqname",
    "target_start",
    "target_end",
    "mapping_confidence",
    "supporting_studies",
    "notes",
}

if not required_obs.issubset(obs_fields):
    raise SystemExit(
        "ERROR: benchmark_observations.tsv schema changed"
    )

if not required_locus.issubset(locus_fields):
    raise SystemExit(
        "ERROR: benchmark_loci.tsv schema changed"
    )

# Idempotent replacement of Gaidukov observations.
existing_obs = [
    x for x in existing_obs
    if x.get("study_id") != "GAIDUKOV2018"
]

# Remove previously generated Gaidukov canonical loci,
# if script is rerun.
existing_loci = [
    x for x in existing_loci
    if not x.get(
        "canonical_locus_id", ""
    ).startswith("CHO_BM_GAIDUKOV_")
]

# ------------------------------------------------------------
# Observation table: preserve ALL 21 published observations.
#
# The four unresolved loci remain source-level experimental
# observations but are excluded from the harmonized locus table.
# ------------------------------------------------------------

new_obs = []

for row in rows:

    locus_id = row["source_locus_id"]
    suffix = locus_id.rsplit("_L", 1)[1]

    mapped_ok = (
        row["mapping_status"] == "mapped_unique"
    )

    source_orientation = {
        "forward": "+",
        "reverse": "-",
    }.get(row["orientation"], "")

    # Exact lentiviral insertion sites are screen evidence.
    # Nearby CRISPR reconstruction is retained as support but
    # does not convert the precise lentiviral nucleotide into
    # a targeted CRISPR coordinate.
    if locus_id == "GAIDUKOV2018_L021":
        evidence_class = "targeted_validated"
        integration_method = (
            "CRISPR/Cas9 targeted integration at "
            "putative CHO Rosa26 locus"
        )
    else:
        evidence_class = "screen_validated"
        integration_method = (
            "Lentiviral random-integration screen"
        )

    notes_parts = [
        f"paper_scaffold=scaffold{row['source_scaffold_number']}",
        f"alias={row['alias']}",
        f"gene={row['gene']}",
        f"region={row['region']}",
        f"discovery_mode={row['discovery_mode']}",
        (
            "targeted_reconstruction="
            f"{row['targeted_reconstruction']}"
        ),
    ]

    if row["targeted_clone"]:
        notes_parts.append(
            f"targeted_clone={row['targeted_clone']}"
        )

    if row["targeted_site_relation"]:
        notes_parts.append(
            "targeted_site_relation="
            f"{row['targeted_site_relation']}"
        )

    if not mapped_ok:
        notes_parts.append(
            "excluded_from_harmonized_locus_table="
            "strict_liftover_failure"
        )

    record = {
        field: ""
        for field in obs_fields
    }

    record.update(
        {
            "observation_id": (
                f"GAIDUKOV2018_OBS{suffix}"
            ),
            "study_id": "GAIDUKOV2018",
            "source_locus_id": locus_id,
            "locus_name": row["lp_clone"],
            "cell_line": "CHO-K1",
            "integration_method": integration_method,
            "transgene_or_product": (
                "Landing-pad reporter cassette"
            ),
            "phenotype": (
                "Stable heterologous transgene expression"
            ),
            "evidence_class": evidence_class,
            # Conservative because the CRISPR reconstruction
            # can be nearby rather than at the exact viral base.
            "evidence_tier": "Tier 2",
            "benchmark_role": "positive",
            "expression_evidence": (
                "Stable reporter expression"
            ),
            "productivity_evidence": "",
            "stability_evidence": "Stable",
            "stability_duration": "",
            "source_table_or_figure": (
                "Gaidukov et al. 2018 Table 1; "
                "Figures 1-2"
            ),
            "source_assembly": (
                "GCF_000223135.1 (CriGri_1.0)"
            ),
            "source_coordinate_system": (
                "Published scaffold point position; "
                "interpreted as 1-based"
            ),
            "source_seqname": row["source_seqname"],
            "source_start": row["source_start"],
            "source_end": row["source_end"],
            "source_strand": source_orientation,
            "target_assembly": TARGET,
            "target_seqname": (
                row["target_seqname"]
                if mapped_ok
                else ""
            ),
            "target_start": (
                row["target_start"]
                if mapped_ok
                else ""
            ),
            "target_end": (
                row["target_end"]
                if mapped_ok
                else ""
            ),
            "mapping_method": row["mapping_method"],
            "mapping_status": row["mapping_status"],
            "mapping_confidence": (
                "high"
                if mapped_ok
                else "unresolved"
            ),
            "notes": "; ".join(notes_parts),
        }
    )

    new_obs.append(record)

if len(new_obs) != 21:
    raise SystemExit(
        "ERROR: expected 21 new observations"
    )

# ------------------------------------------------------------
# Harmonized loci: ONLY the 17 mapped_unique loci.
# ------------------------------------------------------------

mapped_unique = [
    x for x in rows
    if x["mapping_status"] == "mapped_unique"
]

if len(mapped_unique) != 17:
    raise SystemExit(
        "ERROR: expected 17 mapped_unique loci"
    )

# Ensure Gaidukov coordinates themselves are unique.
coords = [
    (
        x["target_seqname"],
        int(x["target_start"]),
        int(x["target_end"]),
    )
    for x in mapped_unique
]

if len(coords) != len(set(coords)):
    raise SystemExit(
        "ERROR: duplicate target coordinate "
        "within Gaidukov loci"
    )

# Existing exact coordinate lookup.
existing_by_coord = {}

for row in existing_loci:

    key = (
        row["target_seqname"],
        int(row["target_start"]),
        int(row["target_end"]),
    )

    if key in existing_by_coord:
        raise SystemExit(
            f"ERROR: pre-existing duplicate benchmark "
            f"coordinate {key}"
        )

    existing_by_coord[key] = row


new_loci = []
merged_existing = 0

for row in mapped_unique:

    locus_id = row["source_locus_id"]
    suffix = locus_id.rsplit("_L", 1)[1]

    coord = (
        row["target_seqname"],
        int(row["target_start"]),
        int(row["target_end"]),
    )

    notes = "; ".join(
        [
            f"source_locus_id={locus_id}",
            (
                "source_assembly="
                "GCF_000223135.1"
            ),
            (
                "source_scaffold="
                f"scaffold{row['source_scaffold_number']}"
            ),
            f"source_position={row['source_position']}",
            f"gene={row['gene']}",
            f"region={row['region']}",
            (
                "targeted_reconstruction="
                f"{row['targeted_reconstruction']}"
            ),
            (
                "targeted_clone="
                f"{row['targeted_clone']}"
            ),
            (
                "mapping=unique sequence liftover "
                "validated independently with paftools.js"
            ),
        ]
    )

    if coord in existing_by_coord:

        old = existing_by_coord[coord]

        if old["benchmark_role"] != "positive":
            raise SystemExit(
                "ERROR: cross-study role conflict at "
                f"{coord}: existing="
                f"{old['benchmark_role']} "
                "vs Gaidukov=positive"
            )

        studies = {
            x
            for x in old[
                "supporting_studies"
            ].split(";")
            if x
        }

        studies.add("GAIDUKOV2018")

        old["supporting_studies"] = ";".join(
            sorted(studies)
        )

        old["n_observations"] = str(
            int(old["n_observations"]) + 1
        )

        old["notes"] = (
            old["notes"]
            + "; exact_coordinate_support="
            + locus_id
        )

        merged_existing += 1
        continue

    record = {
        field: ""
        for field in locus_fields
    }

    record.update(
        {
            "canonical_locus_id": (
                f"CHO_BM_GAIDUKOV_{suffix}"
            ),
            "locus_name": row["lp_clone"],
            "benchmark_role": "positive",
            "best_evidence_tier": "Tier 2",
            "n_observations": "1",
            "target_assembly": TARGET,
            "target_seqname": row["target_seqname"],
            "target_start": row["target_start"],
            "target_end": row["target_end"],
            "mapping_confidence": "high",
            "supporting_studies": "GAIDUKOV2018",
            "notes": notes,
        }
    )

    new_loci.append(record)

all_obs = existing_obs + new_obs

all_obs.sort(
    key=lambda x: (
        x["study_id"],
        x["source_locus_id"],
        x["observation_id"],
    )
)

all_loci = existing_loci + new_loci

all_loci.sort(
    key=lambda x: (
        x["target_seqname"],
        int(x["target_start"]),
        int(x["target_end"]),
        x["canonical_locus_id"],
    )
)

write_tsv(
    OBS_FILE,
    obs_fields,
    all_obs,
)

write_tsv(
    LOCI_FILE,
    locus_fields,
    all_loci,
)

# ------------------------------------------------------------
# Gold binary = harmonized positive + negative only.
# support_only remains outside.
# ------------------------------------------------------------

gold = [
    x for x in all_loci
    if x["benchmark_role"] in {
        "positive",
        "negative",
    }
]

write_tsv(
    GOLD_FILE,
    locus_fields,
    gold,
)

# ------------------------------------------------------------
# Regenerate BEDs from canonical 1-based closed coordinates.
# ------------------------------------------------------------

def write_bed(path, records):

    ordered = sorted(
        records,
        key=lambda x: (
            x["target_seqname"],
            int(x["target_start"]),
            int(x["target_end"]),
            x["canonical_locus_id"],
        )
    )

    with path.open("w") as fh:

        for row in ordered:

            start0 = (
                int(row["target_start"]) - 1
            )

            end0 = int(
                row["target_end"]
            )

            fh.write(
                "\t".join(
                    [
                        row["target_seqname"],
                        str(start0),
                        str(end0),
                        row["canonical_locus_id"],
                        "0",
                        ".",
                    ]
                )
                + "\n"
            )


write_bed(
    BED_FILE,
    all_loci,
)

write_bed(
    GOLD_BED,
    gold,
)

# ------------------------------------------------------------
# Unresolved report with strict failure reasons.
# ------------------------------------------------------------

diag_fields, diag = read_tsv(
    ASM5_DIAG
)

reasons = {}

for row in diag:

    locus = row["source_locus_id"]

    if locus not in unresolved_expected:
        continue

    reason = row["failure_reasons"]

    reasons.setdefault(
        locus,
        [],
    )

    if reason not in reasons[locus]:
        reasons[locus].append(reason)

unresolved_report = []

for row in rows:

    if row["mapping_status"] != "unresolved":
        continue

    unresolved_report.append(
        {
            "source_locus_id": row["source_locus_id"],
            "lp_clone": row["lp_clone"],
            "source_seqname": row["source_seqname"],
            "source_position": row["source_position"],
            "benchmark_role": "positive",
            "mapping_status": "unresolved",
            "strict_failure_reasons": " | ".join(
                reasons.get(
                    row["source_locus_id"],
                    ["not_available"],
                )
            ),
            "canonical_inclusion": "FALSE",
        }
    )

unresolved_fields = [
    "source_locus_id",
    "lp_clone",
    "source_seqname",
    "source_position",
    "benchmark_role",
    "mapping_status",
    "strict_failure_reasons",
    "canonical_inclusion",
]

write_tsv(
    UNRESOLVED_FILE,
    unresolved_fields,
    unresolved_report,
)

# ------------------------------------------------------------
# Final QC / summary
# ------------------------------------------------------------

gaidukov_obs = [
    x for x in all_obs
    if x["study_id"] == "GAIDUKOV2018"
]

mapped_obs = [
    x for x in gaidukov_obs
    if x["mapping_status"] == "mapped_unique"
]

unresolved_obs = [
    x for x in gaidukov_obs
    if x["mapping_status"] == "unresolved"
]

targeted_mapped = [
    x for x in mapped_unique
    if x["targeted_reconstruction"] == "TRUE"
]

roles = Counter(
    x["benchmark_role"]
    for x in all_loci
)

summary = [
    ("study", "GAIDUKOV2018"),
    ("published_physical_loci", 21),
    ("experimental_observations_added", len(gaidukov_obs)),
    ("mapped_unique_observations", len(mapped_obs)),
    ("unresolved_observations", len(unresolved_obs)),
    ("targeted_CRISPR_regions_total", 7),
    ("targeted_CRISPR_regions_mapped", len(targeted_mapped)),
    ("new_canonical_loci_added", len(new_loci)),
    ("exact_existing_loci_merged", merged_existing),
    ("total_benchmark_observations", len(all_obs)),
    ("total_harmonized_loci", len(all_loci)),
    ("total_positive_loci", roles["positive"]),
    ("total_negative_loci", roles["negative"]),
    ("total_support_only_loci", roles["support_only"]),
    ("total_gold_binary_loci", len(gold)),
    ("target_assembly", TARGET),
]

with SUMMARY_FILE.open(
    "w",
    newline="",
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(
        ["metric", "value"]
    )

    writer.writerows(summary)

if len(gaidukov_obs) != 21:
    raise SystemExit(
        "ERROR: Gaidukov observations != 21"
    )

if len(mapped_obs) != 17:
    raise SystemExit(
        "ERROR: mapped observations != 17"
    )

if len(unresolved_obs) != 4:
    raise SystemExit(
        "ERROR: unresolved observations != 4"
    )

if len(targeted_mapped) != 6:
    raise SystemExit(
        "ERROR: expected 6/7 CRISPR-supported regions mapped"
    )

print("GAIDUKOV CURATION: PASS")

for key, value in summary:
    print(f"{key:34s} {value}")

print()
print(
    "Unresolved loci retained as source-level observations:"
)

for row in unresolved_report:
    print(
        " ",
        row["source_locus_id"],
        row["lp_clone"],
        "->",
        row["strict_failure_reasons"],
    )
