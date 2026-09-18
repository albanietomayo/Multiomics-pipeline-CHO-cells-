#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]

OBS_INPUT = (
    ROOT / "intermediate" / "dhiman2020" /
    "dhiman_observed_sites.tsv"
)

MAP_INPUT = (
    ROOT / "intermediate" / "dhiman2020" /
    "liftover" / "dhiman_loci_mapped.tsv"
)

OBS_OUTPUT = (
    ROOT / "curated" / "benchmark_observations.tsv"
)

LOCI_OUTPUT = (
    ROOT / "curated" / "benchmark_loci.tsv"
)

GOLD_OUTPUT = (
    ROOT / "curated" / "benchmark_gold_binary.tsv"
)

BED_OUTPUT = (
    ROOT / "curated" / "benchmark_loci.bed"
)

GOLD_BED_OUTPUT = (
    ROOT / "curated" / "benchmark_gold_binary.bed"
)

SUMMARY_OUTPUT = (
    ROOT / "reports" / "dhiman2020" /
    "curation_summary.tsv"
)

TARGET = "GCF_003668045.3"


def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


observed = read_tsv(OBS_INPUT)
mapped = read_tsv(MAP_INPUT)

# ------------------------------------------------------------
# Fail-closed scientific QC
# ------------------------------------------------------------

if len(observed) != 17:
    raise SystemExit(
        f"ERROR: expected 17 experimental observations, "
        f"found {len(observed)}"
    )

if len(mapped) != 13:
    raise SystemExit(
        f"ERROR: expected 13 unique loci, found {len(mapped)}"
    )

if any(x["mapping_status"] != "mapped_unique" for x in mapped):
    raise SystemExit(
        "ERROR: not all Dhiman loci are mapped_unique"
    )

if any(x["target_assembly_accession"] != TARGET for x in mapped):
    raise SystemExit(
        "ERROR: unexpected target assembly"
    )

for row in mapped:
    if row["span_ratio"] != "1.000000":
        raise SystemExit(
            "ERROR: target/source span discrepancy for "
            + row["source_locus_id"]
        )

# ------------------------------------------------------------
# Source-coordinate → harmonized-locus lookup
# ------------------------------------------------------------

mapped_by_source = {}

for row in mapped:
    key = (
        row["source_seqname"],
        int(row["source_start"]),
        int(row["source_end"]),
    )

    if key in mapped_by_source:
        raise SystemExit(
            f"ERROR: duplicated source coordinate {key}"
        )

    mapped_by_source[key] = row

# ------------------------------------------------------------
# Official observation table
# One row = one biological experimental observation.
# ------------------------------------------------------------

obs_fields = [
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
]

new_observations = []

for source in observed:
    key = (
        source["source_seqname"],
        int(source["source_start"]),
        int(source["source_end"]),
    )

    if key not in mapped_by_source:
        raise SystemExit(
            f"ERROR: observation coordinate has no liftover: {key}"
        )

    target = mapped_by_source[key]

    role = source["benchmark_role_observation"]

    if role not in {"positive", "negative"}:
        raise SystemExit(
            f"ERROR: unexpected observation role: {role}"
        )

    notes = "; ".join(
        [
            f"copy_number={source['copy_number']}",
            f"selected_passage={source['selected_passage']}",
            f"variant_profile={source['variant_profile']}",
            (
                "expression_peak_distance="
                f"{source['expression_peak_distance']}"
            ),
            f"chromatin_state={source['chromatin_state']}",
            f"genic_context={source['genic_context']}",
        ]
    )

    new_observations.append(
        {
            "observation_id": source["observation_id"],
            "study_id": "DHIMAN2020",
            "source_locus_id": target["source_locus_id"],
            "locus_name": target["source_locus_id"],
            "cell_line": source["cell_line"],
            "integration_method": (
                "Random transgene integration; "
                "integration site resolved by TLA-seq"
            ),
            "transgene_or_product": "",
            "phenotype": (
                f"{source['stability']} transgene-expression "
                "stability phenotype"
            ),
            "evidence_class": source["evidence_class"],
            "evidence_tier": "Tier 2",
            "benchmark_role": role,
            "expression_evidence": "",
            "productivity_evidence": "",
            "stability_evidence": source["stability"],
            "stability_duration": "",
            "source_table_or_figure": (
                "Supplementary Tables 1 and 4"
            ),
            "source_assembly": (
                "GCF_003668045.1 (CriGri-PICR)"
            ),
            "source_coordinate_system": (
                "Published genomic coordinates; "
                "treated as 1-based closed for sequence extraction"
            ),
            "source_seqname": source["source_seqname"],
            "source_start": source["source_start"],
            "source_end": source["source_end"],
            "source_strand": "",
            "target_assembly": TARGET,
            "target_seqname": target["target_seqname"],
            "target_start": target["target_start"],
            "target_end": target["target_end"],
            "mapping_method": target["mapping_method"],
            "mapping_status": "mapped_unique",
            "mapping_confidence": "high",
            "notes": notes,
        }
    )

# Idempotent upsert.
existing_obs = []

if OBS_OUTPUT.exists():
    existing_obs = read_tsv(OBS_OUTPUT)

existing_obs = [
    x for x in existing_obs
    if x.get("study_id") != "DHIMAN2020"
]

all_obs = existing_obs + new_observations

all_obs.sort(
    key=lambda x: (
        x["study_id"],
        x["source_locus_id"],
        x["observation_id"],
    )
)

write_tsv(
    OBS_OUTPUT,
    all_obs,
    obs_fields,
)

# ------------------------------------------------------------
# Canonical locus table
#
# Important:
# stable_only   -> positive
# unstable_only -> negative
# context_conflict -> support_only
#
# We do NOT force contradictory loci into a binary label.
# ------------------------------------------------------------

locus_fields = [
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
]

new_loci = []

for row in sorted(
    mapped,
    key=lambda x: x["source_locus_id"],
):
    suffix = row["source_locus_id"].split("_L")[-1]

    canonical_id = f"CHO_BM_DHIMAN_{suffix}"

    status = row["locus_status"]

    if status == "stable_only":
        role = "positive"
    elif status == "unstable_only":
        role = "negative"
    elif status == "context_conflict":
        role = "support_only"
    else:
        raise SystemExit(
            f"ERROR: unexpected locus_status={status}"
        )

    notes = "; ".join(
        [
            f"source_locus_id={row['source_locus_id']}",
            f"source_assembly=GCF_003668045.1",
            f"cell_lines={row['cell_lines']}",
            f"observed_stabilities={row['stabilities']}",
            f"locus_status={status}",
            (
                "binary_gold_eligible="
                f"{row['gold_binary_eligible']}"
            ),
            (
                "mapping="
                "unique boundary-flank sequence alignment; "
                "source/target span ratio 1.0"
            ),
        ]
    )

    new_loci.append(
        {
            "canonical_locus_id": canonical_id,
            "locus_name": row["source_locus_id"],
            "benchmark_role": role,
            "best_evidence_tier": "Tier 2",
            "n_observations": row["n_observations"],
            "target_assembly": TARGET,
            "target_seqname": row["target_seqname"],
            "target_start": row["target_start"],
            "target_end": row["target_end"],
            "mapping_confidence": "high",
            "supporting_studies": "DHIMAN2020",
            "notes": notes,
        }
    )

existing_loci = []

if LOCI_OUTPUT.exists():
    existing_loci = read_tsv(LOCI_OUTPUT)

existing_loci = [
    x for x in existing_loci
    if not x.get(
        "canonical_locus_id", ""
    ).startswith("CHO_BM_DHIMAN_")
]

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
    LOCI_OUTPUT,
    all_loci,
    locus_fields,
)

# ------------------------------------------------------------
# Gold binary subset
# ------------------------------------------------------------

gold = [
    x for x in all_loci
    if x["benchmark_role"] in {
        "positive",
        "negative",
    }
]

write_tsv(
    GOLD_OUTPUT,
    gold,
    locus_fields,
)

# ------------------------------------------------------------
# BED conversion
# Internal TSV: 1-based closed
# BED:          0-based half-open
#
# Strand "." intentionally:
# liftover orientation is assembly-mapping orientation,
# NOT a biological strand of the benchmark locus.
# ------------------------------------------------------------

def write_bed(path, rows):
    ordered = sorted(
        rows,
        key=lambda x: (
            x["target_seqname"],
            int(x["target_start"]),
            int(x["target_end"]),
        )
    )

    with path.open("w") as fh:
        for row in ordered:
            start0 = int(row["target_start"]) - 1
            end0 = int(row["target_end"])

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


write_bed(BED_OUTPUT, all_loci)
write_bed(GOLD_BED_OUTPUT, gold)

# ------------------------------------------------------------
# Final QC
# ------------------------------------------------------------

dhiman_loci = [
    x for x in all_loci
    if x["canonical_locus_id"].startswith(
        "CHO_BM_DHIMAN_"
    )
]

roles = Counter(
    x["benchmark_role"]
    for x in dhiman_loci
)

obs_roles = Counter(
    x["benchmark_role"]
    for x in new_observations
)

expected_locus_roles = {
    "positive": 6,
    "negative": 4,
    "support_only": 3,
}

if dict(roles) != expected_locus_roles:
    raise SystemExit(
        f"ERROR: unexpected locus role counts {dict(roles)}"
    )

if obs_roles != Counter(
    {
        "positive": 10,
        "negative": 7,
    }
):
    raise SystemExit(
        f"ERROR: unexpected observation counts "
        f"{dict(obs_roles)}"
    )

dhiman_gold = [
    x for x in dhiman_loci
    if x["benchmark_role"] in {
        "positive",
        "negative",
    }
]

if len(dhiman_gold) != 10:
    raise SystemExit(
        f"ERROR: expected 10 binary gold loci, "
        f"found {len(dhiman_gold)}"
    )

summary = [
    ("study", "DHIMAN2020"),
    ("experimental_observations", len(new_observations)),
    ("positive_observations", obs_roles["positive"]),
    ("negative_observations", obs_roles["negative"]),
    ("unique_harmonized_loci", len(dhiman_loci)),
    ("mapped_unique_loci", len(dhiman_loci)),
    ("positive_loci", roles["positive"]),
    ("negative_loci", roles["negative"]),
    ("context_conflict_loci", roles["support_only"]),
    ("gold_binary_loci", len(dhiman_gold)),
    ("target_assembly", TARGET),
]

with SUMMARY_OUTPUT.open("w", newline="") as fh:
    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(["metric", "value"])
    writer.writerows(summary)

print("DHIMAN CURATION: PASS")
for key, value in summary:
    print(f"{key:28s} {value}")

print()
print(f"Observations: {OBS_OUTPUT}")
print(f"All loci:     {LOCI_OUTPUT}")
print(f"Gold binary:  {GOLD_OUTPUT}")
print(f"All BED:      {BED_OUTPUT}")
print(f"Gold BED:     {GOLD_BED_OUTPUT}")
