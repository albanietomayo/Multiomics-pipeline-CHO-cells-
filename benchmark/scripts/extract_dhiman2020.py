#!/usr/bin/env python3

import csv
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path, PurePosixPath

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
OFFICE_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"

NS = {"m": MAIN}

SOURCE_ASSEMBLY_NAME = "CriGri-PICR"
SOURCE_ASSEMBLY_ACCESSION = "GCF_003668045.1"

TARGET_ASSEMBLY_NAME = "CriGri-PICRH-1.0"
TARGET_ASSEMBLY_ACCESSION = "GCF_003668045.3"


def all_text(node):
    return "".join(
        x.text or ""
        for x in node.iter()
        if x.tag == f"{{{MAIN}}}t"
    )


def shared_strings(zf):
    path = "xl/sharedStrings.xml"
    if path not in zf.namelist():
        return []

    root = ET.fromstring(zf.read(path))
    return [
        all_text(x)
        for x in root.findall(f"{{{MAIN}}}si")
    ]


def workbook_sheets(zf):
    root = ET.fromstring(zf.read("xl/workbook.xml"))

    relroot = ET.fromstring(
        zf.read("xl/_rels/workbook.xml.rels")
    )

    rels = {
        x.attrib["Id"]: x.attrib["Target"]
        for x in relroot.findall(
            f"{{{PKG_REL}}}Relationship"
        )
    }

    result = {}

    for sheet in root.findall("m:sheets/m:sheet", NS):
        name = sheet.attrib["name"]
        rid = sheet.attrib[f"{{{OFFICE_REL}}}id"]
        target = rels[rid]

        if target.startswith("/"):
            path = target.lstrip("/")
        else:
            path = str(PurePosixPath("xl") / target)

        result[name] = path

    return result


def cell_value(cell, shared):
    typ = cell.attrib.get("t")

    if typ == "inlineStr":
        node = cell.find("m:is", NS)
        return all_text(node) if node is not None else ""

    node = cell.find("m:v", NS)

    if node is None or node.text is None:
        return ""

    value = node.text

    if typ == "s":
        return shared[int(value)]

    return value


def read_sheet(zf, path, shared):
    root = ET.fromstring(zf.read(path))

    rows = {}

    for row in root.findall(".//m:sheetData/m:row", NS):
        rnum = int(row.attrib["r"])
        cells = {}

        for cell in row.findall("m:c", NS):
            ref = cell.attrib["r"]
            col = re.match(r"[A-Z]+", ref).group(0)
            cells[col] = cell_value(cell, shared).strip()

        rows[rnum] = cells

    return rows


def clean_cell_line(value):
    value = value.strip()

    # Remove typographic quote found in Supplementary Table 4.
    value = value.lstrip("‘’'\"")

    return value.strip()


def valid_coord_field(value):
    value = value.strip()

    return value not in {
        "",
        "_",
        "N.A.",
        "NA",
        "N/A",
    }


COORD_RE = re.compile(
    r"^(NW_\d+\.\d+):(\d+)(?:-(\d+))?$"
)


def parse_coordinates(value):
    out = []

    for item in value.split(";"):
        item = item.strip()

        if not item:
            continue

        match = COORD_RE.fullmatch(item)

        if not match:
            raise ValueError(
                f"Unrecognised coordinate: {item!r}"
            )

        seqname = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3) or match.group(2))

        out.append(
            {
                "raw": item,
                "seqname": seqname,
                "start": start,
                "end": end,
            }
        )

    return out


def stability_from_description(description):
    x = description.lower()

    if "unstable" in x:
        return "Unstable"

    if "stable" in x:
        return "Stable"

    raise ValueError(
        f"Cannot determine stability from {description!r}"
    )


def copy_number_from_description(description):
    x = description.lower()

    if "low copy" in x:
        return "Low"

    if "medium copy" in x:
        return "Medium"

    if "high copy" in x:
        return "High"

    if "unknown copy" in x:
        return "Unknown"

    return "Unknown"


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: extract_dhiman2020.py mmc2.xlsx output_dir"
        )

    xlsx = Path(sys.argv[1])
    outdir = Path(sys.argv[2])

    with zipfile.ZipFile(xlsx) as zf:
        shared = shared_strings(zf)
        sheets = workbook_sheets(zf)

        table1 = read_sheet(
            zf,
            sheets["Suppl Table 1"],
            shared,
        )

        table4 = read_sheet(
            zf,
            sheets["Suppl Table 4"],
            shared,
        )

    # ---------------------------------------------------------
    # Supplementary Table 1:
    # identify phenotyped cell lines and P1/P10 coordinates.
    # ---------------------------------------------------------

    cells = {}

    for rownum, row in sorted(table1.items()):
        # Supplementary Table 1 contains the 13 phenotyped
        # experimental cell lines in Excel rows 7-19.
        # Restrict parsing explicitly to this documented block
        # so headers, notes and auxiliary rows can never be
        # interpreted as biological observations.
        if rownum < 7 or rownum > 19:
            continue

        cell_line = clean_cell_line(row.get("A", ""))
        description = row.get("C", "").strip()

        p1 = row.get("E", "").strip()
        p10 = row.get("F", "").strip()

        if not cell_line or not description:
            continue

        if not valid_coord_field(p1):
            continue

        stability = stability_from_description(description)
        copy_number = copy_number_from_description(description)

        # P10 represents the late-passage set where available.
        # If P10 is not available, preserve the P1 set.
        if valid_coord_field(p10):
            selected_passage = "P10"
            selected_raw = p10
        else:
            selected_passage = "P1"
            selected_raw = p1

        selected_coordinates = parse_coordinates(selected_raw)

        cells[cell_line] = {
            "table1_row": rownum,
            "cell_line": cell_line,
            "vector": row.get("B", ""),
            "description": description,
            "copy_number": copy_number,
            "stability": stability,
            "p1_raw": p1,
            "p10_raw": p10,
            "selected_passage": selected_passage,
            "selected_raw": selected_raw,
            "selected_coordinates": selected_coordinates,
            "gene_description": row.get("M", ""),
            "gene_notes": row.get("N", ""),
        }

    if len(cells) != 13:
        raise SystemExit(
            f"ERROR: expected 13 phenotyped cell lines, found {len(cells)}"
        )

    # ---------------------------------------------------------
    # Supplementary Table 4:
    # one row per integration-site profile.
    # ---------------------------------------------------------

    dhiman_profiles = defaultdict(list)
    gaidukov_profiles = []

    for rownum, row in sorted(table4.items()):
        if rownum < 4:
            continue

        cell_line = clean_cell_line(row.get("A", ""))

        if not cell_line:
            continue

        profile = {
            "table4_row": rownum,
            "cell_line": cell_line,
            "copy_number": row.get("B", ""),
            "stability": row.get("C", ""),
            "variant_profile": row.get("D", ""),
            "expression_peak_distance": row.get("E", ""),
            "chromatin_state": row.get("F", ""),
            "genic_context": row.get("G", ""),
        }

        if cell_line in cells:
            dhiman_profiles[cell_line].append(profile)
        else:
            gaidukov_profiles.append(profile)

    n_dhiman_profiles = sum(
        len(v)
        for v in dhiman_profiles.values()
    )

    if n_dhiman_profiles != 17:
        raise SystemExit(
            f"ERROR: expected 17 Dhiman integration-site profiles, "
            f"found {n_dhiman_profiles}"
        )

    if len(gaidukov_profiles) != 20:
        raise SystemExit(
            f"ERROR: expected 20 Gaidukov profiles, "
            f"found {len(gaidukov_profiles)}"
        )

    # ---------------------------------------------------------
    # Match Table 4 site profiles to Table 1 coordinates.
    #
    # The late-passage coordinate set is used when available.
    # Order is required to match exactly.
    # ---------------------------------------------------------

    observations = []

    for cell_line, info in cells.items():
        coords = info["selected_coordinates"]
        profiles = dhiman_profiles[cell_line]

        if len(coords) != len(profiles):
            raise SystemExit(
                "ERROR: coordinate/profile count mismatch for "
                f"{cell_line}: {len(coords)} coordinates vs "
                f"{len(profiles)} profiles"
            )

        for index, (coord, profile) in enumerate(
            zip(coords, profiles),
            start=1,
        ):
            profile_stability = profile["stability"].strip()

            if profile_stability != info["stability"]:
                raise SystemExit(
                    f"ERROR: stability mismatch for {cell_line}: "
                    f"Table1={info['stability']} "
                    f"Table4={profile_stability}"
                )

            role = (
                "positive"
                if info["stability"] == "Stable"
                else "negative"
            )

            evidence_class = (
                "observed_stable"
                if info["stability"] == "Stable"
                else "observed_unstable"
            )

            observation_id = (
                f"DHIMAN2020_{cell_line.replace('.', '_')}_"
                f"SITE{index:02d}"
            )

            observations.append(
                {
                    "observation_id": observation_id,
                    "study_id": "DHIMAN2020",
                    "cell_line": cell_line,
                    "copy_number": info["copy_number"],
                    "stability": info["stability"],
                    "benchmark_role_observation": role,
                    "evidence_class": evidence_class,
                    "selected_passage": info["selected_passage"],
                    "source_assembly_name": SOURCE_ASSEMBLY_NAME,
                    "source_assembly_accession": SOURCE_ASSEMBLY_ACCESSION,
                    "source_coordinate_raw": coord["raw"],
                    "source_seqname": coord["seqname"],
                    "source_start": coord["start"],
                    "source_end": coord["end"],
                    "source_coordinate_convention": (
                        "reported genomic interval; "
                        "base convention not explicitly asserted here"
                    ),
                    "p1_coordinates_raw": info["p1_raw"],
                    "p10_coordinates_raw": info["p10_raw"],
                    "variant_profile": profile["variant_profile"],
                    "expression_peak_distance": (
                        profile["expression_peak_distance"]
                    ),
                    "chromatin_state": profile["chromatin_state"],
                    "genic_context": profile["genic_context"],
                    "table1_row": info["table1_row"],
                    "table4_row": profile["table4_row"],
                    "target_assembly_name": TARGET_ASSEMBLY_NAME,
                    "target_assembly_accession": TARGET_ASSEMBLY_ACCESSION,
                    "target_seqname": "",
                    "target_start": "",
                    "target_end": "",
                    "mapping_status": "unresolved",
                }
            )

    # ---------------------------------------------------------
    # Detect repeated loci and stable/unstable conflicts.
    # ---------------------------------------------------------

    grouped = defaultdict(list)

    for row in observations:
        key = (
            row["source_seqname"],
            row["source_start"],
            row["source_end"],
        )
        grouped[key].append(row)

    loci = []

    for index, (key, rows) in enumerate(
        sorted(grouped.items()),
        start=1,
    ):
        seqname, start, end = key

        stabilities = sorted(
            {row["stability"] for row in rows}
        )

        cell_lines = sorted(
            {row["cell_line"] for row in rows}
        )

        if stabilities == ["Stable"]:
            locus_status = "stable_only"
            gold_binary_eligible = "TRUE"

        elif stabilities == ["Unstable"]:
            locus_status = "unstable_only"
            gold_binary_eligible = "TRUE"

        else:
            locus_status = "context_conflict"
            gold_binary_eligible = "FALSE"

        loci.append(
            {
                "source_locus_id": f"DHIMAN2020_L{index:03d}",
                "source_assembly_accession": SOURCE_ASSEMBLY_ACCESSION,
                "source_seqname": seqname,
                "source_start": start,
                "source_end": end,
                "n_observations": len(rows),
                "cell_lines": ";".join(cell_lines),
                "stabilities": ";".join(stabilities),
                "locus_status": locus_status,
                "gold_binary_eligible": gold_binary_eligible,
            }
        )

    if len(loci) != 13:
        raise SystemExit(
            f"ERROR: expected 13 unique source loci, found {len(loci)}"
        )

    n_conflicts = sum(
        row["locus_status"] == "context_conflict"
        for row in loci
    )

    if n_conflicts != 3:
        raise SystemExit(
            f"ERROR: expected 3 context-conflicted loci, "
            f"found {n_conflicts}"
        )

    # ---------------------------------------------------------
    # Output
    # ---------------------------------------------------------

    observation_columns = [
        "observation_id",
        "study_id",
        "cell_line",
        "copy_number",
        "stability",
        "benchmark_role_observation",
        "evidence_class",
        "selected_passage",
        "source_assembly_name",
        "source_assembly_accession",
        "source_coordinate_raw",
        "source_seqname",
        "source_start",
        "source_end",
        "source_coordinate_convention",
        "p1_coordinates_raw",
        "p10_coordinates_raw",
        "variant_profile",
        "expression_peak_distance",
        "chromatin_state",
        "genic_context",
        "table1_row",
        "table4_row",
        "target_assembly_name",
        "target_assembly_accession",
        "target_seqname",
        "target_start",
        "target_end",
        "mapping_status",
    ]

    locus_columns = [
        "source_locus_id",
        "source_assembly_accession",
        "source_seqname",
        "source_start",
        "source_end",
        "n_observations",
        "cell_lines",
        "stabilities",
        "locus_status",
        "gold_binary_eligible",
    ]

    gaidukov_columns = [
        "table4_row",
        "cell_line",
        "copy_number",
        "stability",
        "variant_profile",
        "expression_peak_distance",
        "chromatin_state",
        "genic_context",
    ]

    write_tsv(
        outdir / "dhiman_observed_sites.tsv",
        observations,
        observation_columns,
    )

    write_tsv(
        outdir / "dhiman_unique_loci.tsv",
        loci,
        locus_columns,
    )

    write_tsv(
        outdir / "gaidukov_profiles_from_dhiman.tsv",
        gaidukov_profiles,
        gaidukov_columns,
    )

    stable_obs = sum(
        row["stability"] == "Stable"
        for row in observations
    )

    unstable_obs = sum(
        row["stability"] == "Unstable"
        for row in observations
    )

    stable_loci = sum(
        row["locus_status"] == "stable_only"
        for row in loci
    )

    unstable_loci = sum(
        row["locus_status"] == "unstable_only"
        for row in loci
    )

    print("DHIMAN 2020 EXTRACTION: PASS")
    print(f"Phenotyped cell lines:            {len(cells)}")
    print(f"Integration-site observations:    {len(observations)}")
    print(f"  Stable observations:            {stable_obs}")
    print(f"  Unstable observations:          {unstable_obs}")
    print(f"Unique source loci:               {len(loci)}")
    print(f"  Stable-only loci:               {stable_loci}")
    print(f"  Unstable-only loci:             {unstable_loci}")
    print(f"  Context-conflicted loci:        {n_conflicts}")
    print(
        "Binary gold-eligible loci:       "
        f"{stable_loci + unstable_loci}"
    )
    print(f"Gaidukov profiles captured:       {len(gaidukov_profiles)}")
    print(f"Source assembly:                  {SOURCE_ASSEMBLY_ACCESSION}")
    print(f"Target assembly:                  {TARGET_ASSEMBLY_ACCESSION}")


if __name__ == "__main__":
    main()
