#!/usr/bin/env python3

import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import PurePosixPath

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
OFFICE_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"

NS = {"m": MAIN}

KEYWORDS = (
    "integration",
    "tla",
    "stable",
    "unstable",
    "gaidukov",
    "chromosome",
    "scaffold",
    "position",
    "coordinate",
    "site",
    "locus",
)

def all_text(node):
    return "".join(
        t.text or ""
        for t in node.iter()
        if t.tag == f"{{{MAIN}}}t"
    )

def shared_strings(zf):
    path = "xl/sharedStrings.xml"
    if path not in zf.namelist():
        return []

    root = ET.fromstring(zf.read(path))
    return [all_text(si) for si in root.findall(f"{{{MAIN}}}si")]

def relationships(zf):
    root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    out = {}

    for rel in root.findall(f"{{{PKG_REL}}}Relationship"):
        out[rel.attrib["Id"]] = rel.attrib["Target"]

    return out

def sheets(zf):
    root = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = relationships(zf)

    for sheet in root.findall("m:sheets/m:sheet", NS):
        name = sheet.attrib["name"]
        rid = sheet.attrib[f"{{{OFFICE_REL}}}id"]
        target = rels[rid]

        if target.startswith("/"):
            path = target.lstrip("/")
        else:
            path = str(PurePosixPath("xl") / target)

        yield name, path

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
        try:
            return shared[int(value)]
        except Exception:
            return value

    return value

def inspect_sheet(zf, name, path, shared):
    root = ET.fromstring(zf.read(path))

    dimension = root.find("m:dimension", NS)
    dim = dimension.attrib.get("ref", "?") if dimension is not None else "?"

    print()
    print("=" * 78)
    print(f"SHEET: {name}")
    print(f"DIMENSION: {dim}")
    print("=" * 78)

    rows = root.findall(".//m:sheetData/m:row", NS)

    preview_count = 0
    matches = []

    for row in rows:
        values = []

        for cell in row.findall("m:c", NS):
            ref = cell.attrib.get("r", "?")
            val = cell_value(cell, shared).strip()

            if val:
                values.append((ref, val))

        if not values:
            continue

        if preview_count < 6:
            text = " | ".join(f"{ref}={val}" for ref, val in values[:20])
            print(f"PREVIEW row {row.attrib.get('r', '?')}: {text}")
            preview_count += 1

        joined = " ".join(v.lower() for _, v in values)

        if any(keyword in joined for keyword in KEYWORDS):
            matches.append(
                (
                    row.attrib.get("r", "?"),
                    " | ".join(
                        f"{ref}={val}"
                        for ref, val in values[:30]
                    ),
                )
            )

    if matches:
        print()
        print("KEYWORD MATCHES:")
        for rownum, text in matches[:120]:
            print(f"  row {rownum}: {text}")

def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: inspect_xlsx.py FILE.xlsx")

    path = sys.argv[1]

    with zipfile.ZipFile(path) as zf:
        shared = shared_strings(zf)
        workbook_sheets = list(sheets(zf))

        print(f"Workbook: {path}")
        print(f"Number of sheets: {len(workbook_sheets)}")
        print("Sheet names:")

        for name, _ in workbook_sheets:
            print(f"  - {name}")

        for name, sheet_path in workbook_sheets:
            inspect_sheet(
                zf,
                name,
                sheet_path,
                shared,
            )

if __name__ == "__main__":
    main()
