#!/usr/bin/env python3

import argparse
import hashlib
import shutil
import urllib.parse
import urllib.request
from pathlib import Path


def sha256sum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_fasta(accession):
    params = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": accession,
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}"

    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nuclear-fasta", required=True)
    parser.add_argument("--mitochondrial-accession", required=True)
    parser.add_argument("--mitochondrial-fasta", required=True)
    parser.add_argument("--mapping-fasta", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()

    nuclear_fasta = Path(args.nuclear_fasta)
    mitochondrial_fasta = Path(args.mitochondrial_fasta)
    mapping_fasta = Path(args.mapping_fasta)
    metadata_file = Path(args.metadata)
    sha256_file = Path(args.sha256)

    for path in [
        mitochondrial_fasta,
        mapping_fasta,
        metadata_file,
        sha256_file,
    ]:
        path.parent.mkdir(parents=True, exist_ok=True)

    fasta_text = fetch_fasta(args.mitochondrial_accession)

    if not fasta_text.startswith(">"):
        raise RuntimeError("NCBI EFetch did not return a valid FASTA record.")

    header = fasta_text.splitlines()[0]

    if args.mitochondrial_accession not in header:
        raise RuntimeError(
            f"Expected accession {args.mitochondrial_accession} "
            f"was not found in FASTA header: {header}"
        )

    mitochondrial_fasta.write_text(fasta_text, encoding="utf-8")

    with open(mapping_fasta, "wb") as output:
        with open(nuclear_fasta, "rb") as nuclear:
            output.write(nuclear.read())

        with open(nuclear_fasta, "rb") as nuclear:
            nuclear.seek(-1, 2)
            if nuclear.read(1) != b"\n":
                output.write(b"\n")

        with open(mitochondrial_fasta, "rb") as mitochondrial:
            output.write(mitochondrial.read())

    metadata_file.write_text(
        "\t".join(
            [
                "resource",
                "accession",
                "source",
                "purpose",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "mitochondrial_sequence",
                args.mitochondrial_accession,
                "NCBI_EUtils_EFetch",
                "ATAC-seq mitochondrial mapping and QC",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with open(sha256_file, "w", encoding="utf-8") as handle:
        for path in [nuclear_fasta, mitochondrial_fasta, mapping_fasta]:
            handle.write(f"{sha256sum(path)}  {path}\n")


if __name__ == "__main__":
    main()
