#!/usr/bin/env python3

"""
Retrieve a fixed NCBI reference genome and its annotation
and generate auditable reference resources.

TFM multi-omics pipeline.
Reference accession and output paths are defined externally
in config/config.yaml.
"""

from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import shlex
import shutil
import subprocess
import sys
import tempfile
import zipfile

import yaml


# ============================================================
# Paths - Localiza la raiz del proyecto y config.yaml
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

CONFIG_FILE = PROJECT_ROOT / "config" / "config.yaml"


# ============================================================
# Configuration - Lee los parametros de config.yaml
# ============================================================

def load_config(config_file):
    """Load pipeline configuration from YAML."""

    with open(config_file, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)



# ============================================================
# Reference configuration - Valida la referencia del proyecto
# ============================================================

def get_reference_config(config):
    """Validate and return the reference genome configuration."""

    reference = config.get("reference")

    if not isinstance(reference, dict):
        raise ValueError(
            "Missing or invalid 'reference' section in config/config.yaml."
        )

    required_fields = [
        "source",
        "species",
        "taxid",
        "assembly_name",
        "refseq_accession",
        "genbank_accession",
        "annotation_release",
        "dir",
        "fasta",
        "gff3",
        "gtf",
        "sequence_report",
        "metadata",
        "sha256",
    ]

    missing_fields = [
        field
        for field in required_fields
        if field not in reference or reference[field] in (None, "")
    ]

    if missing_fields:
        raise ValueError(
            "Missing required reference configuration fields: "
            + ", ".join(missing_fields)
        )

    return reference


# ============================================================
# Reference paths - Resuelve las rutas de salida
# ============================================================

def get_reference_paths(reference):
    """Resolve reference output paths relative to the project root."""

    return {
        "dir": PROJECT_ROOT / reference["dir"],
        "fasta": PROJECT_ROOT / reference["fasta"],
        "gff3": PROJECT_ROOT / reference["gff3"],
        "gtf": PROJECT_ROOT / reference["gtf"],
        "sequence_report": PROJECT_ROOT / reference["sequence_report"],
        "metadata": PROJECT_ROOT / reference["metadata"],
        "sha256": PROJECT_ROOT / reference["sha256"],
    }


# ============================================================
# NCBI Datasets command - Construye la descarga de referencia
# ============================================================

def build_datasets_command(reference, zip_path):
    """Build the NCBI Datasets command for the fixed reference assembly."""

    return [
        "datasets",
        "download",
        "genome",
        "accession",
        reference["refseq_accession"],
        "--include",
        "genome,gff3,gtf,seq-report",
        "--filename",
        str(zip_path),
        "--no-progressbar",
    ]

# ============================================================
# Reference download - Descarga el paquete de NCBI Datasets
# ============================================================

def download_reference_package(reference, temp_path):
    """Download the NCBI reference genome data package."""

    datasets_executable = shutil.which("datasets")

    if datasets_executable is None:
        raise RuntimeError(
            "NCBI Datasets executable 'datasets' was not found in PATH."
        )

    zip_path = temp_path / "reference.zip"

    command = build_datasets_command(
        reference,
        zip_path,
    )

    print(f"[INFO] Downloading reference assembly: {reference['refseq_accession']}")
    print(f"[INFO] Command: {shlex.join(command)}")

    subprocess.run(
        command,
        check=True,
    )

    if not zip_path.is_file():
        raise RuntimeError(
            f"NCBI Datasets did not create the expected ZIP file: {zip_path}"
        )

    if zip_path.stat().st_size == 0:
        raise RuntimeError(
            f"Downloaded ZIP file is empty: {zip_path}"
        )

    print(
        f"[INFO] Reference package downloaded: "
        f"{zip_path.stat().st_size:,} bytes"
    )

    return zip_path

# ============================================================
# Package extraction - Extrae y localiza los recursos NCBI
# ============================================================

def extract_reference_package(reference, zip_path, temp_path):
    """Extract the NCBI data package and locate reference files."""

    extract_dir = temp_path / "extracted"

    with zipfile.ZipFile(zip_path, "r") as archive:
        archive.extractall(extract_dir)

    assembly_dir = (
        extract_dir
        / "ncbi_dataset"
        / "data"
        / reference["refseq_accession"]
    )

    if not assembly_dir.is_dir():
        raise RuntimeError(
            f"Expected assembly directory was not found: {assembly_dir}"
        )

    fasta_candidates = list(
        assembly_dir.glob("*_genomic.fna")
    )

    if len(fasta_candidates) != 1:
        raise RuntimeError(
            "Expected exactly one genomic FASTA file, "
            f"but found {len(fasta_candidates)}."
        )

    source_files = {
        "fasta": fasta_candidates[0],
        "gff3": assembly_dir / "genomic.gff",
        "gtf": assembly_dir / "genomic.gtf",
        "sequence_report": assembly_dir / "sequence_report.jsonl",
    }

    for label, file_path in source_files.items():

        if not file_path.is_file():
            raise RuntimeError(
                f"Expected {label} file was not found: {file_path}"
            )

        if file_path.stat().st_size == 0:
            raise RuntimeError(
                f"Expected {label} file is empty: {file_path}"
            )

        print(
            f"[INFO] Located {label}: "
            f"{file_path.name} "
            f"({file_path.stat().st_size:,} bytes)"
        )

    return source_files


# ============================================================
# NCBI checksum validation - Verifica la integridad del paquete
# ============================================================

def validate_ncbi_md5(temp_path):
    """Validate extracted NCBI package files against md5sum.txt."""

    extract_dir = temp_path / "extracted"
    md5_path = extract_dir / "md5sum.txt"

    if not md5_path.is_file():
        raise RuntimeError(
            f"NCBI MD5 manifest was not found: {md5_path}"
        )

    checked_files = 0

    with open(md5_path, "r", encoding="utf-8") as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue

            expected_md5, relative_name = line.split(maxsplit=1)

            relative_name = relative_name.lstrip("*")
            file_path = extract_dir / relative_name

            if not file_path.is_file():
                raise RuntimeError(
                    f"File listed in NCBI MD5 manifest was not found: "
                    f"{file_path}"
                )

            md5 = hashlib.md5()

            with open(file_path, "rb") as file_handle:
                for chunk in iter(
                    lambda: file_handle.read(1024 * 1024),
                    b"",
                ):
                    md5.update(chunk)

            observed_md5 = md5.hexdigest()

            if observed_md5 != expected_md5:
                raise RuntimeError(
                    f"NCBI MD5 validation failed for {relative_name}: "
                    f"expected {expected_md5}, observed {observed_md5}."
                )

            checked_files += 1

    if checked_files == 0:
        raise RuntimeError(
            "NCBI MD5 manifest did not contain any files to validate."
        )

    print(
        f"[INFO] NCBI MD5 validation passed for "
        f"{checked_files} files."
    )

# ============================================================
# Assembly validation - Verifica la identidad de la referencia
# ============================================================

def validate_assembly_report(reference, temp_path):
    """Validate NCBI assembly metadata against config.yaml."""

    report_path = (
        temp_path
        / "extracted"
        / "ncbi_dataset"
        / "data"
        / "assembly_data_report.jsonl"
    )

    if not report_path.is_file():
        raise RuntimeError(
            f"Assembly data report was not found: {report_path}"
        )

    records = []

    with open(report_path, "r", encoding="utf-8") as handle:
        for line in handle:

            line = line.strip()

            if line:
                records.append(json.loads(line))

    if len(records) != 1:
        raise RuntimeError(
            "Expected exactly one assembly record, "
            f"but found {len(records)}."
        )

    record = records[0]

    observed_accession = record.get("accession")
    observed_organism = record.get("organism", {}).get("organismName")
    observed_taxid = record.get("organism", {}).get("taxId")
    observed_assembly = record.get("assemblyInfo", {}).get("assemblyName")

    paired_assembly = (
        record
        .get("assemblyInfo", {})
        .get("pairedAssembly", {})
        .get("accession")
    )

    expected_values = {
        "RefSeq accession": (
            observed_accession,
            reference["refseq_accession"],
        ),
        "species": (
            observed_organism,
            reference["species"],
        ),
        "taxid": (
            observed_taxid,
            reference["taxid"],
        ),
        "assembly name": (
            observed_assembly,
            reference["assembly_name"],
        ),
    }

    for label, (observed, expected) in expected_values.items():

        if observed != expected:
            raise RuntimeError(
                f"Reference validation failed for {label}: "
                f"expected {expected!r}, observed {observed!r}."
            )

    if paired_assembly and paired_assembly != reference["genbank_accession"]:
        raise RuntimeError(
            "Unexpected paired GenBank accession: "
            f"expected {reference['genbank_accession']!r}, "
            f"observed {paired_assembly!r}."
        )

    print("[INFO] Assembly metadata validated successfully.")
    print(f"[INFO] Verified accession: {observed_accession}")
    print(f"[INFO] Verified organism: {observed_organism}")
    print(f"[INFO] Verified taxid: {observed_taxid}")
    print(f"[INFO] Verified assembly name: {observed_assembly}")

    if paired_assembly:
        print(f"[INFO] Verified paired GenBank accession: {paired_assembly}")

    return record


# ============================================================
# Final reference files - Instala los recursos validados
# ============================================================

def install_reference_files(source_files, paths):
    """Copy validated reference files to their final project paths."""

    labels = [
        "fasta",
        "gff3",
        "gtf",
        "sequence_report",
    ]

    temporary_outputs = {}

    for label in labels:

        source = source_files[label]
        destination = paths[label]

        temporary_destination = Path(
            str(destination) + ".part"
        )

        if temporary_destination.exists():
            temporary_destination.unlink()

        print(
            f"[INFO] Copying {label}: "
            f"{source.name} -> {destination.name}"
        )

        shutil.copy2(
            source,
            temporary_destination,
        )

        if not temporary_destination.is_file():
            raise RuntimeError(
                f"Temporary {label} output was not created: "
                f"{temporary_destination}"
            )

        if temporary_destination.stat().st_size == 0:
            raise RuntimeError(
                f"Temporary {label} output is empty: "
                f"{temporary_destination}"
            )

        temporary_outputs[label] = temporary_destination

    for label in labels:

        temporary_destination = temporary_outputs[label]
        destination = paths[label]

        temporary_destination.replace(destination)

        print(
            f"[INFO] Installed {label}: "
            f"{destination} "
            f"({destination.stat().st_size:,} bytes)"
        )


# ============================================================
# Reference metadata - Registra la procedencia de la referencia
# ============================================================

def write_reference_metadata(reference, paths):
    """Write auditable reference metadata using portable project paths."""

    metadata_path = paths["metadata"]
    temporary_path = Path(str(metadata_path) + ".part")

    if temporary_path.exists():
        temporary_path.unlink()

    retrieval_timestamp = datetime.now(
        timezone.utc
    ).isoformat()

    metadata_rows = [
        ("source", reference["source"]),
        ("species", reference["species"]),
        ("taxid", reference["taxid"]),
        ("assembly_name", reference["assembly_name"]),
        ("refseq_accession", reference["refseq_accession"]),
        ("genbank_accession", reference["genbank_accession"]),
        (
            "configured_annotation_release",
            reference["annotation_release"],
        ),
        ("retrieval_timestamp_utc", retrieval_timestamp),
        ("fasta", reference["fasta"]),
        ("gff3", reference["gff3"]),
        ("gtf", reference["gtf"]),
        ("sequence_report", reference["sequence_report"]),
    ]

    with open(
        temporary_path,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "field",
                "value",
            ]
        )

        writer.writerows(metadata_rows)

    if not temporary_path.is_file():
        raise RuntimeError(
            f"Reference metadata file was not created: {temporary_path}"
        )

    if temporary_path.stat().st_size == 0:
        raise RuntimeError(
            f"Reference metadata file is empty: {temporary_path}"
        )

    temporary_path.replace(metadata_path)

    print(f"[INFO] Reference metadata written: {metadata_path}")


# ============================================================
# SHA-256 checksums - Registra la integridad de los recursos
# ============================================================

def calculate_sha256(file_path):
    """Calculate the SHA-256 checksum of a file."""

    sha256 = hashlib.sha256()

    with open(file_path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha256.update(chunk)

    return sha256.hexdigest()


def write_sha256_manifest(paths):
    """Write SHA-256 checksums for the installed reference resources."""

    checksum_path = paths["sha256"]
    temporary_path = Path(str(checksum_path) + ".part")

    if temporary_path.exists():
        temporary_path.unlink()

    files_to_hash = [
        paths["fasta"],
        paths["gff3"],
        paths["gtf"],
        paths["sequence_report"],
        paths["metadata"],
    ]

    for file_path in files_to_hash:

        if not file_path.is_file():
            raise RuntimeError(
                f"Cannot calculate SHA-256; file not found: {file_path}"
            )

        if file_path.stat().st_size == 0:
            raise RuntimeError(
                f"Cannot calculate SHA-256; file is empty: {file_path}"
            )

    with open(
        temporary_path,
        "w",
        encoding="utf-8",
    ) as handle:

        for file_path in files_to_hash:

            checksum = calculate_sha256(file_path)

            relative_path = file_path.relative_to(PROJECT_ROOT)

            handle.write(
                f"{checksum}  {relative_path.as_posix()}\n"
            )

    if not temporary_path.is_file():
        raise RuntimeError(
            f"SHA-256 manifest was not created: {temporary_path}"
        )

    if temporary_path.stat().st_size == 0:
        raise RuntimeError(
            f"SHA-256 manifest is empty: {temporary_path}"
        )

    temporary_path.replace(checksum_path)

    print(f"[INFO] SHA-256 manifest written: {checksum_path}")


# ============================================================
# Output directory - Prepara el directorio de referencia
# ============================================================

def prepare_reference_directory(paths):
    """Create the reference output directory if it does not exist."""

    reference_dir = paths["dir"]

    reference_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(f"[INFO] Reference directory ready: {reference_dir}")


# ============================================================
# Temporary workspace - Prepara un espacio de trabajo seguro
# ============================================================

def test_temporary_workspace(reference, paths):
    """Create a temporary workspace and preview the NCBI download command."""

    with tempfile.TemporaryDirectory(
        prefix=".tmp_reference_",
        dir=paths["dir"],
    ) as temp_dir:

        temp_path = Path(temp_dir)

        if not temp_path.is_dir():
            raise RuntimeError(
                f"Temporary directory could not be created: {temp_path}"
            )

        zip_path = temp_path / "reference.zip"

        command = build_datasets_command(
            reference,
            zip_path,
        )

        print(f"[INFO] Temporary workspace created: {temp_path}")
        print("[INFO] Planned NCBI Datasets command:")
        print(f"[INFO] {shlex.join(command)}")

    print("[INFO] Temporary workspace cleaned successfully.")


# ============================================================
# Reference acquisition - Coordina la adquisicion completa
# ============================================================

def acquire_reference(reference, paths):
    """Download, validate, and install the configured reference genome."""

    with tempfile.TemporaryDirectory(
        prefix=".tmp_reference_",
        dir=paths["dir"],
    ) as temp_dir:

        temp_path = Path(temp_dir)

        print(f"[INFO] Temporary workspace created: {temp_path}")

        zip_path = download_reference_package(
            reference,
            temp_path,
        )

        source_files = extract_reference_package(
            reference,
            zip_path,
            temp_path,
        )

        validate_ncbi_md5(temp_path)

        validate_assembly_report(
            reference,
            temp_path,
        )

        install_reference_files(
            source_files,
            paths,
        )

    print("[INFO] Temporary workspace cleaned successfully.")

    write_reference_metadata(
        reference,
        paths,
    )

    write_sha256_manifest(paths)

    print("[INFO] Reference genome acquisition completed successfully.")


# ============================================================
# Main
# ============================================================

def main():
    """Validate reference genome configuration."""

    config = load_config(CONFIG_FILE)
    reference = get_reference_config(config)
    paths = get_reference_paths(reference)
    prepare_reference_directory(paths)
    acquire_reference(reference, paths)

    print("[INFO] Reference configuration validated.")
    print(f"[INFO] Species: {reference['species']}")
    print(f"[INFO] Assembly: {reference['assembly_name']}")
    print(f"[INFO] RefSeq accession: {reference['refseq_accession']}")
    print(f"[INFO] Annotation release: {reference['annotation_release']}")
    print(f"[INFO] Reference directory: {paths['dir']}")
    print(f"[INFO] FASTA output: {paths['fasta']}")
    print(f"[INFO] GFF3 output: {paths['gff3']}")
    print(f"[INFO] GTF output: {paths['gtf']}")



if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print(
            f"[ERROR] {error}",
            file=sys.stderr,
        )

        sys.exit(1)
