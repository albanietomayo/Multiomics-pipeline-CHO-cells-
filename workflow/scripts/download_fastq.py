#!/usr/bin/env python3

"""
Download and validate one FASTQ file described in a FASTQ manifest.

The file is first written to a temporary .part file. During download,
its size and MD5 checksum are calculated and compared with the values
reported by ENA. The final output is created only when both checks pass.
"""

import argparse
import hashlib
from pathlib import Path

import pandas as pd
import requests
import time

def normalize_url(url):
    """
    Convert ENA FTP-style paths to HTTPS URLs when necessary.
    """
    url = str(url).strip()

    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]

    if not url.startswith(("http://", "https://")):
        return "https://" + url

    return url


def calculate_md5(path):
    """
    Calculate the MD5 checksum of an existing file.
    """
    md5 = hashlib.md5()

    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(chunk)

    return md5.hexdigest()


def validate_existing_file(path, expected_bytes, expected_md5):
    """
    Check whether an existing file already matches the ENA metadata.
    """
    if not path.exists():
        return False

    observed_bytes = path.stat().st_size

    if observed_bytes != expected_bytes:
        return False

    observed_md5 = calculate_md5(path)

    return observed_md5.lower() == expected_md5.lower()


def download_and_validate(
    url,
    output,
    expected_bytes,
    expected_md5,
    max_attempts=5,
    retry_delay=5,
):
    """
    Download one FASTQ and validate its size and MD5 checksum.

    Failed transfers are retried from the beginning. Incomplete
    temporary files are removed before each new attempt.
    """

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    if validate_existing_file(
        output,
        expected_bytes,
        expected_md5,
    ):
        print(f"[OK] Existing validated file: {output}")
        return

    temporary = Path(str(output) + ".part")
    normalized_url = normalize_url(url)

    last_error = None

    for attempt in range(1, max_attempts + 1):

        if temporary.exists():
            temporary.unlink()

        print(
            f"[INFO] Download attempt "
            f"{attempt}/{max_attempts}: {normalized_url}"
        )
        print(f"[INFO] Output: {output}")

        md5 = hashlib.md5()
        downloaded_bytes = 0

        try:
            with requests.get(
                normalized_url,
                stream=True,
                timeout=(30, 300),
            ) as response:

                response.raise_for_status()

                with open(temporary, "wb") as handle:
                    for chunk in response.iter_content(
                        chunk_size=1024 * 1024
                    ):
                        if not chunk:
                            continue

                        handle.write(chunk)
                        md5.update(chunk)
                        downloaded_bytes += len(chunk)

            observed_md5 = md5.hexdigest()

            if downloaded_bytes != expected_bytes:
                raise ValueError(
                    f"Size mismatch: expected "
                    f"{expected_bytes} bytes, downloaded "
                    f"{downloaded_bytes} bytes"
                )

            if observed_md5.lower() != expected_md5.lower():
                raise ValueError(
                    f"MD5 mismatch: expected {expected_md5}, "
                    f"observed {observed_md5}"
                )

            temporary.replace(output)

            print(
                f"[OK] Download validated: {output} "
                f"({downloaded_bytes} bytes)"
            )

            return

        except Exception as error:

            last_error = error

            if temporary.exists():
                temporary.unlink()

            print(
                f"[WARN] Download attempt "
                f"{attempt}/{max_attempts} failed: {error}"
            )

            if attempt < max_attempts:
                wait_seconds = retry_delay * (2 ** (attempt - 1))

                print(
                    f"[INFO] Retrying in "
                    f"{wait_seconds} seconds..."
                )

                time.sleep(wait_seconds)

    raise RuntimeError(
        f"Download failed after {max_attempts} attempts: {output}"
    ) from last_error


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Download and validate one FASTQ file "
            "from a FASTQ manifest."
        )
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="FASTQ manifest TSV.",
    )

    parser.add_argument(
        "--run-accession",
        required=True,
        help="ENA/SRA run accession.",
    )

    parser.add_argument(
        "--filename",
        required=True,
        help="FASTQ filename to download.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Destination FASTQ path.",
    )

    args = parser.parse_args()

    manifest = pd.read_csv(
        args.manifest,
        sep="\t",
        dtype=str,
    ).fillna("")

    selected = manifest[
        (manifest["run_accession"] == args.run_accession)
        & (manifest["filename"] == args.filename)
    ]

    if len(selected) != 1:
        raise ValueError(
            "Expected exactly one manifest row for "
            f"{args.run_accession} / {args.filename}, "
            f"found {len(selected)}"
        )

    row = selected.iloc[0]

    download_and_validate(
        url=row["fastq_ftp"],
        output=args.output,
        expected_bytes=int(row["fastq_bytes"]),
        expected_md5=row["fastq_md5"],
    )


if __name__ == "__main__":
    main()
