#!/usr/bin/env python3
"""Safely preflight or explicitly submit shared-reference provisioning."""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROVISION_PATH = ROOT / "workflow/scripts/chipseq_reference_provision.py"
SPEC = importlib.util.spec_from_file_location("chipseq_reference_provision_support", PROVISION_PATH)
PROVISION = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(PROVISION)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--submit",
        action="store_true",
        help="Explicitly submit; without this flag the command is check-only",
    )
    result.add_argument(
        "--development-dirty-check",
        action="store_true",
        help="Allow a dirty worktree only for development check-only validation",
    )
    return result


def main() -> None:
    args = parser().parse_args()
    if args.submit and args.development_dirty_check:
        raise ValueError("Dirty-source override is forbidden for submission")
    report = PROVISION.preflight(allow_dirty=args.development_dirty_check)
    report["mode"] = "submit" if args.submit else "check"
    print(json.dumps(report, indent=2, sort_keys=True))
    if not args.submit:
        print("CHECK_ONLY_NO_SBATCH=PASS")
        return

    head = report["repository"]["head"]
    command = [
        "sbatch",
        "--parsable",
        f"--export=ALL,CHIP_PROVISION_EXPECTED_COMMIT={head}",
        "workflow/slurm/chipseq_provision_shared_reference.sbatch",
    ]
    job_id = subprocess.check_output(command, cwd=ROOT, text=True).strip()
    if not job_id or not job_id.split(";", 1)[0].isdigit():
        raise ValueError(f"Unexpected sbatch response: {job_id!r}")
    print(f"SUBMITTED={job_id}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error
