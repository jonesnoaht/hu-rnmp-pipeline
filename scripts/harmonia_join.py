#!/usr/bin/env python3
"""Run Harmonia when the binary is on PATH, otherwise write a skip record.

The pipeline image does not ship the Harmonia executable. The stage still
emits the two files the report expects so Nextflow can finish.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description="Harmonia join, or a recorded skip")
    parser.add_argument("--enrichment", required=True)
    parser.add_argument("--rules", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--complement", required=True)
    parser.add_argument(
        "--run",
        default="false",
        help="true attempts to execute harmonia; anything else skips",
    )
    args = parser.parse_args()

    requested = str(args.run).lower() == "true"
    binary = shutil.which("harmonia")
    if requested and binary:
        subprocess.run(
            [
                binary,
                "harmonize",
                args.enrichment,
                "--rules",
                args.rules,
                "--output",
                args.output,
                "--complement",
                args.complement,
            ],
            check=True,
        )
        return

    if not requested:
        reason = (
            "Harmonia was not run because --run is false. "
            "The harmonia binary is not installed in the pipeline image."
        )
    else:
        reason = "Harmonia was not run because the harmonia binary was not found on PATH."

    with open(args.output, "w") as handle:
        handle.write("status\tskipped\n")
        handle.write(f"reason\t{reason}\n")
    with open(args.complement, "w") as handle:
        json.dump({"status": "skipped", "reason": reason}, handle, indent=2)
        handle.write("\n")
    print(reason, file=sys.stderr)


if __name__ == "__main__":
    main()
