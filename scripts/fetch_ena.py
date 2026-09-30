#!/usr/bin/env python3
"""Download public FASTQ from the ENA filereport API.

Study accessions are not expanded into every run: a multi-run report is an
error unless ``--max-runs 1`` (the smallest FASTQ). GEO series (GSE) are
rejected. Optional reservoir subsampling uses a fixed seed so a test profile
can keep a reproducible subset without fetching a private file.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import shutil
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = (
    "run_accession,sample_title,library_layout,scientific_name,"
    "fastq_ftp,fastq_md5,fastq_bytes,read_count,study_accession,study_title"
)
USER_AGENT = "hu-rnmp-pipeline/test"


def parse_filereport(text: str) -> list[dict[str, str]]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("ENA filereport was empty")
    header = lines[0].split("\t")
    if "run_accession" not in header:
        raise ValueError(f"ENA filereport has no run_accession column: {lines[0][:300]}")
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < len(header):
            parts = parts + [""] * (len(header) - len(parts))
        rows.append(dict(zip(header, parts)))
    return rows


def ftp_to_https(path: str) -> str:
    value = path.strip()
    if value.startswith("ftp://"):
        value = value[len("ftp://") :]
    elif value.startswith("http://"):
        return "https://" + value[len("http://") :]
    if value.startswith("https://"):
        return value
    return "https://" + value


def fastq_entries(row: dict[str, str]) -> list[tuple[str, str]]:
    """Return ``(https_url, md5)`` for each FASTQ listed on a filereport row."""
    ftps = [part for part in row.get("fastq_ftp", "").split(";") if part.strip()]
    md5s = [part.strip() for part in row.get("fastq_md5", "").split(";") if part.strip()]
    while len(md5s) < len(ftps):
        md5s.append("")
    return [(ftp_to_https(ftp), md5) for ftp, md5 in zip(ftps, md5s)]


def select_run(rows: list[dict[str, str]], accession: str, max_runs: int) -> dict[str, str]:
    if accession.upper().startswith("GSE"):
        raise ValueError(
            f"{accession} is a GEO series, not a sequencing run. "
            "Pass an SRR/ERR run (see samples_test.csv)."
        )
    if not rows:
        raise ValueError(f"ENA returned no read runs for {accession}")
    if len(rows) == 1:
        return rows[0]

    def _bytes(row: dict[str, str]) -> int:
        try:
            return int(row.get("fastq_bytes") or 0)
        except ValueError:
            return 0

    ordered = sorted(rows, key=_bytes)
    example = ordered[0].get("run_accession", "")
    if max_runs != 1:
        raise ValueError(
            f"{accession} resolves to {len(rows)} runs. Refusing to download the whole study. "
            f"Pass a single run accession (smallest FASTQ is {example}) "
            "or --max-runs 1 to take only that smallest run."
        )
    return ordered[0]


def fetch_filereport(accession: str, timeout: int = 60) -> str:
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ENA_FIELDS,
            "format": "tsv",
        }
    )
    request = urllib.request.Request(
        f"{ENA_FILEREPORT}?{query}",
        headers={"User-Agent": USER_AGENT},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode()


def download(url: str, dest: Path, retries: int = 3, timeout: int = 120) -> None:
    last_error: Exception | None = None
    for _attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=timeout) as response, dest.open("wb") as handle:
                shutil.copyfileobj(response, handle)
            return
        except (OSError, urllib.error.URLError) as exc:
            last_error = exc
            if dest.exists():
                dest.unlink()
    raise RuntimeError(f"download failed for {url}: {last_error}") from last_error


def md5_file(path: Path) -> str:
    digest = hashlib.md5()  # ENA publishes MD5 checksums for these files.
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_fastq(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if record[0] == "":
                return
            if any(part == "" for part in record):
                raise ValueError(f"truncated FASTQ record in {path}")
            yield record


def subsample_fastq(src: Path, dest: Path, max_reads: int, seed: int) -> int:
    """Copy ``src`` or write a reproducible reservoir sample. Returns reads written."""
    if max_reads <= 0:
        shutil.copyfile(src, dest)
        return sum(1 for _ in iter_fastq(dest))

    rng = random.Random(seed)
    reservoir: list[list[str]] = []
    for seen, record in enumerate(iter_fastq(src), start=1):
        if len(reservoir) < max_reads:
            reservoir.append(record)
        else:
            pick = rng.randrange(seen)
            if pick < max_reads:
                reservoir[pick] = record
    with gzip.open(dest, "wt") as handle:
        for record in reservoir:
            handle.writelines(record)
    return len(reservoir)


def subsample_pair(
    read1: Path,
    read2: Path,
    out1: Path,
    out2: Path,
    max_reads: int,
    seed: int,
) -> int:
    if max_reads <= 0:
        shutil.copyfile(read1, out1)
        shutil.copyfile(read2, out2)
        return sum(1 for _ in iter_fastq(out1))

    rng = random.Random(seed)
    reservoir: list[tuple[list[str], list[str]]] = []
    paired = zip(iter_fastq(read1), iter_fastq(read2), strict=True)
    for seen, (left, right) in enumerate(paired, start=1):
        item = (left, right)
        if len(reservoir) < max_reads:
            reservoir.append(item)
        else:
            pick = rng.randrange(seen)
            if pick < max_reads:
                reservoir[pick] = item
    with gzip.open(out1, "wt") as handle1, gzip.open(out2, "wt") as handle2:
        for left, right in reservoir:
            handle1.writelines(left)
            handle2.writelines(right)
    return len(reservoir)


def _verify_md5(path: Path, expected: str) -> None:
    if not expected:
        print(f"warning: no ENA md5 for {path.name}; skipped checksum", file=sys.stderr)
        return
    actual = md5_file(path)
    if actual.lower() != expected.lower():
        raise ValueError(f"md5 mismatch for {path.name}: expected {expected}, got {actual}")


def stage_run(
    row: dict[str, str],
    outdir: Path,
    max_reads: int,
    seed: int,
) -> dict:
    entries = fastq_entries(row)
    if not entries:
        raise ValueError(f"{row.get('run_accession')} has no fastq_ftp URL at ENA")
    run = row["run_accession"]
    outdir.mkdir(parents=True, exist_ok=True)
    tmp = outdir / "_download"
    tmp.mkdir(exist_ok=True)
    downloaded: list[Path] = []
    try:
        for index, (url, digest) in enumerate(entries, start=1):
            dest = tmp / f"{run}_{index}.fastq.gz"
            download(url, dest)
            _verify_md5(dest, digest)
            downloaded.append(dest)
        if len(downloaded) == 1:
            written = subsample_fastq(downloaded[0], outdir / f"{run}_R1.fastq.gz", max_reads, seed)
            layout = "SINGLE"
        elif len(downloaded) == 2:
            written = subsample_pair(
                downloaded[0],
                downloaded[1],
                outdir / f"{run}_R1.fastq.gz",
                outdir / f"{run}_R2.fastq.gz",
                max_reads,
                seed,
            )
            layout = "PAIRED"
        else:
            raise ValueError(f"{run} has {len(downloaded)} FASTQ files; expected 1 or 2")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    manifest = {
        "run_accession": run,
        "sample_title": row.get("sample_title", ""),
        "scientific_name": row.get("scientific_name", ""),
        "study_accession": row.get("study_accession", ""),
        "study_title": row.get("study_title", ""),
        "library_layout": row.get("library_layout", layout),
        "ena_read_count": row.get("read_count", ""),
        "fastq_bytes": row.get("fastq_bytes", ""),
        "urls": [url for url, _md5 in entries],
        "md5": [digest for _url, digest in entries],
        "layout": layout,
        "max_reads": max_reads,
        "subsample_seed": seed,
        "reads_written": written,
    }
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch a public ENA/SRA run as FASTQ")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--max-reads", type=int, default=0, help="0 keeps every read")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--max-runs", type=int, default=0)
    parser.add_argument("--manifest", default=None)
    args = parser.parse_args()

    rows = parse_filereport(fetch_filereport(args.accession))
    row = select_run(rows, args.accession, args.max_runs)
    manifest = stage_run(row, Path(args.outdir), args.max_reads, args.seed)
    text = json.dumps(manifest, indent=2)
    manifest_path = Path(args.manifest) if args.manifest else Path(args.outdir) / "fetch_manifest.json"
    manifest_path.write_text(text + "\n")
    print(text)
    print(
        f"Fetched {manifest['run_accession']} ({manifest['layout']}, "
        f"{manifest['reads_written']} reads written) → {args.outdir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
