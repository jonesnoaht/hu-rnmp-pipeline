#!/usr/bin/env python3
"""
SRA data fetcher using NCBI E-utilities.

Uses NCBI_API_KEY from .env to raise rate limit from 3→10 req/s.
Never prints the API key value.

Usage:
    python fetch_sra.py --accession SRP123456 --outdir /data/raw
"""
import argparse
import os
import subprocess
import sys
import time
import urllib.parse
import urllib.request

# Local imports
sys.path.insert(0, os.path.dirname(__file__))
from credentials import ncbi_api_key, has_key

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _api_key_param() -> str:
    """Return &api_key=... if key is set, else empty string."""
    key = ncbi_api_key()
    if key:
        return f"&api_key={key}"
    return ""


def _rate_limit_sleep():
    """Sleep to respect rate limits."""
    if has_key("NCBI_API_KEY"):
        time.sleep(0.11)   # 10 req/s
    else:
        time.sleep(0.34)   # 3 req/s


def esearch(db: str, term: str, retmax: int = 50) -> list[str]:
    """Search NCBI, return list of IDs."""
    url = (f"{EUTILS}/esearch.fcgi?db={db}"
           f"&term={urllib.parse.quote(term)}"
           f"&retmax={retmax}&retmode=json{_api_key_param()}")
    _rate_limit_sleep()
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.loads(resp.read())
    return data.get("esearchresult", {}).get("idlist", [])


def efetch_sra(accession: str, outdir: str):
    """Prefetch + fasterq-dump an SRA accession."""
    os.makedirs(outdir, exist_ok=True)
    # prefetch
    subprocess.run([
        "prefetch", accession, "-o", os.path.join(outdir, accession)
    ], check=True, capture_output=True)
    # fasterq-dump
    subprocess.run([
        "fasterq-dump", accession,
        "--outdir", outdir,
        "--threads", str(os.cpu_count() or 4),
        "--split-files"
    ], check=True, capture_output=True)
    # gzip
    for f in os.listdir(outdir):
        if f.endswith(".fastq"):
            subprocess.run(["gzip", "-f", os.path.join(outdir, f)], check=True)


def search_hydroxyurea(db: str = "sra", retmax: int = 50) -> list[str]:
    """Search for hydroxyurea studies."""
    return esearch(db, "hydroxyurea", retmax)


def search_rnase_h2(db: str = "sra", retmax: int = 50) -> list[str]:
    """Search for RNase H2 / RNASEH2 studies."""
    terms = ["RNASEH2A", "RNASEH2B", "RNASEH2C", "RNase H2"]
    results = []
    for term in terms:
        results.extend(esearch(db, term, retmax))
    return list(set(results))


def main():
    import json  # noqa: E402

    ap = argparse.ArgumentParser(description="Fetch SRA data with API key")
    ap.add_argument("--accession", required=True, help="SRA accession (SRP/SRR)")
    ap.add_argument("--outdir", default="./data/raw", help="Output directory")
    args = ap.parse_args()

    key_status = "10 req/s (API key set)" if has_key("NCBI_API_KEY") else "3 req/s (no key)"
    print(f"NCBI rate limit: {key_status}", file=sys.stderr)

    efetch_sra(args.accession, args.outdir)
    print(f"Done: {args.accession} → {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    import json  # noqa: F401,E402
    main()
