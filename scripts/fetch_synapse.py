#!/usr/bin/env python3
"""
Synapse dataset fetcher.

Uses SYNAPSE_TOKEN from .env to authenticate with Synapse.
Requires: pip install synapseclient
Never prints the token value.

Usage:
    python fetch_synapse.py --synapse-id syn12345678 --outdir /data/raw
    python fetch_synapse.py --search "ribonucleotide" --outdir /data/raw
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from credentials import has_key, synapse_token


def get_client():
    """Create an authenticated Synapse client. Never prints token."""
    if not has_key("SYNAPSE_TOKEN"):
        print("ERROR: SYNAPSE_TOKEN not set in .env", file=sys.stderr)
        sys.exit(1)
    try:
        import synapseclient
    except ImportError:
        print("Installing synapseclient...", file=sys.stderr)
        os.system("pip install synapseclient")
        import synapseclient

    token = synapse_token()
    client = synapseclient.Synapse()
    client.login(authToken=token, silent=True)
    return client


def download(synapse_id: str, outdir: str):
    """Download a file/entity from Synapse."""
    client = get_client()
    os.makedirs(outdir, exist_ok=True)
    entity = client.get(synapse_id, downloadLocation=outdir)
    print(f"Downloaded: {entity.name} → {outdir}", file=sys.stderr)
    return entity


def search(query: str, outdir: str, limit: int = 20):
    """Search Synapse for datasets matching query."""
    client = get_client()
    results = client.apiServices(
        "/entity?query={}&limit={}".format(
            __import__("urllib.parse", fromlist=["quote"]).quote(query), limit
        )
    )
    return results


def main():
    ap = argparse.ArgumentParser(description="Fetch Synapse data with token")
    ap.add_argument("--synapse-id", help="Synapse entity ID (e.g. syn12345678)")
    ap.add_argument("--search", help="Search query instead of direct download")
    ap.add_argument("--outdir", default="./data/raw", help="Output directory")
    args = ap.parse_args()

    if not has_key("SYNAPSE_TOKEN"):
        print("ERROR: SYNAPSE_TOKEN not set in .env", file=sys.stderr)
        sys.exit(1)

    if args.synapse_id:
        download(args.synapse_id, args.outdir)
    elif args.search:
        results = search(args.search, args.outdir)
        print(f"Found {len(results)} results", file=sys.stderr)
    else:
        ap.error("Provide --synapse-id or --search")


if __name__ == "__main__":
    main()
