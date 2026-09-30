#!/usr/bin/env python3
"""Load API credentials from .env file.

Never prints or logs the values. Other modules import getenv() which
returns the value or None if not set / file missing.
"""
import os
from pathlib import Path

_ENV_FILE = Path(__file__).parent.parent / ".env"

_loaded = False


def _load():
    global _loaded
    if _loaded:
        return
    _loaded = True
    if _ENV_FILE.exists():
        for line in _ENV_FILE.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, val = line.partition("=")
                key = key.strip()
                val = val.strip().strip("'").strip('"')
                if key and key not in os.environ:
                    os.environ[key] = val


def getenv(name: str, default: str | None = None) -> str | None:
    """Return a credential from env or .env file. Never prints the value."""
    _load()
    return os.environ.get(name, default)


def has_key(name: str) -> bool:
    _load()
    return bool(os.environ.get(name))


# Convenience names (values never exposed)
def ncbi_api_key() -> str | None:
    return getenv("NCBI_API_KEY")


def synapse_token() -> str | None:
    return getenv("SYNAPSE_TOKEN")
