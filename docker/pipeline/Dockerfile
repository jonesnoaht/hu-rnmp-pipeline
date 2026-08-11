FROM python:3.11-slim AS base

LABEL maintainer="Noah T. Jones"
LABEL description="HU-rNMP pipeline: hydroxyurea + ribonucleotide incorporation analysis"

# System deps for NGS tools + pysam
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ make zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    libssl-dev samtools bcftools bedtools \
    && rm -rf /var/lib/apt/lists/*

# Python deps
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# Pipeline scripts
COPY scripts/ /app/scripts/
COPY pathways/ /app/pathways/
COPY harmonia/ /app/harmonia/

# .env is NEVER copied into the image — mount at runtime
# Credentials loaded from: /app/.env (mounted) or env vars

WORKDIR /app
ENV PYTHONPATH=/app
ENV PATH="/app/scripts:$PATH"
