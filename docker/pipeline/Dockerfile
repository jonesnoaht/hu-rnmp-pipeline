FROM python:3.11-slim AS base

LABEL maintainer="Noah T. Jones"
LABEL description="HU-rNMP pipeline: hydroxyurea + ribonucleotide incorporation analysis"

# NGS tools for fetch, QC, end-to-end alignment, and the Python callers.
# Build context is the repository root (see .github/workflows/ci.yml).
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ make zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    libssl-dev \
    samtools bcftools bedtools bowtie2 fastqc fastp \
    curl ca-certificates gzip bzip2 procps \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY scripts/ /app/scripts/
COPY pathways/ /app/pathways/
COPY harmonia/ /app/harmonia/

WORKDIR /app
ENV PYTHONPATH=/app
ENV PATH="/app/scripts:${PATH}"
