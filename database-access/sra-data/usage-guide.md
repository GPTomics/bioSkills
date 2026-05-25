# SRA Data Usage Guide

## Overview

Download raw sequencing reads from NCBI SRA. Encodes the source-of-truth decision (SRA-direct vs ENA mirror vs AWS/GCP STRIDES cloud), prefetch `--max-size` trap (silent 20 GB skip), fasterq-dump uncompressed-scratch trap (~3x final size), `--include-technical` for 10x single-cell records, MD5 validation, accession hierarchy navigation (SRR/SRX/SRS/SRP/PRJNA), pysradb for metadata, and Aspera deprecation post-2019.

## Prerequisites

```bash
conda install -c bioconda sra-tools         # 3.0+
pip install pysradb                         # optional, for metadata
fasterq-dump --version
```

For AWS STRIDES, install `aws-cli` and run from EC2 in `us-east-1` for free egress.

## Quick Start

- "Download SRR12345678 as paired-end FASTQ via the ENA mirror with MD5 verification"
- "Resolve GSE123456 to its SRR run accessions using pysradb, then download with prefetch"
- "Download a 10x single-cell record -- include technical reads (barcodes/UMIs/indexes)"
- "Pull SRA data from AWS Open Data (STRIDES) into an EC2 instance in us-east-1"
- "Validate downloaded FASTQ against ENA's md5; fail loudly on mismatch"

## Example Prompts

### Default: ENA mirror

> "Download SRR12345678 via the ENA mirror (https://ftp.sra.ebi.ac.uk/). It's typically faster than SRA-direct because the FASTQ is pre-compressed -- no SRA->FASTQ conversion needed. Get the md5 from the ENA portal API and verify."

### Resolving GSE to SRR

> "I have GSE123456 from a GEO paper. Use pysradb gse_to_srp -> srp_to_srr to get all SRR run accessions. Then batch-download them via ENA."

### Large prefetch without silent skip

> "Use prefetch with --max-size 200G explicitly -- the default 20 GB silently skips larger runs. After prefetch, run vdb-validate and then fasterq-dump."

### 10x single-cell

> "This is a 10x v3 record. fasterq-dump with --include-technical --split-files. Expect 3 files: R1 (28-bp barcode+UMI), R2 (cDNA), I1 (sample index). The default (without --include-technical) only gives R2 which is useless for CellRanger."

### Cloud-native pipeline

> "I'm running a Nextflow pipeline on AWS Batch in us-east-1. Pull SRA data from s3://sra-pub-run-odp/ with --no-sign-request. Same-region transfer is free."

### Scratch disk management

> "fasterq-dump writes uncompressed FASTQ to scratch (~3x final compressed size). My scratch dir has 500 GB free; the run is 200 GB compressed. That's tight -- use fastq-dump --gzip instead, which writes compressed in-place."

## What the Agent Will Do

1. Choose download source: ENA (default), STRIDES (in-cloud), SRA-direct (fallback).
2. Use prefetch with explicit `--max-size` (never the 20 GB default for unknown sizes).
3. Run `vdb-validate` (SRA format) or `md5sum -c` (ENA FASTQ) on every downloaded file.
4. For 10x or single-cell records, always pass `--include-technical`.
5. After fasterq-dump, compress with pigz post-hoc (fasterq-dump doesn't compress).
6. Resolve hierarchical accessions (PRJNA, GSE, SRX -> SRR) via pysradb or ENA portal API.
7. Recommend STRIDES for any in-cloud analysis pipeline.
8. Warn on SRA-direct during US business hours (9 AM-5 PM ET); recommend ENA or off-peak.

## Tips

- ENA mirror is the right default for off-cloud downloads in 2026 -- typically faster and gives FASTQ directly.
- prefetch's `--max-size 20G` default silently skips larger runs. Set it generously (200G) or use `pysradb metadata` to check sizes first.
- fasterq-dump scratch overhead is ~3x final size; on tight scratch use fastq-dump --gzip (slower but lower scratch).
- 10x records require `--include-technical` to get barcode/UMI reads; without it CellRanger and STARsolo will fail.
- AWS STRIDES (`s3://sra-pub-run-odp/sra/`) is free egress in same-region; cross-region pulls incur cost.
- Aspera (`ascp`): NCBI public access retired 2019; ENA public Aspera retired ~2023; institutional accounts still work.
- For finding accessions: `efetch -db sra -id <UID> -rettype runinfo` gives a CSV; pysradb's `metadata` is more ergonomic.
- vdb-config persistence: in Docker, mount `~/.ncbi/` as a volume so user-settings.mkfg survives container rebuilds.
- For raw read downloads do NOT use NCBI Datasets CLI -- that's for genome assemblies, not sequencing reads.

## Related Skills

- entrez-search - Search SRA db for accessions
- geo-data - GEO Series often link to SRA via ELink
- read-qc/quality-reports - QC the downloaded FASTQ
- read-qc/fastp-workflow - Adapter trim downloaded FASTQ
- read-alignment/bwa-alignment - Align downloaded reads
- ncbi-datasets-cli - Modern bulk path for genomes (NOT for raw reads)
