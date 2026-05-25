# NCBI Datasets CLI Usage Guide

## Overview

Use the NCBI Datasets v2 CLI (launched 2023, current as of 2024) for genome and gene-centric bulk workflows. Encodes the defection rule (use Datasets for genome/gene; E-utilities for PubMed/SRA/custom queries), the `--dehydrated` flag for cloud-friendly parallel pulls, JSON-lines output + `dataformat` conversion to TSV, automatic MD5 verification, and the choice between `--reference` (one per species) vs full set.

## Prerequisites

```bash
conda install -c conda-forge ncbi-datasets-cli
datasets --version    # 16.0+
dataformat --version
```

For high-volume bulk: get an NCBI API key from `https://www.ncbi.nlm.nih.gov/account/settings/` and pass `--api-key`.

## Quick Start

- "Download human reference genome (GCF_000001405.40) with genome + GFF3 + protein + CDS"
- "Pull every reference bacterial genome from RefSeq; use --dehydrated + aria2c for parallel transfer"
- "Get a TSV of all reference E. coli assemblies with N50 and assembly level"
- "Find NCBI ortholog set for human BRCA1 across mammals (one rep per species)"
- "Download SARS-CoV-2 genome assemblies released after 2024-01-01"

## Example Prompts

### Single-assembly download

> "Download the human GRCh38 reference assembly (GCF_000001405.40) via Datasets CLI. Include genome,gff3,gtf,protein,cds. Datasets verifies MD5 automatically -- no manual checksum step needed."

### Bulk download via --dehydrated

> "Pull every RefSeq reference bacterial genome with annotation. The serial path takes hours; use --dehydrated to get a fetch.txt of URLs, then aria2c --max-concurrent-downloads=8 for parallel transfer."

### Bulk metadata via dataformat tsv

> "Get a TSV of all reference Salmonella enterica assemblies released since 2024-01-01 with: accession, organism-name, assembly-level, scaffold-n50, submission-date. Use datasets summary + dataformat tsv genome --fields=..."

### NCBI ortholog set

> "Get NCBI's ortholog set for human BRCA1 across Mammalia. Use datasets summary gene symbol BRCA1 --taxon human --ortholog. Note: this is the simple single-representative-per-species view; for tree-reconciled orthology with co-orthologs use Ensembl Compara via ortholog-inference."

### Datasets vs E-utilities decision

> "I need 100 reference genomes. Don't loop EFetch -- use datasets download genome accession ... It's 5-50x faster, handles checksums, and parallelizes within one ZIP. For PubMed or SRA reads, stay with E-utilities -- Datasets doesn't cover those."

## What the Agent Will Do

1. Determine if the question is in Datasets' scope (genome/gene/ortholog/virus) -- if not, defect to E-utilities.
2. Pick `summary` (metadata only, JSON-lines) vs `download` (data files).
3. For metadata, pipe through `dataformat tsv <type> --fields=...` for tabular.
4. For bulk downloads >100 genomes, use `--dehydrated` + parallel transfer with aria2c.
5. Use `--reference` only when one-per-species is desired; drop for full set.
6. Choose `--include` to limit data types (genome, gff3, protein, cds, etc.).
7. Pass `--api-key` for high-volume metadata calls.
8. Document the version of Datasets CLI in pipelines (v16+ as of 2024).

## Tips

- The Datasets CLI is the supported modern path for genome and gene data. assembly_summary.txt scraping is the pre-2023 approach -- works but fragile.
- Datasets verifies MD5 automatically on download; no manual md5sum step needed.
- `--dehydrated` separates discovery from transfer: download metadata stubs, inspect, then parallel-pull. Essential for HPC / cloud.
- `--reference` returns one canonical assembly per species. Drop the flag for full set; combine with `--assembly-level chromosome,complete` for quality filter.
- `dataformat tsv <type> --help` lists valid field names per summary type. Don't guess field names.
- Subcommand and flag names changed between major versions; pin to v16+ for current docs.
- Datasets covers genome/gene/virus/ortholog. PubMed, SRA reads, BLAST, and custom Entrez queries are out of scope -- use E-utilities or other dedicated skills.
- For tree-aware orthology with co-orthologs, Datasets is too simple -- use Ensembl Compara via `ortholog-inference` or de novo via `comparative-genomics/ortholog-inference`.

## Related Skills

- entrez-search - PubMed and non-genome queries
- entrez-fetch - Single-record fetches outside genome/gene scope
- batch-downloads - Bulk E-utilities for non-genome data
- sra-data - Raw sequencing reads (NOT covered by Datasets)
- ensembl-rest - Ensembl REST as alternative
- ortholog-inference - Compara/OMA/OrthoDB for tree-aware orthology
