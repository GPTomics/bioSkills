# Functional Profiling - Usage Guide

## Overview
HUMAnN3 profiles the functional potential of metagenomic communities by quantifying gene families (UniRef90) and inferring pathway abundances (MetaCyc).

## Prerequisites
```bash
conda create -n humann -c bioconda humann
conda activate humann

# Download databases (~16 GB total)
humann_databases --download chocophlan full /db/humann
humann_databases --download uniref uniref90_diamond /db/humann
```

## Quick Start
Tell your AI agent what you want to do:
- "Profile the functional potential of my metagenome"
- "Identify metabolic pathways in my microbial community"
- "Compare pathway abundances across treatment groups"

## Example Prompts
### Basic Profiling
> "Run HUMAnN3 on sample.fastq.gz to get pathway abundances"

> "Profile functional genes in my metagenome with 8 threads"

### Multi-sample Analysis
> "Process all fastq files through HUMAnN3 and merge the results"

> "Normalize pathway abundances to relative abundance and compare groups"

### Gene Family Analysis
> "Regroup gene families to KEGG Orthologs (KO)"

> "Show which species contribute to each pathway in my sample"

### Database Selection
> "Use the smaller uniref50 database for faster processing"

> "Run HUMAnN3 with a pre-computed MetaPhlAn profile to speed things up"

## What the Agent Will Do
1. Run MetaPhlAn for taxonomic profiling (if not provided)
2. Map reads to species-specific pangenomes
3. Translate unmapped reads to protein database
4. Quantify gene families and infer pathway abundances
5. Merge and normalize results across samples if requested

## Tips
- Concatenate paired-end reads before running (HUMAnN handles both orientations)
- Pre-compute MetaPhlAn profile with `--taxonomic-profile` for speed
- High UNMAPPED rate indicates missing database coverage
- Always normalize before cross-sample comparison (use `humann_renorm_table`)
- Stratified output shows species contributions to each pathway

## Database Options

| Database | Size | Speed | Sensitivity |
|----------|------|-------|-------------|
| uniref90_diamond | 16GB | Fast | Standard |
| uniref50_diamond | 5GB | Faster | Lower |
| uniref90_ec_filtered | 0.8GB | Fastest | EC only |

## Key Parameters

| Parameter | Description |
|-----------|-------------|
| `--threads` | Number of CPUs |
| `--memory-use` | Memory limit (minimum/maximum) |
| `--taxonomic-profile` | Pre-computed MetaPhlAn profile |
| `--bypass-nucleotide-search` | Skip pangenome search |

## Resources
- [HUMAnN3 User Manual](https://huttenhower.sph.harvard.edu/humann/)
- [bioBakery Tools](https://github.com/biobakery/biobakery)
