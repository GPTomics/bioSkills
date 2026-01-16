# Abundance Estimation with Bracken - Usage Guide

## Overview

Bracken (Bayesian Reestimation of Abundance with KrakEN) improves Kraken2 abundance estimates by redistributing reads assigned to higher taxonomic levels down to species. It uses a Bayesian approach based on the database's expected k-mer distribution.

## When to Use This Skill

- You have Kraken2 classification output
- You want species-level abundance estimates
- You need quantitative (not just presence/absence) results
- You're comparing abundances across samples

## Installation

```bash
conda install -c bioconda bracken
```

## Requirements

- Kraken2 database with Bracken database files
- Kraken2 report file (not the per-read output)
- Read length must match a Bracken database

## Basic Workflow

```bash
# 1. Run Kraken2 (if not already done)
kraken2 --db kraken_db --report report.txt --paired R1.fq.gz R2.fq.gz

# 2. Run Bracken
bracken -d kraken_db -i report.txt -o abundances.txt -r 150 -l S
```

## Understanding Bracken Output

| Column | Description |
|--------|-------------|
| name | Taxon name |
| taxonomy_id | NCBI taxon ID |
| taxonomy_lvl | Level (S, G, etc.) |
| kraken_assigned_reads | Direct Kraken2 assignments |
| added_reads | Reads redistributed from higher levels |
| new_est_reads | Total estimated reads |
| fraction_total_reads | Proportion of classified reads |

## Key Insight

Kraken2 assigns many reads to genus or higher when it can't distinguish species. Bracken probabilistically redistributes these reads based on which species are present and their expected k-mer composition.

## Common Issues

### Missing Bracken Database

```bash
# Build Bracken database for your read length
bracken-build -d /path/to/kraken_db -t 8 -l 150
```

### Read Length Mismatch

Use the closest available read length. Common options: 100, 150, 200, 250, 300.

## Resources

- [Bracken GitHub](https://github.com/jenniferlu717/Bracken)
- [Bracken Paper](https://peerj.com/articles/cs-104/)
