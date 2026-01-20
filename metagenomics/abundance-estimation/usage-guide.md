# Abundance Estimation with Bracken - Usage Guide

## Overview
Bracken improves Kraken2 abundance estimates by redistributing reads assigned to higher taxonomic levels down to species using a Bayesian approach based on expected k-mer distributions.

## Prerequisites
```bash
conda install -c bioconda bracken

# Bracken database must match your read length
# Build if not included with Kraken2 database
bracken-build -d /path/to/kraken_db -t 8 -l 150
```

## Quick Start
Tell your AI agent what you want to do:
- "Estimate species abundances from my Kraken2 output"
- "Run Bracken on my classification report with 150bp reads"
- "Compare microbial abundances across my samples"

## Example Prompts
### Basic Abundance Estimation
> "Run Bracken on kraken_report.txt to get species-level abundances"

> "Estimate abundances for my 150bp paired-end reads after Kraken2 classification"

### Multi-sample Analysis
> "Process all my Kraken2 reports through Bracken and combine into one table"

> "Generate relative abundance profiles for all samples in this directory"

### Different Taxonomic Levels
> "Get genus-level abundances instead of species"

> "Run Bracken at family level for my low-coverage samples"

## What the Agent Will Do
1. Verify Kraken2 report and Bracken database are compatible
2. Run Bracken with appropriate read length and taxonomic level
3. Generate abundance table with read counts and fractions
4. Combine multiple sample outputs if requested

## Tips
- Use the closest available read length (100, 150, 200, 250, 300)
- Bracken requires a Kraken2 report file, not the per-read output
- Species level (`-l S`) is default; use `-l G` for genus
- Check `added_reads` column to see how much redistribution occurred
- Normalize to relative abundance for cross-sample comparisons

## Output Columns

| Column | Description |
|--------|-------------|
| name | Taxon name |
| taxonomy_id | NCBI taxon ID |
| taxonomy_lvl | Level (S, G, etc.) |
| kraken_assigned_reads | Direct Kraken2 assignments |
| added_reads | Reads redistributed from higher levels |
| new_est_reads | Total estimated reads |
| fraction_total_reads | Proportion of classified reads |

## Resources
- [Bracken GitHub](https://github.com/jenniferlu717/Bracken)
- [Bracken Paper](https://peerj.com/articles/cs-104/)
