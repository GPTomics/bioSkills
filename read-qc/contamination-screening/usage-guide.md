# Contamination Screening - Usage Guide

## Overview
FastQ Screen quickly identifies contamination by aligning a subset of reads against multiple reference genomes. It detects cross-species contamination, bacterial contamination, adapter sequences, and sample swaps.

## Prerequisites
```bash
conda install -c bioconda fastq-screen

# Download pre-built databases
fastq_screen --get_genomes
```

## Quick Start
Tell your AI agent what you want to do:
- "Screen my FASTQ files for contamination"
- "Check if my human sample has any mouse contamination"
- "Verify my samples are not swapped"

## Example Prompts

### Basic Screening
> "Run FastQ Screen on my samples to check for contamination"

> "Screen my FASTQ files against human, mouse, and common contaminants"

### Troubleshooting
> "My RNA-seq has unexpected results, check for contamination or sample swap"

> "Identify the source of contamination in my samples"

### Configuration
> "Set up FastQ Screen with custom genome databases"

> "Add E. coli and yeast to my contamination screening panel"

## What the Agent Will Do
1. Configure FastQ Screen with appropriate reference databases
2. Run screening on a subset of reads from each sample
3. Generate reports showing percentage mapping to each genome
4. Identify problematic samples with contamination or swaps
5. Recommend appropriate remediation steps

## Interpretation Guide

### Good Results (Human Sample)
```
Genome      %One_hit  %Multi_hit
Human       95.0      2.0
Mouse       0.1       0.0
Ecoli       0.0       0.0
Adapters    0.1       0.0
```

### Common Issues

| Issue | Pattern | Solution |
|-------|---------|----------|
| Bacterial contamination | High E.coli/bacteria % | Investigate source |
| Sample swap | Wrong species dominant | Re-sequence or exclude |
| Adapter contamination | High Adapters % | Run adapter trimming |
| PhiX spike-in | High PhiX % | Filter PhiX reads |
| rRNA contamination | High rRNA % | rRNA depletion failed |

## Tips
- Run contamination screening early before investing in alignment and analysis
- Include common contaminants: E. coli, yeast, PhiX, adapters
- For xenograft samples, expect mixed human/mouse signal
- Sample swaps are easier to detect than fix - verify sample identity early
- Use `--subset` to screen a random subset of reads for faster results

## Resources
- [FastQ Screen Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
- [Pre-built Databases](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.15.3.tar.gz)
