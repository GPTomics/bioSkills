# Contamination Screening - Usage Guide

## Overview

FastQ Screen quickly identifies contamination by aligning a subset of reads against multiple reference genomes. It detects cross-species contamination, bacterial contamination, adapter sequences, and sample swaps.

## When to Use This Skill

- New sequencing data from external facility
- Unexpected results in downstream analysis
- Quality control before large-scale processing
- Identifying sample swaps or mix-ups
- Checking for bacterial contamination

## Installation

```bash
conda install -c bioconda fastq-screen

# Download pre-built databases
fastq_screen --get_genomes
```

## Basic Workflow

### 1. Configure Databases

Edit `~/.fastq_screen.conf` or create a project config:

```
DATABASE  Human    /path/to/genomes/Human/GRCh38
DATABASE  Mouse    /path/to/genomes/Mouse/GRCm39
DATABASE  Ecoli    /path/to/genomes/Ecoli/K12
DATABASE  PhiX     /path/to/genomes/PhiX/phix
DATABASE  Adapters /path/to/genomes/Adapters/adapters
BOWTIE2   /path/to/bowtie2
THREADS   8
```

### 2. Run Screening

```bash
fastq_screen --outdir results/ sample.fastq.gz
```

### 3. Interpret Results

Open `sample_screen.html` or view `sample_screen.txt`.

## Interpretation Guide

### Good Results (Human Sample)

```
Genome      %One_hit  %Multi_hit
Human       95.0      2.0
Mouse       0.1       0.0
Ecoli       0.0       0.0
Adapters    0.1       0.0
```

### Problematic Results

**Bacterial Contamination:**
```
Genome      %One_hit
Human       70.0
Ecoli       25.0      <- Problem!
```

**Sample Swap:**
```
Genome      %One_hit
Human       5.0       <- Expected Human
Mouse       90.0      <- Got Mouse!
```

**Adapter Contamination:**
```
Genome      %One_hit
Human       80.0
Adapters    15.0      <- Trim adapters
```

## Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| High adapter % | Run Cutadapt/fastp |
| High PhiX % | Filter PhiX reads |
| High E.coli % | Investigate source |
| Multiple species | Check for contamination/swap |
| High rRNA % | rRNA depletion may have failed |

## Resources

- [FastQ Screen Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
- [Pre-built Databases](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.15.3.tar.gz)
