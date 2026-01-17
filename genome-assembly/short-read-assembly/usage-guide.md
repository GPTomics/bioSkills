# Short-Read Assembly - Usage Guide

## Overview

SPAdes assembles genomes from Illumina short reads using a multi-k-mer de Bruijn graph approach. It's the standard tool for bacterial, fungal, and small genome assembly.

## When to Use This Skill

- Bacterial genome assembly from Illumina data
- Fungal genome assembly
- Plasmid extraction
- Metagenome assembly
- Transcriptome assembly (RNA-seq)
- Hybrid assembly with short + long reads

## Installation

```bash
conda install -c bioconda spades
```

## Quick Start

```bash
# Basic assembly
spades.py -1 R1.fastq.gz -2 R2.fastq.gz -o assembly_output

# Bacterial isolate (recommended)
spades.py --isolate --careful -1 R1.fq.gz -2 R2.fq.gz -o output
```

## Mode Selection

| Mode | Flag | Use Case |
|------|------|----------|
| Default | (none) | General purpose |
| Isolate | `--isolate` | Single organism, uniform coverage |
| Careful | `--careful` | Reduce misassemblies |
| Meta | `--meta` | Metagenomes |
| RNA | `--rna` | Transcriptomes |
| Plasmid | `--plasmid` | Extract plasmids |

## Typical Workflow

1. QC reads with FastQC
2. Trim adapters with fastp/Trimmomatic
3. Assemble with SPAdes
4. Assess quality with QUAST/BUSCO
5. Polish if needed (Pilon)

## Output

Primary output: `scaffolds.fasta`

Header format: `>NODE_1_length_500000_cov_50.5`

## Resources

- [SPAdes Manual](https://github.com/ablab/spades)
