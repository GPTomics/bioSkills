# Genome Assembly Pipeline - Usage Guide

## Overview

This workflow assembles genomes from sequencing reads, supporting short reads, long reads, and hybrid approaches with polishing and quality assessment.

## Prerequisites

```bash
conda install -c bioconda spades flye minimap2 bwa samtools pilon medaka quast busco
```

## Quick Start

Tell your AI agent what you want to do:
- "Assemble my bacterial genome from Illumina reads"
- "Run de novo assembly on my Nanopore data"
- "Polish my assembly with short reads"

## Example Prompts

### Assembly
> "Assemble my bacterial genome with SPAdes"

> "Use Flye to assemble my ONT reads"

> "Create a hybrid assembly from long and short reads"

### QC
> "Run QUAST on my assembly"

> "Check completeness with BUSCO"

> "What's the N50 of my assembly?"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Short reads | FASTQ | Illumina paired-end |
| Long reads | FASTQ | ONT or PacBio |
| Expected size | Number | Approximate genome size |

## What the Workflow Does

1. **Quality Control** - Assess and filter reads
2. **Assembly** - De novo contig/scaffold generation
3. **Polishing** - Error correction
4. **Assessment** - Contiguity and completeness metrics

## Assembly Strategy Selection

| Data Type | Assembler | Polisher |
|-----------|-----------|----------|
| Illumina only | SPAdes | Pilon (optional) |
| ONT raw | Flye | medaka + Pilon |
| ONT HQ | Flye | medaka |
| PacBio HiFi | Flye or Hifiasm | Usually not needed |
| Hybrid | Flye + short polish | Pilon |

## Tips

- **Coverage**: 50-100x for bacteria, 30-50x for larger genomes
- **Long reads**: Improves contiguity dramatically
- **BUSCO lineage**: Choose appropriate lineage for your organism
- **Contamination**: Check for unexpected sequences
