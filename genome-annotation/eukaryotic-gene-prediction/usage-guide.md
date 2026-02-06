# Eukaryotic Gene Prediction - Usage Guide

## Overview

Predict protein-coding genes in eukaryotic genomes using BRAKER3 (RNA-seq + protein evidence), GALBA (protein-only), or Augustus (pre-trained models). Produces gene models in GFF3/GTF format with protein and CDS sequences.

## Prerequisites

```bash
# BRAKER3 (includes Augustus and GeneMark-ETP)
conda install -c bioconda braker3

# GALBA (protein-only alternative)
conda install -c bioconda galba

# Augustus standalone
conda install -c bioconda augustus

# Evidence preparation
conda install -c bioconda hisat2 samtools

# Evaluation
conda install -c bioconda busco

# Python parsing
pip install gffutils pandas
```

Requires a softmasked genome assembly. Run repeat-annotation/RepeatMasker first.

## Quick Start

Tell your AI agent what you want to do:
- "Predict genes in my eukaryotic genome using RNA-seq and protein evidence"
- "Run BRAKER3 on my softmasked plant genome"
- "Predict genes using only protein homology with GALBA"

## Example Prompts

### BRAKER3 with Full Evidence

> "Run BRAKER3 on my softmasked genome with these RNA-seq BAM files and OrthoDB proteins"

> "Predict genes in my fungal genome using BRAKER3 with the --fungus flag"

### Protein-Only

> "I only have protein evidence - run GALBA on my assembly with closely related species proteins"

> "Predict genes using OrthoDB vertebrate proteins without RNA-seq"

### Evaluation

> "Run BUSCO on my predicted proteins to check annotation completeness"

> "Compare gene model statistics between BRAKER3 and GALBA predictions"

## What the Agent Will Do

1. Verify the assembly is softmasked (repeats in lowercase)
2. Prepare evidence (align RNA-seq, download OrthoDB proteins if needed)
3. Run BRAKER3 or GALBA with appropriate parameters
4. Extract protein and CDS sequences from predictions
5. Evaluate with BUSCO and report gene model statistics
6. Flag issues (low gene count, many single-exon genes, poor BUSCO scores)

## Tips

- **Softmasking is essential** - Always run RepeatMasker before gene prediction to avoid false positives in repeats
- **RNA-seq quality matters** - Higher RNA-seq depth and tissue diversity improve predictions significantly
- **OrthoDB clade** - Choose the most specific clade available (Vertebrata > Metazoa for fish)
- **FASTA headers** - Remove special characters from contig names before running BRAKER3
- **Fungal genomes** - Always use `--fungus` flag which adjusts intron size expectations
- **UTR prediction** - Only enable `--UTR` with high RNA-seq coverage (>50x)
- **Multiple BAMs** - Provide RNA-seq from multiple tissues as comma-separated BAM list

## Related Skills

- genome-annotation/repeat-annotation - Softmask repeats before gene prediction
- genome-annotation/functional-annotation - Add functional annotation to predicted genes
- read-alignment/star-alignment - Alternative RNA-seq aligner
- genome-assembly/assembly-qc - Verify assembly quality first
