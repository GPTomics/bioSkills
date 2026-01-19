# Motif Analysis - Usage Guide

## Overview

Motif analysis identifies DNA sequence patterns enriched in ChIP-seq or ATAC-seq peaks. This reveals transcription factor binding sites and regulatory elements.

## When to Use This Skill

- After ChIP-seq peak calling to validate target TF
- Discover co-binding TFs in ChIP-seq data
- Analyze ATAC-seq peaks for accessible TF sites
- Characterize regulatory regions
- Compare motif enrichment between conditions

## Installation

```bash
# HOMER
conda install -c bioconda homer
perl /path/to/homer/configureHomer.pl -install hg38

# MEME Suite
conda install -c bioconda meme
```

## Quick Start

### HOMER (Recommended)

```bash
findMotifsGenome.pl peaks.bed hg38 motif_output/ -size 200
```

### MEME-ChIP

```bash
bedtools getfasta -fi genome.fa -bed peaks.bed -fo peaks.fa
meme-chip -oc meme_output peaks.fa
```

## Key Concepts

### De Novo vs Known Motif Analysis

| Type | Purpose | Tools |
|------|---------|-------|
| De novo | Discover new motifs | MEME, HOMER (homer*) |
| Known | Test enrichment of database motifs | HOMER (known*), CentriMo |

### Motif Representation

- **Consensus**: Simple sequence (e.g., CACGTG)
- **PWM**: Position Weight Matrix (probabilities)
- **Logo**: Visual representation of PWM

## Interpretation

### Good Results
- P-value < 1e-10
- Target % significantly > Background %
- Expected TF motif ranks highly
- Biologically relevant co-factors

### Red Flags
- No significant motifs (check peak quality)
- Only repeat-like motifs (mask repeats)
- Background similar to target (wrong control)

## Resources

- [HOMER Documentation](http://homer.ucsd.edu/homer/motif/)
- [MEME Suite](https://meme-suite.org/)
- [JASPAR Database](https://jaspar.genereg.net/)
