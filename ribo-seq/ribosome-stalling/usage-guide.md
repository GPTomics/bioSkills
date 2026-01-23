# Ribosome Stalling - Usage Guide

## Overview

Detect ribosome pausing and stalling sites at codon resolution to study translational regulation, rare codon effects, and nascent chain interactions.

## Prerequisites

```bash
pip install plastid numpy scipy biopython
```

## Quick Start

Tell your AI agent:
- "Find ribosome pause sites in my Ribo-seq data"
- "Calculate codon-specific ribosome occupancy"
- "Identify stalling at rare codons"
- "Analyze pause motifs"

## Example Prompts

### Pause Site Detection

> "Find positions with elevated ribosome occupancy (z-score > 3)"

> "Identify stalling sites in my genes of interest"

> "How many pause sites are there per gene?"

### Codon Analysis

> "Calculate average ribosome occupancy per codon"

> "Which codons have highest pause frequency?"

> "Correlate pausing with tRNA abundance"

### Motif Analysis

> "What amino acid motifs are enriched at pause sites?"

> "Find polyproline-associated stalling"

> "Analyze context around pause sites"

## What the Agent Will Do

1. Map reads to P-site positions using offset
2. Calculate codon-level ribosome occupancy
3. Identify positions with elevated occupancy (z-score)
4. Aggregate by codon type
5. Analyze sequence/motif context

## Tips

- **Correct P-site offset** is critical for codon assignment
- **Z-score > 3** is a typical threshold for pause sites
- **Normalize by gene expression** to compare across genes
- **Polyproline** (PPP) is a well-known pause motif
- **Rare codons** correlate with pausing (tRNA limitation)
