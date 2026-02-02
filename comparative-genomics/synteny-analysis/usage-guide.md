# Synteny Analysis - Usage Guide

## Overview

Analyze genome collinearity and conserved gene order between species using MCScanX, SyRI, and JCVI for evolutionary and comparative genomics studies.

## Prerequisites

```bash
# MCScanX (compile from source)
git clone https://github.com/wyp1125/MCScanX
cd MCScanX && make
export PATH=$PATH:$(pwd)

# SyRI for structural variants
pip install syri

# JCVI for visualization
pip install jcvi

# BLAST for homology search
conda install -c bioconda blast
```

## Quick Start

Tell your AI agent what you want to do:
- "Find syntenic blocks between human and mouse genomes"
- "Detect whole-genome duplications in my plant genome"
- "Identify chromosomal rearrangements between two assemblies"

## Example Prompts

### Basic Synteny

> "Find collinear gene blocks between Arabidopsis and rice"

> "Create a synteny dot plot comparing my two genomes"

### WGD Detection

> "Look for whole-genome duplication signatures in this genome"

> "Calculate Ks distribution for syntenic gene pairs"

### Structural Variants

> "Identify inversions and translocations between reference and query"

> "Find structural rearrangements between chromosome assemblies"

## What the Agent Will Do

1. Prepare genome and annotation files for analysis
2. Run all-vs-all BLASTP for homology detection
3. Identify collinear gene blocks with MCScanX
4. Classify syntenic relationships (1:1, 1:many)
5. Detect structural variants with SyRI if requested
6. Calculate Ks for dating duplications
7. Generate visualization plots

## Tips

- **Block size** - Minimum 5 genes per block reduces noise; use 10+ for stringent analysis
- **E-value** - Use 1e-10 for close species, 1e-5 for distant comparisons
- **Ks saturation** - Values >2 are unreliable; filter before interpretation
- **WGD peaks** - Multiple Ks peaks suggest multiple WGD events
- **Visualization** - JCVI produces publication-quality synteny plots
