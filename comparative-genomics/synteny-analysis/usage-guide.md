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

1. Verify assembly quality (N50, BUSCO completeness) before analysis
2. Softmask genomes for repeat content if not already masked
3. Run all-vs-all BLASTP for homology detection
4. Identify collinear gene blocks with MCScanX (or JCVI, i-ADHoRe)
5. Classify syntenic relationships (1:1, 1:many, many:many)
6. Detect structural variants with SyRI if requested
7. Calculate Ks for dating duplications; fit mixture models for WGD peaks
8. Generate visualization plots (dot plots, karyotype views)

## Tips

- **Assembly quality first** - Require N50 > 1 Mb minimum; fragmented assemblies underestimate synteny by up to 40%
- **Repeat masking** - Softmask genomes before BLAST to prevent TE-derived spurious hits
- **Block size** - Minimum 5 genes per block reduces noise; use 10+ for stringent analysis
- **E-value** - Use 1e-10 for close species, 1e-5 for distant comparisons
- **JCVI cscore** - Use 0.99 for reciprocal best hits only; 0.70 (default) for more sensitivity
- **Ks saturation** - Values >2 are unreliable; use mixture models (wgd v2) for formal peak fitting
- **Polyploidy** - Assign subgenomes before comparative analysis; use AnchorWave for WGD-aware alignment
- **Reference-guided scaffolding** - Creates false synteny; avoid comparing genomes scaffolded against each other
- **Tool choice** - MCScanX for general use; i-ADHoRe for ancient WGD; SyRI for structural variants

## Related Skills

- comparative-genomics/positive-selection - dN/dS on syntenic gene pairs
- comparative-genomics/ortholog-inference - Identify orthologs for synteny
- phylogenetics/modern-tree-inference - Phylogenetic context for dating
- genome-annotation/annotation-transfer - Synteny-guided annotation transfer
