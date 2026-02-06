# Structure Probing - Usage Guide

## Overview
Process experimental RNA structure probing data from SHAPE-MaP and DMS-MaPseq experiments. ShapeMapper2 converts raw sequencing reads into per-nucleotide reactivity profiles that reflect single-stranded vs base-paired status, which can then constrain computational structure prediction.

## Prerequisites
```bash
# ShapeMapper2 (Linux only; use Docker on macOS)
conda install -c bioconda shapemapper2

# Docker alternative for macOS
docker pull shapemapper2/shapemapper2

# DMS-MaPseq analysis
pip install seismic-rna

# ViennaRNA for SHAPE-constrained folding
conda install -c bioconda viennarna

# Visualization
pip install matplotlib pandas numpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Process my SHAPE-MaP sequencing data to get reactivity profiles"
- "Run ShapeMapper2 on my modified and untreated FASTQ files"
- "Use my SHAPE reactivities to constrain RNA structure prediction"
- "Analyze my DMS-MaPseq data with SEISMIC-RNA"
- "Plot the reactivity profile for my RNA"

## Example Prompts

### SHAPE-MaP Analysis
> "I have SHAPE-MaP paired-end reads for modified and untreated samples targeting my RNA. Run ShapeMapper2 to get reactivities."

> "Process my amplicon SHAPE-MaP data and plot the reactivity profile."

### Structure-Constrained Folding
> "Use my SHAPE reactivity data to constrain RNAfold and predict a more accurate secondary structure."

> "Compare the unconstrained MFE structure with the SHAPE-constrained prediction for my RNA."

### DMS-MaPseq
> "Process my DMS-MaPseq data using SEISMIC-RNA to get mutation rates at A and C positions."

> "Cluster my DMS-MaPseq data to detect multiple RNA conformations."

### Visualization
> "Plot the SHAPE reactivity profile with color coding for paired and unpaired regions."

> "Generate an arc diagram of my RNA structure colored by SHAPE reactivity."

## What the Agent Will Do
1. Run ShapeMapper2 (or SEISMIC-RNA for DMS) with modified and control samples
2. Assess data quality (read depth, mutation rates, modification efficiency)
3. Generate per-nucleotide reactivity profiles
4. Use reactivities to constrain ViennaRNA folding
5. Visualize reactivity profiles and structure predictions

## Tips
- **Read depth** - Minimum 5,000 reads per nucleotide recommended for reliable mutation rate estimation
- **Controls** - An untreated control is required; a denatured control improves normalization
- **macOS** - ShapeMapper2 is Linux-only; use Docker or Singularity on macOS
- **DMS coverage** - DMS only probes A and C residues (~50% of positions), so structure constraints are sparser than SHAPE
- **Multiple conformations** - SEISMIC-RNA can cluster mutations to detect alternative RNA structures; consider this if reactivity profiles seem inconsistent with a single structure
- **Reagent parameters** - Use m=1.8, b=-0.6 for standard SHAPE reagents (1M7, NAI); adjust for DMS (m=1.1, b=-0.3)

## Related Skills
- secondary-structure-prediction - Unconstrained and SHAPE-constrained RNA folding
- clip-seq/binding-site-annotation - Protein-RNA binding site annotation
- epitranscriptomics/m6a-peak-calling - RNA modification detection
