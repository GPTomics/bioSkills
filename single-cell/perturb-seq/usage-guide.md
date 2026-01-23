# Perturb-seq Analysis - Usage Guide

## Overview

Analyze CRISPR screening data (Perturb-seq, CROP-seq) integrated with single-cell RNA-seq.

## Prerequisites

```bash
pip install pertpy scanpy
```

```r
install.packages('Seurat')
# Mixscape is included in Seurat v4+
```

## Quick Start

- "Analyze my Perturb-seq CRISPR screen"
- "Find genes affected by each perturbation"
- "Classify perturbed vs unperturbed cells"

## Example Prompts

### Perturbation Analysis

> "Run differential expression for each CRISPR guide"

> "Identify on-target vs off-target effects"

### Cell Classification

> "Use Mixscape to separate perturbed from non-perturbed cells"

> "Calculate perturbation scores for each cell"

### Downstream

> "Find pathways affected by my gene knockouts"

> "Compare perturbation signatures across conditions"

## What the Agent Will Do

1. Load scRNA-seq data with guide assignments
2. Filter low-quality cells and guides
3. Run differential expression per perturbation
4. Calculate perturbation signatures
5. Visualize perturbation effects

## Tips

- **Non-targeting controls** - Essential for background comparison
- **Multiple guides per gene** - Helps distinguish on-target effects
- **MOI < 1** - Low multiplicity ensures one guide per cell
- **Pertpy** - Modern Python framework for perturbation analysis
- **Mixscape** - Seurat's method for classification and scoring
