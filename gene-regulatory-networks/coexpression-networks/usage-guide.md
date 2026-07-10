# Co-expression Networks - Usage Guide

## Overview

Build weighted gene co-expression networks to identify modules of co-regulated genes and relate them to phenotypes. WGCNA groups genes into modules by co-expression, summarizes each module by its eigengene, and correlates eigengenes with sample traits to surface biologically relevant modules and hub genes. The insight that bounds interpretation: co-expression measures **marginal** correlation, so a module is a descriptive object (genes that vary together), not a regulatory network -- and the scale-free criterion used to pick the soft power is a heuristic, not proof the biology is scale-free. For direct (not indirect) edges, use a Gaussian graphical model; for directed TF-target regulation, use scenic-regulons or grn-inference.

## Prerequisites

```r
install.packages('WGCNA')
BiocManager::install(c('CEMiTool', 'hdWGCNA'))
install.packages('GeneNet')   # Gaussian graphical model (partial correlation)
```

Minimum 20 samples recommended for reliable module detection (absolute floor of 15). Single-cell data needs metacell aggregation (hdWGCNA), not raw counts.

## Quick Start

Tell your AI agent what you want to do:
- "Build a signed co-expression network from my RNA-seq data"
- "Find gene modules correlated with treatment response"
- "Identify hub genes by module membership in the most significant module"
- "Test whether my modules are preserved in an independent cohort"
- "Find co-expression modules in my single-cell data with hdWGCNA"

## Example Prompts

### Module Detection
> "I have normalized RNA-seq counts from 30 samples. Build a signed WGCNA network with bicor and find gene modules."

> "Run CEMiTool on my expression data for an automated first-pass co-expression analysis."

### Module-Trait Relationships
> "Correlate my WGCNA module eigengenes with survival time and treatment group."

> "Which modules are significantly associated with disease status, and what are their hub genes by kME?"

### Module Preservation
> "I have a discovery and a validation cohort. Test which of my modules are preserved using Zsummary."

### Direct vs Indirect Edges
> "I want only direct gene-gene dependencies, not indirect co-expression. Build a Gaussian graphical model."

### Single-Cell
> "Run hdWGCNA on my Seurat object to find co-expression modules per cell type, capping metacell overlap."

## What the Agent Will Do

1. Load and quality-filter the expression matrix (`goodSamplesGenes`); confirm sample size is adequate
2. Select the soft-thresholding power on the SAME network type used for construction (signed)
3. Construct a signed network with bicor and detect modules, keeping everything in one block
4. Compute module eigengenes and correlate them with sample traits
5. Identify hub genes by module membership (kME), not raw connectivity
6. Test module preservation in independent data (Zsummary / medianRank) where a validation cohort exists

## Tips

- Co-expression is not regulation - a module is "genes that vary together"; an edge can be entirely indirect (a shared upstream driver). Reserve regulatory language for motif/perturbation-supported networks.
- Use signed networks - unsigned merges activators with their anti-correlated repressors into one module. Set `networkType='signed'` at BOTH pickSoftThreshold and blockwiseModules.
- Scale-free R^2 is a heuristic - choose the lowest power where R^2 plateaus above ~0.8; do not interpret it as evidence the network is biologically scale-free.
- Define hubs by kME - signed, bounded, comparable across modules, with a p-value; preferable to raw intramodular connectivity which tracks mean expression.
- Preservation beats detection - any dendrogram yields modules; the scientific claim is reproduction in independent data.
- Batch effects manufacture modules - correct batch before network construction; the top module is often the sequencing run.
- Single-cell needs metacells - naive WGCNA fails on dropout; hdWGCNA aggregates cells but cap `max_shared` to avoid pseudo-replication.

## Related Skills

- differential-networks - compare co-expression structure between conditions (rewiring)
- scenic-regulons - directed TF regulons from single-cell data
- grn-inference - bulk directed GRN inference and TF protein-activity with VIPER
- differential-expression/batch-correction - remove batch effects before network construction
- pathway-analysis/go-enrichment - functional enrichment of gene modules
- single-cell/preprocessing - QC and normalization for hdWGCNA inputs
