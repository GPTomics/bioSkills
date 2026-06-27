# Multimodal Integration - Usage Guide

## Overview

Multimodal single-cell assays measure several biological layers per cell (RNA + surface protein in CITE-seq, RNA + chromatin accessibility in 10x Multiome) or pair independent single-modality datasets (unpaired scRNA + scATAC). This skill classifies the integration task by anchor structure, denoises each modality natively, and selects a joint method (WNN, totalVI, MultiVI, MOFA+, GLUE, Seurat v5 bridge) with explicit failure modes.

## Prerequisites

```r
install.packages(c('Seurat', 'dsb'))
BiocManager::install('Signac')          # Multiome ATAC
```

```bash
pip install muon mudata scanpy anndata scvi-tools
pip install scglue                       # unpaired/diagonal integration
```

## Quick Start

Tell your AI agent what you want to do:
- "Classify my multimodal task as paired or unpaired and pick a method"
- "Denoise my CITE-seq ADT with DSB before clustering"
- "Run WNN on my CITE-seq data and show me the modality weights"
- "Integrate independent scRNA and scATAC that have no shared cells"
- "Decide between totalVI and WNN for my CITE-seq experiment"

## Example Prompts

### CITE-seq (RNA + Protein)
> "My ADT clusters look like background smears; denoise with DSB using empty droplets, then run WNN"
> "Train totalVI on my CITE-seq MuData and give me denoised protein plus foreground probabilities"
> "WNN clustering seems driven by CD markers only; show the per-cell modality weight distribution and test stability if ADT is down-weighted"

### Multiome (RNA + ATAC, same cell)
> "Process RNA with PCA and ATAC with TF-IDF/LSI, drop depth-correlated LSI components, then join with WNN"
> "I merged two multiome runs and see batch structure; re-quantify against a unified peak set"

### Unpaired / Diagonal / Mosaic
> "Align my independent scRNA and scATAC with GLUE using a peak-near-gene guidance graph"
> "Map my scATAC query onto an scRNA reference using a multiome bridge dataset"
> "Integrate a mosaic design where one batch has RNA+ATAC and another has RNA only"

## What the Agent Will Do

1. Classify the task by anchor structure (vertical/paired, diagonal/unpaired, mosaic) to narrow the method class.
2. Run per-modality QC and intersect cells across modalities (paired) or keep them separate (diagonal).
3. Denoise ADT with DSB (or totalVI's built-in background mixture) before any joint embedding.
4. Reduce each modality in its native pipeline (PCA for RNA, TF-IDF/LSI for ATAC).
5. Select and run the joint method that matches the anchor structure, verifying defaults against installed docs.
6. Report per-cell modality weights and test clustering stability against down-weighting a suspect modality.
7. Flag any imputed modality (gene activity, MultiVI imputation) as inference, not measurement.

## Decision Guidance

### Anchor structure first
- Paired (same cells): WNN, totalVI, MultiVI, MOFA+, mojitoo.
- Unpaired/diagonal (no shared cells or features): GLUE, Seurat v5 bridge.
- Mosaic (some modalities only): MultiVI, StabMap, Cobolt.

### Paired CITE-seq method choice
- totalVI: need denoised protein, principled DE, batch integration, or to merge different antibody panels (GPU helps).
- WNN: fast assumption-light joint clustering of one well-normalized dataset; pair with DSB because WNN does not denoise protein.
- MOFA+: interpret shared vs modality-specific axes of variation, not for clustering or denoising.

### ADT normalization
- DSB when raw/unfiltered matrices (empty droplets) are available; it removes background.
- CLR (margin=2) as a quick fallback; it rescales but does not remove background.

## Tips

- **Classify by anchor first** - the anchor structure (shared cells vs features vs neither) decides the method class before any code.
- **ADT background is three parts** - ambient (in empties), cell-intrinsic non-specific binding (not in empties), and spillover; one per-cell factor cannot separate all three.
- **Denoise protein before embedding** - raw or CLR-only ADT carries background into the joint graph; use DSB or totalVI.
- **Watch for modality domination** - WNN rewards local predictability, which saturating ADT features can fake; always inspect per-cell weights.
- **Unified peaks for merged Multiome** - re-quantify all cells against one peak set or peak-boundary mismatch creates fake batches.
- **Gene activity is not RNA** - the ATAC gene-activity matrix is an approximation; do not treat it as a measured transcriptome.
- **Imputed is not measured** - MultiVI/StabMap imputation and gene-activity scores are model inferences; flag DE on them as model-dependent.
- **RNA-protein discordance can be biology** - protein lags transcription and has its own half-life; single-gene mismatch is often informative, not an artifact.
- **GLUE needs matched genome coordinates** - a build mismatch silently empties the guidance graph and ruins alignment.

## Related Skills

single-cell/scatac-analysis - ATAC QC, TF-IDF/LSI, gene-activity caveats for the Multiome ATAC half
single-cell/preprocessing - per-modality RNA QC and normalization before integration
single-cell/clustering - clustering and UMAP on the joint graph
single-cell/batch-integration - horizontal (same-modality, cross-sample) correction
single-cell/markers-annotation - marker-based interpretation of joint clusters
atac-seq/motif-deviation - chromVAR TF activity on the Multiome ATAC modality
pathway-analysis/go-enrichment - functional interpretation of modality-specific factors
