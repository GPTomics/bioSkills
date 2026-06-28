# Spatial Statistics - Usage Guide

## Overview

This skill quantifies spatial structure in spatial transcriptomics data: spatially variable genes (SVG), spatial autocorrelation (Moran's I, Geary's C, Getis-Ord, local LISA), and cell-type colocalization (neighborhood enrichment, co-occurrence). It is decision-first, not recipe-first. Two traps drive every choice: a spatially variable gene is usually just a marker of a spatially clustered cell type rather than a spatially regulated gene, and the neighbor graph plus the permutation null -- both researcher choices -- silently determine the result. The skill helps choose an SVG method by its null and scaling, choose the right autocorrelation statistic for the question (clustering vs hot/cold vs local), and choose a colocalization null strong enough to defeat the abundance and tissue-compartment confounds.

## Prerequisites

```bash
pip install squidpy scanpy anndata esda libpysal
# SVG methods that scale or test alternative nulls are R/Bioconductor:
# R: BiocManager::install(c('SPARK', 'nnSVG'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Find spatially variable genes and tell me which are just cell-type markers"
- "Compute Moran's I for my spatial data"
- "Map hot spots of this gene, not just whether it clusters"
- "Test whether these two cell types specifically colocalize, controlling for abundance"
- "Run local LISA and handle the multiple-testing problem"

## Example Prompts

### Spatially Variable Genes
> "Rank genes by spatial autocorrelation, but separate genes that are spatially variable because of cell-type composition from genes regulated within a cell type"

> "Pick an SVG method for my Xenium data with 300k cells -- I need something that scales"

### Choosing the Statistic
> "I want to know where this gene is HIGH, not just whether it clusters -- which statistic do I use"

> "My global Moran's I is near zero but the tissue is obviously structured -- what am I missing"

### Cell-Type Colocalization
> "Test whether tumor-associated macrophages specifically associate with tumor cells, controlling for the fact that both are abundant in the tumor bed"

> "My neighborhood-enrichment z-score is positive -- is that a real interaction or a compartment artifact"

## What the Agent Will Do

1. Establish the neighbor graph and confirm coordinate units (microns vs pixels), since the graph defines every statistic.
2. For SVG, choose a method by null and scaling, then warn that sample-wide SVGs are largely cell-type markers and compare to the HVG list.
3. For autocorrelation, select Moran/Geary (clustering), Getis-Ord (hot/cold), or local LISA (non-stationary tissue) to match the question.
4. For colocalization, run neighborhood enrichment or co-occurrence and interrogate the null, demanding survival under a stronger null before claiming specific affinity.
5. Apply FDR correctly, with stricter cutoffs for local statistics, and threshold on effect size, not p alone.

## Tips

- A spatially variable gene is not necessarily spatially regulated; sample-wide SVG lists mostly re-derive cell-type markers and overlap heavily with HVGs. The payoff is the SVG-not-HVG subset (modest-amplitude gradients).
- Moran's I and Geary's C cannot tell a hot spot from a cold spot; use Getis-Ord Gi* when the sign matters.
- A non-significant global Moran's I does not mean "no spatial structure" -- opposite-sign regions cancel; use local LISA/Gi* on non-stationary tissue.
- Local statistics carry a double trap: n tests AND the statistics are spatially autocorrelated, so naive BH-FDR is wrong. Use stricter cutoffs (0.001) and treat the map as exploratory.
- The default Squidpy `nhood_enrichment` null only beats complete randomness; abundant co-compartment cell types pass trivially. Demand a conditional/within-compartment or toroidal null for a specific-affinity claim.
- Be most skeptical of enrichment involving rare cell types: few edges produce large, unstable z-scores.
- SVG methods disagree because they test different nulls; that is expected. Report graph-robust hits and expect cross-method intersection to be small.
- Confirm coordinates are in microns before setting any radius or length-scale parameter; Visium array coordinates are not microns.

## Related Skills

- spatial-neighbors - builds the graph W that every statistic here inherits
- spatial-domains - region-level structure, distinct from a colocalization result
- spatial-communication - ligand-receptor in space, where the same confounds recur
- single-cell/markers-annotation - cell-type labels and the marker overlap that confounds SVG
