# Spatial Analysis - Usage Guide

## Overview

Tests cell-cell co-location, discovers tissue niches, and quantifies spatial structure from segmented IMC/MIBI cell tables. The point that shapes every result here: a spatial interaction or niche is a hypothesis test against a null model, and the null chosen (label-shuffle vs CSR vs density-conditioned vs patient-level) silently decides whether the result is real biology, a density gradient, a segmentation artifact, or upstream clustering -- so the skill forces two questions onto every analysis: against which null, and at which unit.

## Prerequisites

```bash
pip install squidpy scanpy anndata scikit-learn numpy
# R for inhomogeneous Ripley's K and imcRtools graphs (optional)
# install.packages('spatstat'); BiocManager::install('imcRtools')
```

## Quick Start

Tell your AI agent what you want to do:
- "Test whether CD8 T cells co-locate with tumor cells, with the right null"
- "Build a contact graph and prune the long edges"
- "Discover cellular neighborhoods and check they are stable to window size"
- "Tell me if this interaction is real or just a density artifact"
- "Compute per-image neighborhood enrichment to aggregate across patients"

## Example Prompts

### Neighborhood enrichment
> "Run neighborhood enrichment between all cell types per image. Tell me which null this uses and what it cannot control for."

### Density confound
> "My CD8 and Treg cells both live in the tumor margin. Is their 'interaction' real, or are they just sharing a compartment? Cross-check with inhomogeneous Ripley's K."

### Niche discovery
> "Find cellular neighborhoods in my cohort and show me whether the niches survive changing the window from 10 to 30 neighbors."

### Graph choice
> "I want to test contact-dependent (juxtacrine) interactions. Which graph should I build and how do I avoid edges across empty tissue?"

## What the Agent Will Do

1. Build a spatial graph matching the biological claim -- pruned Delaunay for contact, fixed radius for paracrine proximity -- per image.
2. Run neighborhood enrichment per image, naming the null and treating the z-score as a per-image summary (never thresholded across images of different size).
3. Cross-check density-driven hits with inhomogeneous cross-K.
4. Discover cellular neighborhoods, sweeping the window size and reporting stability.
5. Hand cross-condition comparison (does a niche/interaction differ between groups) to differential-analysis at the patient unit, with FDR across pairs.

## Tips

- The label-shuffle null controls abundance but not local density -- two types sharing a compartment will look "interacting".
- The squidpy z-score is unbounded and graph-degree-dependent; it is not comparable across images of different cell counts.
- Delaunay invents long edges across lumen/necrosis; prune by a biological max distance before counting.
- kNN graphs mix contact and proximity by density; state the micron scale of any "interaction".
- A cellular neighborhood has no built-in significance test; it is exploratory until compared across patients.
- Verify that niche-defining markers are not one cell's signal bleeding into a neighbor (spillover).

## Related Skills

- phenotyping - cell-type labels are the input to every spatial test
- differential-analysis - testing whether interactions/niches differ between conditions at the patient level
- cell-segmentation - over-segmentation and lateral spillover create fake niches
- data-preprocessing - uncompensated spillover manufactures false cell-cell interactions
- spatial-transcriptomics/spatial-statistics - shared squidpy neighborhood and autocorrelation methods
