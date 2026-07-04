# Temporal Gene Clustering - Usage Guide

## Overview

Groups PRE-SELECTED temporally variable genes by expression-profile SHAPE into candidate co-expression programs. Clustering is descriptive and unsupervised - it has no null and always returns clusters - so it is strictly downstream of gene selection (differential-expression/timeseries-de or a variance filter), NOT a test of which genes are dynamic. Supports soft (fuzzy) clustering with Mfuzz, hard/soft clustering with TCseq, automatic-k hierarchical clustering with DEGreport, and DTW/soft-DTW clustering with tslearn for phase-shifted patterns.

## Prerequisites

### R
```r
BiocManager::install(c('Mfuzz', 'TCseq', 'DEGreport'))
```

### Python
```bash
pip install tslearn scikit-learn numpy matplotlib
```

### Data Requirements
- Expression matrix (genes x timepoints), typically averaged across replicates
- Pre-filtered for temporally variable genes (e.g., via limma or DESeq2 LRT)
- At least 4 timepoints; 6-12 timepoints is ideal for meaningful clustering

## Quick Start

Tell the AI agent what to cluster:
- "Cluster my time-course genes by expression profile shape using Mfuzz"
- "Group my differentially expressed genes into temporal response patterns"
- "Find genes with similar trajectories using DTW clustering"
- "Run degPatterns on my RNA-seq time-course data"

## Example Prompts

### Soft Clustering
> "I have 2000 temporally variable genes across 8 timepoints. Cluster them with Mfuzz soft clustering and filter for core members."

> "Run fuzzy c-means clustering on my time-series expression data and show me the cluster centroids."

### Automatic Clustering
> "Use DEGreport degPatterns to automatically determine the number of temporal clusters in my RNA-seq data."

> "Cluster my time-course genes and automatically pick the best number of clusters."

### DTW-Based Clustering
> "Some of my genes have phase-shifted responses. Cluster them using DTW distance in Python."

> "Run time-series k-means with dynamic time warping on my expression profiles."

### Post-Clustering Analysis
> "After clustering, run GO enrichment on each temporal cluster to identify pathway themes."

> "Show me which transcription factors are enriched in each temporal expression cluster."

## What the Agent Will Do

1. Confirm the input is pre-selected temporally variable genes; if it is the full matrix, prefilter (DE hits or top-variance) before clustering
2. Standardize expression profiles (z-score per gene across timepoints)
3. Choose a distance metric (Euclidean-on-zscore / correlation / constrained DTW) appropriate to whether phase shifts are expected
4. Estimate and VALIDATE the fuzzifier (Mfuzz) and triangulate k across indices, biology, and bootstrap stability (not a single index)
5. Run clustering (Mfuzz, TCseq, DEGreport, or tslearn) and filter genes by membership
6. Interpret centroids as candidate programs; generate profile plots
7. Run per-cluster enrichment using the INPUT gene set as background, avoiding the double-dipping trap
8. Export cluster assignments for downstream analysis

## Tips

- Prefilter FIRST: cluster only temporally variable genes (timeseries-DE hits or a variance filter). Clustering all genes is the #1 error - z-scoring amplifies flat-gene noise into fake programs, and clustering has no null so it always returns clusters.
- Always standardize (z-score per gene) before clustering; raw values let high-expression genes dominate on magnitude instead of shape.
- Start with Mfuzz (soft) for most tasks; its continuous membership exposes ambiguous genes rather than hiding them in a hard label.
- Do not hardcode m=2 and do not trust mestimate() blindly: it implements Schwaemmle & Jensen (2010), is dominated by the number of timepoints, and can go degenerate at extreme D. Inspect the returned m and the fraction of genes clearing the acore cutoff.
- Membership threshold 0.5 for Mfuzz is a convention; relax to 0.3 for exploratory work (admits more noise). Always report the retained fraction.
- The distance metric matters more than the algorithm. Default to Euclidean-on-zscore or correlation. Use DTW only when phase shifts are biologically expected (e.g., signaling cascades), always with a Sakoe-Chiba band - unconstrained DTW invents structure. Soft-DTW is not faster (still quadratic); its value is a differentiable, well-defined average.
- Pick k by triangulation, not one index: min centroid distance, silhouette (scored under the clustering metric), biology, and above all bootstrap/consensus STABILITY.
- Do not test clusters for the temporal signal used to select them (double-dipping inflates p-values); test only against independent annotations.
- Run per-cluster GO/GSEA with the INPUT gene set as background, never the whole genome - genome-as-background rediscovers the selection.
- For large gene sets (>10,000), pre-filter to the top variable genes to keep clustering tractable (and correct).

## Related Skills

- temporal-genomics/circadian-rhythms - Rhythm-specific clustering by phase
- temporal-genomics/trajectory-modeling - Continuous trajectory fitting
- differential-expression/timeseries-de - Upstream temporal DE for gene selection
- pathway-analysis/go-enrichment - Per-cluster functional enrichment
