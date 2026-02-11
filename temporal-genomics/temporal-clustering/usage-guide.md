# Temporal Gene Clustering - Usage Guide

## Overview

Groups genes with similar temporal expression profiles into clusters, revealing coordinated regulatory programs across time-course experiments. Supports soft (fuzzy) clustering with Mfuzz, automatic cluster selection with DEGreport, and DTW-based clustering with tslearn for phase-shifted patterns.

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

1. Load the expression matrix and verify timepoint structure
2. Filter for temporally variable genes if not already filtered
3. Standardize expression profiles (z-score across timepoints)
4. Estimate optimal cluster number (silhouette, gap statistic, or degPatterns auto-selection)
5. Run clustering (Mfuzz, TCseq, DEGreport, or tslearn)
6. Filter genes by membership score or cluster quality
7. Generate cluster profile plots and membership heatmaps
8. Export cluster assignments for downstream enrichment analysis

## Tips

- Always standardize before clustering; raw expression values cause high-expression genes to dominate cluster assignments
- Start with Mfuzz for most temporal clustering tasks; its soft membership handles ambiguous genes gracefully
- The fuzzifier parameter m should be estimated from data with mestimate(); do not hardcode a value
- Use DTW clustering only when phase shifts are biologically expected (e.g., signaling cascades); it is much slower than Euclidean
- Membership threshold of 0.5 for Mfuzz is standard but can be relaxed to 0.3 for exploratory analyses
- Run per-cluster GO enrichment to validate that clusters capture distinct biological processes
- For large gene sets (>10,000), pre-filter to the top 2000-5000 variable genes to keep clustering tractable

## Related Skills

- temporal-genomics/circadian-rhythms - Rhythm-specific clustering by phase
- temporal-genomics/trajectory-modeling - Continuous trajectory fitting
- differential-expression/timeseries-de - Upstream temporal DE for gene selection
- pathway-analysis/go-enrichment - Per-cluster functional enrichment
