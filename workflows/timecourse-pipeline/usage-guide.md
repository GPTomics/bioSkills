# Time-Course Pipeline - Usage Guide

## Overview
Complete time-course expression analysis workflow from expression matrix to temporal patterns and pathway enrichment. Covers temporal differential expression (limma splines or DESeq2 LRT), Mfuzz soft clustering of expression profiles, optional circadian rhythm detection with MetaCycle or CosinorPy, GAM trajectory fitting with mgcv, and per-cluster pathway enrichment with clusterProfiler. Supports both R and Python alternatives at each step.

## Prerequisites
```bash
# R packages
install.packages(c('BiocManager'))
BiocManager::install(c('limma', 'DESeq2', 'Mfuzz', 'clusterProfiler', 'org.Hs.eg.db'))
install.packages(c('mgcv', 'cluster', 'splines'))

# Optional: rhythm detection
BiocManager::install('MetaCycle')

# Python alternative
pip install pandas numpy scipy statsmodels patsy scikit-learn tslearn pygam gseapy CosinorPy
```

**Input data:**
- Expression matrix (genes x samples) - normalized counts for limma, raw counts for DESeq2
- Sample metadata with time point column
- For circadian analysis: sampling over 24h+ cycles with 2-4h resolution

## Quick Start
Tell your AI agent what you want to do:
- "I have a time-course RNA-seq experiment with 6 time points - find temporal expression patterns"
- "Cluster my time-series gene expression data and find enriched pathways per cluster"
- "Run temporal differential expression and soft clustering on my microarray time course"
- "Check if any genes in my 48-hour time course show circadian rhythms"
- "Fit smooth trajectories to my temporal gene clusters and run GO enrichment"

## Example Prompts

### Temporal DE and Clustering
> "I have normalized counts from a 5-day differentiation experiment sampled every 12 hours. Find genes with significant temporal changes and cluster them into expression profiles."

> "Run DESeq2 LRT on my raw count matrix to identify time-dependent genes, then use Mfuzz to group them into soft clusters."

### Circadian Analysis
> "My experiment sampled liver tissue every 4 hours over 48 hours. Detect circadian genes and estimate their phases."

> "Run MetaCycle on my time-course data to find genes with 24-hour periodicity."

### Trajectory and Enrichment
> "Fit GAM curves to each temporal cluster and run GO enrichment to interpret the biological meaning of each pattern."

> "I have 8 Mfuzz clusters from my time-course experiment. Run pathway enrichment on each cluster and summarize the results."

## What the Agent Will Do
1. Load expression matrix and time metadata
2. Run temporal differential expression (limma splines or DESeq2 LRT depending on input type)
3. Filter significant temporal genes at FDR <0.05
4. Validate: check sufficient temporal genes detected (>100)
5. Perform Mfuzz soft clustering with automatic fuzzifier estimation
6. Validate: check no empty clusters, membership quality >0.5
7. Optionally run rhythm detection if circadian design detected
8. Fit GAM trajectories per cluster using mgcv
9. Run per-cluster GO/KEGG enrichment with clusterProfiler
10. Validate: check at least 3 clusters have significant enrichment terms
11. Export cluster assignments, trajectory fits, and enrichment tables

## Tips
- Use normalized counts (e.g., voom, rlog) for limma and Mfuzz; raw counts for DESeq2 LRT
- Mfuzz soft clustering assigns membership scores rather than hard labels; filter with membership >0.5 for core genes
- The fuzzifier parameter m controls cluster softness; mestimate() calculates it automatically
- Start with 4-12 clusters and refine using gap statistic or silhouette scores
- Rhythm detection only applies when sampling covers at least one full cycle (24h for circadian)
- GAM basis dimension k should not exceed (n_timepoints - 1); k=5 is a safe default
- Per-cluster enrichment reveals distinct biological processes; compare across clusters for a systems view
- For organisms other than human/mouse, swap org.Hs.eg.db for the appropriate OrgDb package

## Related Skills
- differential-expression/timeseries-de - Temporal DE methods
- temporal-genomics/temporal-clustering - Cluster analysis details
- temporal-genomics/circadian-rhythms - Rhythm detection details
- temporal-genomics/trajectory-modeling - GAM fitting details
- pathway-analysis/go-enrichment - Enrichment analysis details
- temporal-genomics/periodicity-detection - Unknown-period discovery
