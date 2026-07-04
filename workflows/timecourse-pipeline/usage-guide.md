# Time-Course Pipeline - Usage Guide

## Overview
Complete bulk time-course expression analysis from an expression matrix to temporal gene modules and pathway enrichment. The result is dominated by the SAMPLING DESIGN, not the algorithm: clustering, GAM fitting, and enrichment are descriptive steps DOWNSTREAM of temporal-DE gene selection, so they run on the temporally variable genes only, never the full matrix. The pipeline covers temporal differential expression (limma splines or DESeq2 LRT), soft clustering of expression-profile shapes (Mfuzz or tslearn), GAM trajectory fitting (mgcv or pygam), and per-cluster GO enrichment against a temporal-gene background (clusterProfiler or gseapy). Circadian rhythm detection (MetaCycle or CosinorPy) is an OPTIONAL, gated branch that runs only when the design covers at least two full cycles with 6-8 evenly spaced samples per cycle and a randomized collection order; under any other design it is skipped, because a rhythm test on a short or aliased series is uninterpretable regardless of p-value.

## Prerequisites
```bash
# R packages
install.packages('BiocManager')
BiocManager::install(c('limma', 'DESeq2', 'Mfuzz', 'clusterProfiler', 'org.Hs.eg.db'))
install.packages(c('mgcv', 'cluster'))   # splines ships with base R

# Optional: gated rhythm-detection branch
install.packages('MetaCycle')

# Python alternative
pip install pandas numpy scipy statsmodels patsy scikit-learn tslearn pygam gseapy CosinorPy
```

**Input data:**
- Expression matrix (genes x samples) - normalized/voom/vst counts for limma and clustering, raw counts for DESeq2 LRT
- Sample metadata with a numeric time column
- For the optional rhythm branch only: sampling over >=2 full cycles with >=6-8 samples/cycle at roughly even spacing, and randomized collection order

## Quick Start
Tell your AI agent what you want to do:
- "I have a 6-timepoint RNA-seq time course - find temporal patterns and cluster them"
- "Cluster my time-series gene expression and find enriched pathways per cluster"
- "Run temporal differential expression and soft clustering on my microarray time course"
- "Fit smooth trajectories to my temporal gene clusters and run GO enrichment"
- "My liver samples span 48h every 4h with randomized collection - is anything circadian?"

## Example Prompts

### Temporal DE and Clustering
> "I have normalized counts from a 5-day differentiation experiment sampled every 12 hours. Find genes with significant temporal changes and cluster them into expression profiles."

> "Run DESeq2 LRT on my raw count matrix to identify time-dependent genes, then use Mfuzz to group them into soft clusters and report the membership fraction retained."

### Circadian Analysis (only for a circadian design)
> "My experiment sampled liver tissue every 4 hours over 48 hours with randomized processing order. Check the design gate, then detect circadian genes and estimate their phases."

> "I sampled every 12 hours for 24 hours - can I test for circadian rhythms?" (Expected answer: no; the design covers one cycle at two samples per cycle, below the gate.)

### Trajectory and Enrichment
> "Fit GAM curves to each temporal cluster and run GO enrichment against the temporal genes as background to interpret each pattern."

> "I have 8 Mfuzz clusters from my time-course experiment. Run pathway enrichment per cluster and summarize the distinct biological programs."

## What the Agent Will Do
1. Load the expression matrix and time metadata.
2. Run temporal differential expression (limma splines or DESeq2 LRT depending on input type).
3. Filter to significant temporal genes at FDR <0.05 (the clustering input).
4. Validate that enough temporal genes were detected (>100).
5. Soft-cluster the z-scored profiles (Mfuzz or tslearn) with fuzzifier estimation and a k sweep.
6. Validate no empty clusters and adequate membership (>0.5 core genes).
7. Evaluate the rhythm-detection GATE; run MetaCycle/CosinorPy ONLY if the circadian design holds, otherwise skip.
8. Fit GAM trajectories per cluster on standardized cluster-mean profiles (REML smoothing).
9. Run per-cluster GO enrichment with the temporal genes as background (not the genome).
10. Validate that at least 3 clusters carry significant enrichment terms.
11. Export cluster assignments, trajectory fits, and enrichment tables.

## Tips
- Use normalized counts (voom, rlog, vst) for limma and clustering; raw counts for DESeq2 LRT.
- Cluster only the temporal-DE hits; clustering the full matrix returns confident-looking clusters of noise once profiles are z-scored.
- Mfuzz assigns membership scores, not hard labels; filter with membership >0.5 for core genes and report the retained fraction.
- The fuzzifier m from mestimate() is dominated by the number of timepoints; inspect the returned value rather than hardcoding m=2.
- Cluster number k is a resolution choice, not a result; sweep it and report the criterion (silhouette, gap, or bootstrap stability).
- Rhythm detection is gated: it needs >=2 full cycles, >=6-8 samples/cycle at even spacing, and randomized collection order. A rhythm found under a light-dark cycle may be light-driven masking (diurnal), not endogenous (circadian).
- For GAMs, k is a flexibility ceiling, not the number of bends; keep k < the number of unique timepoints and let REML pick the penalty. Read edf (not k) for realized complexity.
- Always set the enrichment background to the temporal genes, not the full genome, or every cluster lights up for generic dynamic-gene biology.
- Use simplify() on GO results to collapse redundant parent-child terms before interpreting clusters.
- For organisms other than human, swap org.Hs.eg.db for the appropriate OrgDb package.

## Related Skills
- differential-expression/timeseries-de - Temporal DE methods
- temporal-genomics/temporal-clustering - Cluster analysis details
- temporal-genomics/circadian-rhythms - Rhythm detection and sampling design
- temporal-genomics/differential-rhythmicity - Comparing rhythms between conditions
- temporal-genomics/trajectory-modeling - GAM fitting details
- temporal-genomics/periodicity-detection - Unknown-period discovery
- pathway-analysis/go-enrichment - Enrichment analysis details
