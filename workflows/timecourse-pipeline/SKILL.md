---
name: bio-workflows-timecourse-pipeline
description: End-to-end bulk time-course analysis from an expression matrix to temporal gene modules and per-cluster pathway enrichment. Orchestrates temporal DE (limma splines or DESeq2 LRT), Mfuzz/tslearn soft clustering of expression-profile shapes, GAM trajectory fitting, per-cluster GO enrichment against a temporal-gene background, and an OPTIONAL circadian rhythm-detection branch (MetaCycle/CosinorPy) that runs only when the design covers >=2 full cycles with >=6-8 evenly spaced samples per cycle. Use when analyzing a bulk time-series expression experiment from any omics platform and deciding limma-splines vs DESeq2-LRT for temporal DE, soft vs hard clustering, whether the sampling design even licenses rhythm detection, and which background to use for enrichment. Not for single-cell pseudotime (see temporal-genomics/trajectory-modeling for the bulk-vs-pseudotime boundary) or unknown-period discovery (see temporal-genomics/periodicity-detection).
tool_type: mixed
primary_tool: Mfuzz
goal_approach_exempt: true
workflow: true
depends_on:
  - differential-expression/timeseries-de
  - temporal-genomics/temporal-clustering
  - temporal-genomics/circadian-rhythms
  - temporal-genomics/trajectory-modeling
  - pathway-analysis/go-enrichment
qc_checkpoints:
  - after_de: "Significant temporal genes >100 at FDR <0.05; model-fit residuals reasonable"
  - after_clustering: "Membership >0.5 for soft clustering; no empty clusters; k validated by silhouette/gap or bootstrap stability (typical 4-20)"
  - before_rhythm_detection: "GATE: design covers >=2 full cycles AND >=6-8 samples/cycle at ~even spacing AND collection order was randomized; else SKIP rhythm detection"
  - after_enrichment: "At least 3 clusters with significant GO terms at FDR <0.05; background = temporal genes, not the genome"
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, limma 3.58+, splines (R base), Mfuzz 2.62+, MetaCycle 1.2+, mgcv 1.9+, clusterProfiler 4.10+, CosinorPy 3.1+, tslearn 0.6+, pygam 0.9+, gseapy 1.2+, statsmodels 0.14+, patsy 1.0+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: the whole pipeline is dominated by the SAMPLING DESIGN, not the algorithm. Temporal DE, clustering, and trajectory fitting need genes pre-selected for temporal change and enough timepoints to resolve a shape; rhythm detection additionally requires >=2 full cycles and >=6-8 evenly spaced samples per cycle. A small p-value on 6 timepoints over a single cycle is not evidence of a rhythm.

# Time-Course Analysis Pipeline

**"Analyze my bulk time-course expression data end-to-end"** -> Orchestrate temporal differential expression, soft clustering of expression-profile shapes, GAM trajectory fitting, per-cluster pathway enrichment, and (only under a circadian sampling design) rhythm detection.
- Python: statsmodels/patsy spline F-test -> tslearn TimeSeriesKMeans -> pygam LinearGAM -> gseapy; CosinorPy for the optional rhythm branch
- R: limma splines or DESeq2 LRT -> Mfuzz -> mgcv -> clusterProfiler; MetaCycle for the optional rhythm branch

## Pipeline principles (read before running)

- Clustering, GAM fitting, and enrichment are DESCRIPTIVE steps DOWNSTREAM of the temporal-DE gene selection; they add no inference of their own. Cluster only the temporally variable genes, never the full matrix, or every method returns confident-looking clusters of noise.
- Do NOT test the clusters for the temporal signal used to select the genes (circular / double-dipping). Enrichment is valid only against an INDEPENDENT annotation (GO/KEGG), with the temporal genes as background.
- Rhythm detection is an OPTIONAL branch, not a routine default. It is licensed only by a circadian sampling design (>=2 full cycles, >=6-8 samples/cycle, roughly even spacing, randomized collection order). Under any other design it is SKIPPED, not run with a warning.
- Cluster number k and Mfuzz fuzzifier m are analyst CHOICES, not results; sweep and report the criterion.

## Pipeline Overview

```
Expression matrix + time metadata
    |
    v
[1. Temporal DE] ---------> limma splines / DESeq2 LRT / statsmodels spline F-test
    |                            (selects temporally variable genes; FDR <0.05)
    v
[2. Filter] --------------> Significant temporal genes only (clustering input)
    |
    v
[3. Soft Clustering] -----> Mfuzz (R) / tslearn TimeSeriesKMeans (Python) on z-scored profiles
    |                            +---> QC: membership >0.5, no empty clusters, sweep k
    |
    v
[4b. GAM Trajectory] -----> mgcv / pygam GAM on standardized cluster-mean profiles
    |
    v
[5. Pathway Enrichment] --> clusterProfiler / gseapy per cluster
    |                            background = temporal genes (NOT the genome)
    v
Temporal gene modules + enriched pathways + trajectory fits

    OPTIONAL branch, gated (NOT on the default path):
    IF design covers >=2 full cycles AND >=6-8 samples/cycle (even spacing)
       AND collection order was randomized:
        [4a. Rhythm Detection] --> MetaCycle meta2d / CosinorPy fit_group
                                   period/phase/amplitude + BH q across genes
    ELSE: SKIP (design does not license a rhythm test)
```

## Step 1: Temporal Differential Expression

### R (limma splines)

```r
library(limma)
library(splines)

expr <- as.matrix(read.csv('counts_normalized.csv', row.names = 1))
meta <- read.csv('metadata.csv')

# Natural cubic spline on time; df=3 is enough for most courses, raise to 4-5 for >10 timepoints
design <- model.matrix(~ ns(meta$time, df = 3))
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# Joint F-test on all spline coefficients = "expression changes over time"
temporal_results <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = 'F')
# topTable already returns adj.P.Val (BH-corrected); use it directly (not $FDR, which does not exist)
```

### R (DESeq2 LRT)

```r
library(DESeq2)

counts <- as.matrix(read.csv('raw_counts.csv', row.names = 1))
meta <- read.csv('metadata.csv')
meta$time <- factor(meta$time)

# LRT: full model (time as factor) vs reduced (intercept) = any between-timepoint change
dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ time)
dds <- DESeq(dds, test = 'LRT', reduced = ~ 1)
res <- results(dds)   # BH-adjusted padj
```

Choose limma-splines when the time axis is continuous and a smooth trend is expected (normalized/voom or vst input); choose DESeq2-LRT for raw counts and few, discrete timepoints treated as a factor.

### Python (statsmodels spline F-test)

```python
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from patsy import dmatrix
from scipy import stats

expr = pd.read_csv('counts_normalized.csv', index_col=0)
meta = pd.read_csv('metadata.csv')

spline_basis = dmatrix('bs(time, df=3)', data=meta, return_type='dataframe')
# patsy already includes an Intercept column; do NOT prepend another np.ones (that duplicates the
# intercept, inflates model df, and miscalibrates the F-test -> wrong temporal-DE FDR). Use it directly.
design_full = spline_basis.values                          # [Intercept, bs1, bs2, bs3] = 4 cols
design_reduced = np.ones((len(meta), 1))
df_diff = design_full.shape[1] - design_reduced.shape[1]   # 4 - 1 = 3
df_resid = len(meta) - design_full.shape[1]                # n - 4

pvals = []
for gene in expr.index:
    y = expr.loc[gene].values
    ss_full = np.sum((y - design_full @ np.linalg.lstsq(design_full, y, rcond=None)[0]) ** 2)
    ss_red = np.sum((y - design_reduced @ np.linalg.lstsq(design_reduced, y, rcond=None)[0]) ** 2)
    f_stat = ((ss_red - ss_full) / df_diff) / (ss_full / df_resid)
    pvals.append(1 - stats.f.cdf(f_stat, df_diff, df_resid))

# multipletests default is Holm-Sidak; force BH explicitly
_, fdr, _, _ = multipletests(pvals, method='fdr_bh')
temporal_genes = expr.index[fdr < 0.05].tolist()
```

### QC Checkpoint: Temporal DE

```r
sig_genes <- temporal_results[temporal_results$adj.P.Val < 0.05, ]
n_sig <- nrow(sig_genes)
message(sprintf('Significant temporal genes: %d', n_sig))
# <100: underpowered clustering; >10000: check batch/normalization before proceeding
if (n_sig < 100) message('WARNING: Few temporal genes. Check timepoint spacing or relax FDR.')
if (n_sig > 10000) message('WARNING: Many temporal genes. Inspect batch effects / normalization.')
```

## Step 2: Filter Significant Genes

```r
# FDR <0.05 standard; 0.1 acceptable for exploratory clustering only
sig_genes <- rownames(temporal_results[temporal_results$adj.P.Val < 0.05, ])
expr_sig <- expr[sig_genes, ]
message(sprintf('Genes passing FDR <0.05: %d', length(sig_genes)))
```

## Step 3: Soft Clustering (of expression-profile shapes)

Clustering groups genes by SHAPE, not magnitude, so profiles are z-scored per gene first; otherwise abundance dominates and high-expression genes cluster together regardless of dynamics. Cluster only `expr_sig` (the temporal genes), never the full matrix.

### R (Mfuzz)

```r
library(Mfuzz)

eset <- ExpressionSet(assayData = as.matrix(expr_sig))
eset <- standardise(eset)   # per-gene mean 0, sd 1: makes distance shape-based, not magnitude-based

# mestimate() implements Schwaemmle & Jensen (2010): the smallest m that stops fuzzy c-means
# from clustering RANDOMIZED data. It is dominated by the number of timepoints; inspect the
# returned value and the membership distribution rather than trusting it blindly (or hardcoding m=2).
m <- mestimate(eset)
message(sprintf('Estimated fuzzifier m = %.2f', m))

# k is a resolution CHOICE, not a result: start ~sqrt(n_genes/2), then sweep and validate (below)
n_clusters <- 8
cl <- mfuzz(eset, c = n_clusters, m = m)

# Membership >0.5 = core (confident) genes; lower to 0.3 only for exploratory overlap
core_genes <- acore(eset, cl, min.acore = 0.5)
```

### Python (tslearn)

```python
from tslearn.clustering import TimeSeriesKMeans

# Row-wise z-score: normalize each gene across its own timepoints (shape, not level)
expr_scaled = (expr_sig.values - expr_sig.values.mean(axis=1, keepdims=True)) / expr_sig.values.std(axis=1, keepdims=True)

# soft-DTW tolerates phase-shifted profiles Euclidean would split; gamma smooths the DTW geometry.
# Use plain 'euclidean' when absolute phase is biologically meaningful (morning vs evening genes).
model = TimeSeriesKMeans(n_clusters=8, metric='softdtw', metric_params={'gamma': 0.01},
                         max_iter=50, random_state=42)
labels = model.fit_predict(expr_scaled.reshape(expr_scaled.shape[0], expr_scaled.shape[1], 1))
```

### QC Checkpoint: Clustering

```r
library(cluster)
cluster_sizes <- table(cl$cluster)
print(cluster_sizes)
if (any(cluster_sizes == 0)) message('WARNING: Empty clusters. Reduce n_clusters.')

for (i in seq_along(core_genes)) {
    message(sprintf('Cluster %d: %d core genes (membership >0.5)', i, nrow(core_genes[[i]])))
}

# Silhouette on the z-scored profiles; triangulate k with a sweep, do not crown one index
sil <- silhouette(cl$cluster, dist(exprs(eset)))
message(sprintf('Mean silhouette: %.3f', mean(sil[, 3])))
```

## Step 4a: Rhythm Detection (OPTIONAL branch - GATED)

**This branch runs only when the sampling design licenses a rhythm test. Otherwise it is SKIPPED, not run with a warning.** Rhythm detection is not a routine step in a general time-course pipeline.

Hard precondition (all must hold), a design gate, not a soft aside:
- **>=2 full cycles** of the target period (48h+ for a 24h circadian rhythm). One cycle cannot distinguish an oscillation from a monotone trend or a single transient bump.
- **>=6-8 samples per cycle at roughly even spacing** (a lenient practical convention; Hughes 2017 recommends denser -- every 2h / >=12 per cycle -- and calls 4h intervals underpowered). Nyquist's 2/cycle is a mathematical floor with zero robustness to noise, no phase, no amplitude, no waveform shape.
- **Collection order randomized.** Harvest/collection order aliases directly onto circadian time: any drift (reagent lot, RIN, lane) is perfectly confounded with the rhythm axis and CANNOT be removed analytically. The only fix is design (randomize processing order).

Even when the gate passes, a rhythm found under a light-dark (LD) cycle may be light/feeding-DRIVEN masking, not endogenous clock output: diurnal != circadian. Endogeneity requires persistence under constant conditions (constant darkness, DD). Report ZT for entrained (LD) data, CT for free-running (DD) data.

```python
CIRCADIAN_DESIGN = False  # the un-computable precondition (circadian design + randomized order); set True only if it holds

n_cycles = (meta['time'].max() - meta['time'].min()) / 24.0
samples_per_cycle = meta['time'].nunique() / max(n_cycles, 1e-9)
gate_design = n_cycles >= 2 and samples_per_cycle >= 6   # the COMPUTABLE part of the gate
if CIRCADIAN_DESIGN and gate_design:
    pass  # run rhythm detection (CosinorPy/MetaCycle below)
elif not gate_design:
    print(f'Rhythm detection SKIPPED: inadequate design ({n_cycles:.1f} cycles, {samples_per_cycle:.1f}/cycle; need >=2 and >=6-8).')
else:
    print('Rhythm detection SKIPPED: design meets the cycle/sampling floor but CIRCADIAN_DESIGN is not set (randomized-order / circadian precondition unconfirmed).')
```

### R (MetaCycle)

```r
library(MetaCycle)
expr_for_meta <- expr_sig
colnames(expr_for_meta) <- meta$time
write.csv(expr_for_meta, 'expr_for_metacycle.csv')

# Circadian search window 20-28h; ARS/JTK require EVEN integer sampling and drop out silently
# (analysisStrategy='auto') on uneven/replicated data, leaving LS only.
meta2d('expr_for_metacycle.csv', filestyle = 'csv',
       minper = 20, maxper = 28,
       timepoints = sort(unique(meta$time)),
       outdir = 'metacycle_results')
# Filter on meta2d_BH.Q (BH FDR) AND meta2d_rAMP (relative amplitude); significance alone over-detects.
```

### Python (CosinorPy)

```python
from CosinorPy import cosinor   # note: import name is capitalized CosinorPy, not cosinorpy

# fit_group expects long-format columns 'x' (time), 'y' (expression), 'test' (gene id).
# period=24 for circadian; n_components=2 adds the 12h harmonic for non-sinusoidal shapes.
results = cosinor.fit_group(expr_long, period=24, n_components=1)

# fit_group ALREADY returns a BH-adjusted 'q' column across the fitted group; use it, do NOT
# threshold the raw per-gene 'p'. Add a RELATIVE-amplitude filter (fit_group has no rAMP column,
# so compute rAMP = amplitude/mesor): significance alone over-detects rhythms. rAMP>0.1 = >=10% of baseline.
results['rAMP'] = results['amplitude'] / results['mesor']
rhythmic = results[(results['q'] < 0.05) & (results['rAMP'] > 0.1)]
```

## Step 4b: GAM Trajectory Fitting

GAMs here summarize each cluster's temporal shape by fitting a penalized smooth to the STANDARDIZED cluster-mean profile. Because the input is a z-scored mean (not raw counts), a Gaussian family is appropriate; NB-family/offsets are needed only when fitting raw counts directly (see temporal-genomics/trajectory-modeling).

### R (mgcv)

```r
library(mgcv)

cluster_trajectories <- list()
for (cl_id in 1:n_clusters) {
    cl_genes <- names(cl$cluster[cl$cluster == cl_id])
    mean_profile <- colMeans(expr_sig[cl_genes, ])
    df_gam <- data.frame(time = meta$time, expr = mean_profile)

    # k is a flexibility CEILING (max basis dimension), NOT the number of knots/bends to fit.
    # REML (not GCV) then picks the wiggliness penalty; realized complexity is reported as edf.
    # Keep k < number of unique timepoints (identifiability); k=5 suits >=6-8 timepoints.
    gam_fit <- gam(expr ~ s(time, k = 5), data = df_gam, method = 'REML')

    cluster_trajectories[[cl_id]] <- list(fit = gam_fit,
                                          r_squared = summary(gam_fit)$r.sq,
                                          edf = summary(gam_fit)$edf)
    message(sprintf('Cluster %d: R^2 = %.3f, EDF = %.2f (edf~1 => linear; edf~k-1 => highly non-linear)',
                    cl_id, summary(gam_fit)$r.sq, summary(gam_fit)$edf))
}
```

### Python (pygam)

```python
from pygam import LinearGAM, s

for cl_id in range(n_clusters):
    mean_profile = expr_scaled[labels == cl_id].mean(axis=0)
    # n_splines is the basis-dimension ceiling (like mgcv k); the penalty picks realized wiggliness.
    gam = LinearGAM(s(0, n_splines=5)).fit(meta['time'].values.reshape(-1, 1), mean_profile)
    print(f'Cluster {cl_id}: GCV = {gam.statistics_["GCV"]:.4f}, edof = {gam.statistics_["edof"]:.2f}')
```

## Step 5: Per-Cluster Pathway Enrichment

The enrichment BACKGROUND (universe) must be the temporal genes that were clustered, NOT the whole genome. Genome background re-detects the generic biology of being a dynamic gene (the selection step), so every cluster lights up; temporal-gene background isolates what makes THIS shape distinct.

### R (clusterProfiler)

```r
library(clusterProfiler)
library(org.Hs.eg.db)

all_temporal_entrez <- bitr(rownames(expr_sig), fromType = 'SYMBOL', toType = 'ENTREZID',
                            OrgDb = org.Hs.eg.db)

enrichment_results <- list()
for (i in seq_along(core_genes)) {
    entrez <- bitr(core_genes[[i]]$NAME, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene = entrez$ENTREZID,
                    universe = all_temporal_entrez$ENTREZID,   # background = temporal genes
                    OrgDb = org.Hs.eg.db, ont = 'BP', pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    if (nrow(as.data.frame(ego)) > 0) {
        ego <- simplify(ego, cutoff = 0.7, by = 'p.adjust')   # collapse redundant parent-child GO terms
    }
    enrichment_results[[i]] <- ego
    message(sprintf('Cluster %d: %d significant GO terms', i, nrow(as.data.frame(ego))))
}
```

### Python (gseapy)

```python
import gseapy as gp

all_temporal_genes = list(expr_sig.index)   # background = temporal genes, not the genome

for cl_id in range(n_clusters):
    cl_genes = [g for g, l in zip(expr_sig.index, labels) if l == cl_id]
    # enrichr hits the Enrichr web API; pass background=temporal genes and outdir=None (no files).
    # For a strictly offline hypergeometric test with a custom background, gp.enrich(gene_sets=<gmt/dict>,
    # background=all_temporal_genes) computes it locally instead.
    enr = gp.enrichr(gene_list=cl_genes, gene_sets='GO_Biological_Process_2023',
                     organism='human', background=all_temporal_genes, outdir=None)
    sig_terms = enr.results[enr.results['Adjusted P-value'] < 0.05]
    print(f'Cluster {cl_id}: {len(sig_terms)} significant GO terms')
```

### QC Checkpoint: Enrichment

```r
clusters_with_terms <- sum(sapply(enrichment_results, function(x) nrow(as.data.frame(x)) > 0))
message(sprintf('Clusters with significant GO terms: %d / %d', clusters_with_terms, length(enrichment_results)))
if (clusters_with_terms < 3) message('WARNING: Few clusters enriched. Check gene ID mapping or thresholds.')
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| Temporal DE | Spline df | 3 (default); 4-5 for >10 timepoints |
| Temporal DE | FDR | 0.05 (standard); 0.1 exploratory clustering only |
| Clustering | fuzzifier m | Use mestimate(); inspect returned value + membership distribution |
| Clustering | n_clusters (k) | 4-20; a CHOICE, not a result; sweep + validate (silhouette/gap/bootstrap) |
| Clustering | min membership | 0.5 (core); 0.3 (exploratory) |
| Rhythm (gated) | design gate | >=2 cycles AND >=6-8 samples/cycle AND randomized order; else skip |
| Rhythm (gated) | period window | 20-28h circadian; filter on BH q AND rAMP |
| GAM | k (basis ceiling) | 5 for >=6-8 timepoints; keep k < #unique timepoints; REML picks the penalty |
| Enrichment | background | temporal genes (NOT genome); pvalueCutoff 0.05 |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Clusters look clean but mean nothing | Clustered the full matrix / did not z-score | Cluster only temporal-DE hits; standardise per gene first |
| Cluster enrichment lights up everywhere with generic terms | Genome used as enrichment background | Use temporal genes as universe/background |
| "Rhythmic" hits on a short/1-cycle design | Rhythm test run without the design gate | Enforce >=2 cycles + >=6-8 samples/cycle; else skip the branch |
| A 24h rhythm appears at ~12h or ~36h | Aliasing from sub-Nyquist sampling | Sample >=2x per target period; report the interval |
| Implausibly many rhythmic genes | Significance-only threshold; undetrended trend | Filter on rAMP/amplitude too; require >=2 cycles |
| CosinorPy import fails | Wrong import name | `from CosinorPy import cosinor` (capitalized) |
| MetaCycle "wrote files but read.csv fails" | Output is under `outdir/` as `meta2d_<infile>` | Read the actual emitted path; ARS/JTK silently drop on uneven sampling |
| gam.check k-index < 1 | Basis ceiling k too low (or residual autocorrelation) | Double k and refit; if edf barely moves, suspect autocorrelation, not k |
| GAM p=1e-30 over-trusted | Smooth-term p-values are approximate | Treat as categorical significant/not; apply BH across genes |

## Related Skills

- differential-expression/timeseries-de - Temporal DE methods (limma splines, DESeq2 LRT, maSigPro)
- temporal-genomics/temporal-clustering - Soft/hard clustering, k selection, DTW details
- temporal-genomics/circadian-rhythms - Single-condition rhythm detection, sampling design, phase/amplitude
- temporal-genomics/differential-rhythmicity - Comparing rhythms between conditions (gain/loss/phase/amplitude)
- temporal-genomics/trajectory-modeling - GAM fitting, k-vs-edf, bulk-vs-pseudotime boundary
- temporal-genomics/periodicity-detection - Unknown-period discovery (Lomb-Scargle, wavelets)
- pathway-analysis/go-enrichment - Enrichment, background choice, GO term simplification

## References

- Hughes ME, Abruzzi KC, Allada R, et al. 2017. Guidelines for Genome-Scale Analysis of Biological Rhythms. *J Biol Rhythms* 32(5):380-393. doi:10.1177/0748730417728663. (Recommends >=2 cycles and dense sampling -- at least every 2 h / >=12 timepoints per cycle, explicitly calling 4 h intervals underpowered -- plus biological replicates. The lenient 6-8/cycle floor used above is a stated practical convention, denser is better.)
- Wu G, Anafi RC, Hughes ME, Kornacker K, Hogenesch JB. 2016. MetaCycle: an integrated R package to evaluate periodicity in large scale data. *Bioinformatics* 32(21):3351-3353. doi:10.1093/bioinformatics/btw405. (meta2d integrating ARS/JTK/LS; output columns.)
- Moškon M. 2020. CosinorPy: a python package for cosinor-based rhythmometry. *BMC Bioinformatics* 21(1):485. doi:10.1186/s12859-020-03830-w. (fit_group / cosinor rhythmometry.)
- Futschik ME, Carlisle B. 2005. Noise-robust soft clustering of gene expression time-course data. *J Bioinform Comput Biol* 3(4):965-988. doi:10.1142/S0219720005001375. (Fuzzy c-means noise-robustness rationale for Mfuzz.)
- Schwämmle V, Jensen ON. 2010. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. *Bioinformatics* 26(22):2841-2848. doi:10.1093/bioinformatics/btq534. (The mestimate() fuzzifier estimator.)
- Wood SN. 2011. Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. *J R Stat Soc B* 73(1):3-36. doi:10.1111/j.1467-9868.2010.00749.x. (Why REML over GCV for the smoothing penalty.)
- Laloum D, Robinson-Rechavi M. 2020. Methods detecting rhythmic gene expression are biologically relevant only for strong signal. *PLoS Comput Biol* 16(3):e1007666. doi:10.1371/journal.pcbi.1007666. (Amplitude filtering; significance alone over-detects rhythms.)
