---
name: bio-atac-seq-motif-deviation
description: Analyze TF motif accessibility variability across samples or single cells using chromVAR. Use when identifying TF motifs whose accessibility correlates with conditions, computing per-sample motif z-scores after matched background correction, comparing to ArchR / Signac equivalents, or distinguishing motif-accessibility signal from per-site footprinting.
tool_type: r
primary_tool: chromVAR
---

## Version Compatibility

Reference examples tested with: chromVAR 1.24+, motifmatchr 1.24+, JASPAR2024 0.99+, TFBSTools 1.40+, BSgenome.Hsapiens.UCSC.hg38 1.4+, SummarizedExperiment 1.32+, limma 3.58+, ggplot2 3.5+, Matrix 1.6+, ArchR 1.0.2+, Signac 1.13+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# Motif Deviation (chromVAR)

**"Which TF motifs explain accessibility variation across my samples or cells?"** -> Compute per-sample (or per-cell) deviation z-scores: how many standard deviations above expectation each TF motif's accessibility falls, controlling for GC content and overall accessibility via matched background peak sets.

- R: `chromVAR::computeDeviations(counts, motifs)` -> per-sample z-scores
- R: `chromVAR::computeVariability(dev)` -> per-motif variance ranking
- Single-cell alternative: `Signac::RunChromVAR()` (wrapper with matched defaults) or `ArchR::addDeviationsMatrix()`

chromVAR answers a different question than footprinting: footprinting asks "is this specific motif site bound?", chromVAR asks "do peaks containing this motif have systematically more or less accessibility than expected?" The two are complementary.

## What chromVAR Computes

For each (motif, sample) pair:
- **Raw deviation** = Sum of accessibility at peaks containing the motif - expected from a matched-GC, matched-accessibility background.
- **Bias-corrected deviation** = Raw deviation / SD of background deviations.
- **Z-score** = (corrected deviation - mean across cells) / SD across cells. Reported as the principal output.

Z-scores are signed: positive = motif more accessible in this sample than population average; negative = less. Magnitudes 2-5 are typical for biologically interesting motifs; >5 indicates strong covariation with sample state.

## Algorithmic Taxonomy

| Tool | Input | Background | Output | Best for | Fails when |
|------|-------|------------|--------|----------|------------|
| chromVAR | Peak count matrix + motif annotations | Matched GC + accessibility (50 peaks per match by default) | Per-sample motif z-score | Bulk + single-cell (sparse-aware); cross-sample variability | < 1500 reads/sample (bulk) or < 500 cells/cluster (sc); too few peaks (< 5000) |
| Signac::RunChromVAR | Seurat single-cell ATAC object | Same as chromVAR (delegated) | Motif assay in Seurat object | Standard single-cell workflows in Seurat ecosystem | Same as chromVAR; needs Seurat object setup |
| ArchR::addDeviationsMatrix | ArrowFile + tile/peak matrix | ArchR's getBgdPeaks (matched on GC + log accessibility) | Per-cell deviation matrix in ArchR project | ArchR ecosystem; faster on large scATAC | ArchR-specific format; not portable to chromVAR objects |
| chromVAR (multinomial) | Peak counts + motif annotations | Same | Same z-scores via different model | When read-depth is highly variable across samples | Original poisson-based variant; multinomial in newer versions |
| Signac::FindMarkers (with motifs as features) | Motif accessibility matrix from RunChromVAR | Per-cell-cluster | Differential motifs per cluster | Cluster-level differential | Differential test must be on z-scores; raw counts will mislead |
| TF activity inference (DecoupleR / SCENIC+) | Gene expression + motif activity | Multi-modal | TF activity score | Multi-omics integration | Requires paired RNA-seq; chromVAR alone is insufficient |

Methodology evolves; verify against current chromVAR (Schep 2017), ArchR (Granja 2021), and Signac (Stuart 2021) benchmarks before locking pipelines.

## chromVAR vs Footprinting -- Different Questions

| Question | Tool |
|----------|------|
| Does the bulk pattern of motif-containing peaks vary with condition? | chromVAR |
| Is THIS specific motif site bound by a TF? | TOBIAS / HINT-ATAC |
| Per-cell TF activity in scATAC | chromVAR (via Signac/ArchR) |
| Per-cell TF binding at specific sites | scprinter |
| TF activity correlated with gene expression | chromVAR + co-expression OR SCENIC+ |
| Which TF families distinguish my cell clusters? | chromVAR per-cluster z-scores |
| Differential bound vs unbound between conditions | TOBIAS BINDetect |

chromVAR is fundamentally a *summary statistic* over many motif sites. Footprinting is per-site classification. Use chromVAR when motif site count >> 100; use footprinting when specific sites matter.

## Per-Tool Failure Modes

### chromVAR -- Too few peaks or too few reads

**Trigger:** ATAC peakset < 5000 peaks; per-sample read depth < 1500 in peaks.

**Mechanism:** chromVAR's background sampling requires enough peaks to find matched GC + accessibility partners. Sparse sampling at low peak count creates correlated null distributions, inflating both positive and negative z-scores.

**Symptom:** Variability scores all > 5 (suspiciously high); top variable motifs are dominated by AT-rich or GC-rich sequences regardless of biology.

**Fix:** Verify peakset is at full ATAC scale (typically 50k-200k peaks). For sc ATAC, aggregate cells to clusters of >= 500 cells before running.

### chromVAR background peaks -- Default is good, custom requires care

**Trigger:** Calling `getBackgroundPeaks()` with non-default `niterations` or `bias`.

**Mechanism:** Default 50 iterations × 10 bgd peaks per peak generates 500 background sets. Reducing `niterations` increases noise; increasing slows linearly without much accuracy gain.

**Symptom:** Custom backgrounds inflate variability when niterations < 30.

**Fix:** Stick to defaults unless benchmarking. If running on huge cell counts, test on subsample first.

### chromVAR on broadly accessible cell types -- Z-scores compressed

**Trigger:** Multiple cell types in dataset have very different overall accessibility magnitudes.

**Mechanism:** chromVAR's correction normalizes for total accessibility; cell types with high background accessibility have compressed z-scores even if their motif-specific signal is strong.

**Symptom:** PCA on z-scores does not separate cell types as cleanly as raw counts.

**Fix:** Run chromVAR per-cell-type-cluster (separate runs) when global accessibility differs by > 5x. Alternatively use ArchR's per-cluster background.

### chromVAR on bulk samples without enough variation -- All z-scores near zero

**Trigger:** All bulk samples are technical replicates or very similar.

**Mechanism:** Z-scores normalize across the sample population; if there is no across-sample variability, all z-scores collapse to zero.

**Symptom:** Variability ranking is unstable across runs; top motifs change.

**Fix:** chromVAR is designed for variability; if the dataset has only one biological condition replicated, use footprinting or differential accessibility instead. chromVAR needs 6+ samples with biological variation to be informative.

### Signac::RunChromVAR -- Motif matching mismatch

**Trigger:** Motif assay added before peak set finalized; peak coordinates change.

**Mechanism:** RunChromVAR matches motifs to peaks at the time it's called; if peaks change downstream (e.g., after merge), the motif annotations become stale.

**Symptom:** Some peaks have NA motif annotations; deviation matrix has missing entries.

**Fix:** Run `AddMotifs()` -> `RunChromVAR()` AFTER finalizing peakset. Re-run if peaks change.

### ArchR::addDeviationsMatrix -- TileMatrix vs PeakMatrix

**Trigger:** Calling on tile matrix when peak matrix is more appropriate.

**Mechanism:** ArchR can compute deviations on either tiles (regular bins) or peaks. Peaks are biologically meaningful; tiles add noise from intergenic background.

**Fix:** Use `matrixName='PeakMatrix'` after addReproduciblePeakSet. Tile-based deviations are mainly for embedding, not biology.

## Decision Tree by Setting

| Setting | Workflow |
|---------|---------|
| Bulk, 6+ samples, condition contrast | chromVAR + limma differential on z-scores; rank by FDR |
| Bulk, 3-5 samples | chromVAR; report variability ranking; differential underpowered |
| scATAC, Signac ecosystem | Signac AddMotifs + RunChromVAR; FindMarkers on motif assay |
| scATAC, ArchR ecosystem | ArchR addPeakMatrix + addDeviationsMatrix + getMarkerFeatures |
| Multimodal scATAC + scRNA | chromVAR + paired DE; consider SCENIC+ for TF -> target inference |
| Plant / non-model organism | chromVAR with custom motif PFM (from CIS-BP); custom BSgenome |
| Time-course bulk (5+ time points) | chromVAR z-scores -> spline regression on time; identify motifs with non-monotone trajectories |

## chromVAR Workflow (Bulk)

```r
library(chromVAR); library(motifmatchr); library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2024); library(TFBSTools); library(SummarizedExperiment)

peaks <- rtracklayer::import('consensus_peaks.bed')           # GRanges
counts <- as.matrix(read.delim('peak_counts.tsv', row.names=1))    # rows = peaks, cols = samples
se <- SummarizedExperiment(assays=list(counts=counts), rowRanges=peaks)
se <- addGCBias(se, genome=BSgenome.Hsapiens.UCSC.hg38)

# Filter: depth >= 1500 reads/sample, FRiP >= 0.15, peak in >= 10% samples with >= 10 reads
se <- filterSamples(se, min_depth=1500, min_in_peaks=0.15, shiny=FALSE)
se <- filterPeaks(se, non_overlapping=TRUE, min_count=10, n_samples_frac=0.1)

# Motifs: JASPAR vertebrate CORE (default for human/mouse)
# JASPAR2024 + TFBSTools incompatibility (TFBSTools issue #39): getMatrixSet does not dispatch on the
# JASPAR2024 object directly. Open the SQLite handle and pass that to getMatrixSet instead.
library(RSQLite)
jaspar2024 <- JASPAR2024::JASPAR2024()
sq <- dbConnect(SQLite(), db(jaspar2024))
pfm <- getMatrixSet(sq, opts=list(collection='CORE', tax_group='vertebrates'))
# JASPAR2020 (older) accepts the package object directly: getMatrixSet(JASPAR2020, opts=...)
motif_ix <- matchMotifs(pfm, se, genome=BSgenome.Hsapiens.UCSC.hg38, p.cutoff=5e-05)

# Background peaks: matched GC + accessibility (default 50 iterations is fine)
bg <- getBackgroundPeaks(object=se, niterations=50)

# Compute deviations
dev <- computeDeviations(object=se, annotations=motif_ix, background_peaks=bg)
zscores <- deviations(dev)              # motif x sample matrix of z-scores
variability <- computeVariability(dev) # per-motif variability ranking
```

## Differential Motif Activity (limma on z-scores)

```r
library(limma)
groups <- factor(colData(se)$condition, levels=c('control', 'treated'))
design <- model.matrix(~groups)
fit <- lmFit(zscores, design); fit <- eBayes(fit)
diff_motifs <- topTable(fit, coef=2, number=Inf, p.value=0.05)   # adj.P.Val column
```

Use `adj.P.Val` (limma's BH FDR), not `FDR` (which limma does not return). `logFC` is the difference in z-scores between groups; magnitudes 0.5-2 typical.

## chromVAR for Single-Cell ATAC (Signac)

```r
library(Signac); library(Seurat); library(JASPAR2024); library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38); library(RSQLite)

# Assume `seurat_obj` has an ATAC assay with consensus peaks
# JASPAR2024 + TFBSTools workaround (see TFBSTools issue #39):
jaspar2024 <- JASPAR2024::JASPAR2024()
sq <- dbConnect(SQLite(), db(jaspar2024))
pfm <- getMatrixSet(sq, opts=list(collection='CORE', tax_group='vertebrates'))
seurat_obj <- AddMotifs(seurat_obj, genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
seurat_obj <- RunChromVAR(seurat_obj,
                          genome=BSgenome.Hsapiens.UCSC.hg38,
                          new.assay.name='chromvar')
DefaultAssay(seurat_obj) <- 'chromvar'

# Per-cluster differential motifs
markers <- FindAllMarkers(seurat_obj, only.pos=TRUE, mean.fxn=rowMeans, fc.name='avg_diff')
```

`mean.fxn=rowMeans` is required for z-score-style data; the default (`fc.name='avg_log2FC'`) makes no sense here.

## chromVAR for Single-Cell ATAC (ArchR)

```r
library(ArchR)
proj <- addReproduciblePeakSet(proj, groupBy='Clusters', pathToMacs2='/path/macs2')
proj <- addPeakMatrix(proj)
proj <- addMotifAnnotations(proj, motifSet='cisbp', name='Motif')
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(proj, peakAnnotation='Motif')

# Per-cluster deviation summary
markersMotifs <- getMarkerFeatures(proj, useMatrix='MotifMatrix',
                                   groupBy='Clusters', useSeqnames='z')
```

ArchR uses `cisbp` by default (CIS-BP database, ~5000 motifs); switch to `JASPAR2020` for fewer, more curated motifs.

## Reconciling chromVAR vs ArchR vs Signac

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Top variable motifs disagree | Different motif databases (JASPAR vs CIS-BP) | Re-run with matched motif set |
| Z-scores correlate but magnitudes differ | Different background sampling | Inspect per-tool background; defaults are similar but not identical |
| Signac chromvar assay has NA values | Motifs added after peakset finalized | Re-run AddMotifs + RunChromVAR after peaks are stable |
| ArchR per-cluster signature differs from Signac | Different clustering; different cell membership | Standardize clustering before comparison |

**Operational rule:** chromVAR z-scores are tool-specific. For cross-study comparison, recompute on the same peakset with the same motif database; do not use stored z-scores from heterogeneous sources directly.

## Variability Score Interpretation

| Variability | Z-score range typical | Interpretation |
|-------------|----------------------|----------------|
| < 1 | -1 to +1 | Motif activity ~constant; not biologically variable |
| 1-2 | -2 to +2 | Modest variation; condition-driven possible |
| 2-5 | -3 to +5 | Strong cross-sample / cross-cluster variability; biologically interesting |
| > 5 | -5 to +10 | Major driver of cell-state differences; flagship hits |

Variability is the across-sample variance of z-scores; it ranks motifs without requiring condition labels. For unsupervised TF discovery (e.g., trajectory analysis) variability is the primary metric.

## Background Peak Matching Mathematics

**Trigger:** Tuning chromVAR's `getBackgroundPeaks` parameters; benchmarking against published results.

**Mechanism:** chromVAR matches each foreground peak to background peaks by GC content + total accessibility, in 25 GC bins x 25 accessibility bins (default). For each foreground peak, the algorithm samples k_iterations (default 50) replacement peaks from the matching bin. Variance across these matched samples becomes the null reference.

**Threshold tuning:**
- **bins=25** (default): each bin is ~4% GC range; works for >= 5,000 peaks. For very small peaksets (< 2,000), drop to `bins=10` (~10% GC range) to avoid empty bins.
- **niterations=50** (default): with 25 x 25 = 625 bins, each background sample uses ~50 peaks per match. Reducing below 30 inflates noise; increasing above 100 yields diminishing returns.

For non-canonical genomes (mouse mm10 with different GC distribution), consider rebuilding bins manually with `quantile()` to ensure equal-sized bins.

## chromVAR vs scBasset for Single-Cell

| Tool | Approach | Best for | Limitation |
|------|---------|----------|------------|
| chromVAR | Matched-background z-score per motif | Standard sc workflow; integrated in Signac/ArchR | Linear; no sequence context beyond motif PWM |
| scBasset (Yuan & Kelley 2022) | Sequence CNN with per-cell projection | Higher cluster-discrimination accuracy than chromVAR; predicts cell states from sequence | Newer; ecosystem smaller; needs >= 100 cells per cluster for stable projection |
| EnFormer-derived TF activity | Long-context Transformer | Cross-cell-type TF activity prediction; distal regulation | Pre-trained models cell-type-specific |
| DecoupleR ULM/MLM (Badia-i-Mompel 2022) | Multi-method consensus TF activity scoring | Multi-omics integration; aggregation across motif databases | Requires careful cell-x-motif input matrix |

For high-stakes per-cell TF activity, run chromVAR + scBasset and report the intersection. See atac-seq/deep-learning-atac for scBasset details.

## DecoupleR Multi-Method TF Activity

```python
import decoupler as dc
# adata: AnnData with motif_x_cell deviation matrix as input
acts_ulm = dc.run_ulm(mat=adata.obsm['chromvar'], net=collectri_net,
                      source='source', target='target')
acts_mlm = dc.run_mlm(mat=adata.obsm['chromvar'], net=collectri_net,
                      source='source', target='target')
acts_consensus = dc.run_consensus(mat=adata.obsm['chromvar'], net=collectri_net)
```

DecoupleR aggregates multiple TF-activity inference methods (ULM, MLM, viper, GSVA, etc.). The consensus output is more robust than any single method to motif database biases.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error in addGCBias`: missing `seqlengths` | GRanges object lacks chrom sizes | Use `seqlengths(peaks) <- seqlengths(genome)` first |
| All z-scores near zero | Too few samples or too little variation | chromVAR requires biological variation; use footprinting or differential instead |
| `getBackgroundPeaks` slow | Default niterations and large peakset | Default is fine; do not reduce iterations below 30 |
| Differential motifs all significant | FDR not applied; or compared identical samples | Apply BH correction; verify groups are correct |
| Signac chromvar assay all zero | RunChromVAR called before peakset | Re-run after AddMotifs and peakset stability |
| FindAllMarkers reports `avg_log2FC` for chromvar | Default fc method incorrect for z-scores | Use `mean.fxn=rowMeans` and `fc.name='avg_diff'` |
| z-score interpretation flipped | Sign of contrast reversed | Verify factor level order; first level is reference |
| ArchR `cisbp` vs `JASPAR2020` results differ | Different motif databases | Choose one and report the choice |

## References

- Schep AN et al 2017 Nat Methods 14:975 (chromVAR)
- Granja JM et al 2021 Nat Genet 53:403 (ArchR)
- Stuart T et al 2021 Nat Methods 18:1333 (Signac)
- Aibar S et al 2017 Nat Methods 14:1083 (SCENIC; downstream TF target inference)
- Bravo Gonzalez-Blas C et al 2023 Nat Methods 20:1355 (SCENIC+)
- Castro-Mondragon JA et al 2022 NAR 50:D165 (JASPAR 2022)
- Weirauch MT et al 2014 Cell 158:1431 (CIS-BP)
- Vorontsov IE et al 2024 NAR 52:D116 (HOCOMOCO v12)

## Related Skills

- atac-seq/footprinting - Per-site TF binding (different question)
- atac-seq/differential-accessibility - Peak-level DA (alternative approach)
- atac-seq/single-cell-atac - sc workflow integration with Signac/ArchR
- atac-seq/co-accessibility - Cis-regulatory connections
- atac-seq/deep-learning-atac - scBasset / EnFormer alternative
- gene-regulatory-networks/scenic-regulons - Downstream TF -> target inference
- chip-seq/motif-analysis - Alternative motif-enrichment approaches
- single-cell/clustering - Inputs for per-cluster motif activity
