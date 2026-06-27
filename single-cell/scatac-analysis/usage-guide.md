# scATAC-seq Analysis - Usage Guide

## Overview

Single-cell ATAC-seq measures chromatin accessibility per cell, revealing cell-type-specific regulatory elements and TF activity. This skill centers on Signac (R/Seurat) with ArchR (R, on-disk, large data) and SnapATAC2 (Python, >1M cells) as alternatives. It treats the central epistemics of the data: a zero is ambiguous, the matrix is near-binary by sampling not biology, binarization is disfavored, gene activity is a weak proxy, and the peak set is circular with clustering.

## Prerequisites

```r
install.packages('Signac')
BiocManager::install(c('EnsDb.Hsapiens.v86', 'BSgenome.Hsapiens.UCSC.hg38', 'chromVAR', 'motifmatchr', 'JASPAR2020', 'TFBSTools', 'scDblFinder'))
devtools::install_github('GreenleafLab/ArchR')   # large on-disk workflows
```

```bash
pip install snapatac2 scanpy   # Python alternative, >1M cells (the SnapATAC2 example imports scanpy)
```

MACS2 or MACS3 must be on PATH for peak calling. Peaks, fragments, the EnsDb annotation, and the BSgenome must all be on the same genome build; a mismatch is silently wrong (no crash).

## Quick Start

Tell your AI agent what you want to do:
- "Process my 10X scATAC fragments through QC, LSI, and clustering"
- "Diagnose which LSI components track sequencing depth and drop them"
- "Call consensus peaks per cell type and quantify a peak matrix"
- "Run chromVAR with a GC-matched background and rank TFs by z-score"
- "Detect homotypic and heterotypic doublets and combine them"
- "Pick a framework for 2 million cells"

## Example Prompts

### Data Loading and QC
> "Load my 10X filtered_peak_bc_matrix.h5 with the fragments file and add TSS enrichment and nucleosome signal"
> "Filter cells from the joint distribution of TSS enrichment and fragment count, not copied thresholds"

### Dimensionality Reduction
> "Run TF-IDF and SVD, then show DepthCor and drop only the components that track depth"
> "Use ArchR iterative LSI so rare populations are not lost in a single clustering pass"

### Peaks and Differential Accessibility
> "Call peaks per cluster on pseudobulk and merge to a fixed-width consensus set"
> "Find differentially accessible peaks with a logistic-regression test and fragment count as a latent variable"

### Motifs
> "Run chromVAR with GC-matched backgrounds and rank TFs by z-score"
> "This motif is enriched; confirm the actual TF with multiome RNA before claiming it drives the program"

### Doublets and Integration
> "Run AMULET for homotypic doublets and ArchR doublet scores for heterotypic, then combine"
> "Transfer cell-type labels from my scRNA-seq reference onto the ATAC cells"

## What the Agent Will Do

1. Load fragments and build a ChromatinAssay (Signac) or Arrow files (ArchR).
2. Compute QC on chromatin signal (TSS enrichment, nucleosome signal, FRiP) and filter from the data's own distributions.
3. Run TF-IDF + SVD, then diagnose and drop depth-correlated LSI components via DepthCor.
4. Cluster on the retained LSI dimensions and embed with UMAP.
5. Call peaks per cluster on pseudobulk and merge into a consensus set, re-quantifying a peak matrix.
6. Score chromVAR motif deviations against GC-matched backgrounds, ranking with z-scores.
7. Detect homotypic (AMULET) and heterotypic (ArchR/scDblFinder) doublets and combine.
8. Build gene-activity scores for cluster-level annotation and scRNA integration only.

## ArchR Workflow (Large, On-Disk)

ArchR keeps data in chunked HDF5 Arrow files and is the R choice for ~1M cells. The flow: `createArrowFiles(filterTSS=4, filterFrags=1000, addTileMat=TRUE, addGeneScoreMat=TRUE)` -> `ArchRProject()` -> `filterDoublets()` -> `addIterativeLSI(useMatrix='TileMatrix')` -> `addClusters()` -> `addUMAP()` -> `addGroupCoverages()` -> `addReproduciblePeakSet(pathToMacs2=...)` -> `addPeakMatrix()`. Motifs: `addMotifAnnotations(motifSet='cisbp')` + `peakAnnoEnrichment()` (set `background='bgdPeaks'` for a GC-fair test; default `'all'` is not GC-matched). Deviations: `addBgdPeaks()` + `addDeviationsMatrix()`. ArchR `filterTSS=4`/`filterFrags=1000` are human-tuned defaults whose numeric value depends on the TSS set and are not transferable. HDF5 file-locking fails on networked filesystems.

## SnapATAC2 Workflow (Python, scverse)

SnapATAC2 stores data in backed AnnData and scales matrix-free past 1M cells. The flow: `pp.import_data()` -> `metrics.tsse()` -> `pp.add_tile_matrix(bin_size=500, counting_strategy='paired-insertion')` (operationalizes the anti-binarization evidence) -> `pp.select_features()` -> `tl.spectral()` (graph-Laplacian, SD-weighted, sidesteps the LSI "drop component 1" step) -> `pp.knn()` -> `tl.leiden()` -> `tl.macs3()` + `tl.merge_peaks()`. Use `distance_metric='cosine'` for spectral; `'jaccard'` without subsampling is a memory trap.

## Decision Guidance

### Framework
- Signac for Seurat-integrated multimodal (WNN), up to ~10^5 cells.
- ArchR for large R workflows (~1M cells) and the integrated peak/GRN/trajectory suite.
- SnapATAC2 for >1M cells and the scverse/scvi-tools Python stack.

### Binarize or not
- Do not binarize; model fragment counts (paired-insertion / PoissonVI). The practical gain concentrates in deeply sequenced data where the count=2 tier is populated.

### Embedding
- LSI (TF-IDF + SVD) with depth-correlated components dropped by DepthCor diagnostic.
- SnapATAC2 spectral embedding (SD-weighted) avoids the manual component-drop step.
- cisTopic/LDA when interpretable cis-regulatory topics are wanted; PeakVI/PoissonVI for deep generative with explicit depth modeling.

## Key Differences from scRNA-seq

| Aspect | scRNA-seq | scATAC-seq |
|--------|-----------|------------|
| Features | ~20,000 genes | ~100,000+ peaks |
| Per-cell non-zeros | 10-45% of genes | ~1-10% of peaks (near-binary by sampling) |
| Normalization | log-normalize / SCTransform | TF-IDF (model fragment counts) |
| Dim reduction | PCA | LSI (drop depth components) or spectral |
| Cell-type markers | gene expression | accessibility + motif activity (gene activity only a proxy) |

## Tips

- **A zero is ambiguous** - closed vs uncaptured; interpret accessibility at cluster/pseudobulk level, not per single cell.
- **Do not binarize** - model fragment counts; the count=2 tier carries real signal that grows with depth.
- **Diagnose the depth component** - run DepthCor and drop components that track depth, not blindly LSI_1.
- **Gene activity is cluster-level only** - repressed/bivalent promoters and misassigned distal enhancers make it a weak single-cell proxy.
- **Peak set is circular with clustering** - call peaks per cluster iteratively; never test DA on a peak set called from the same comparison.
- **chromVAR needs GC-matched backgrounds** - report z-scores for cross-TF ranking, never raw deviations.
- **Motif is not a TF** - a motif implicates a family; confirm with TF expression or a footprint before causal claims.
- **Run both doublet strategies** - AMULET (homotypic) and ArchR/scDblFinder (heterotypic) catch different classes.
- **Thresholds are not portable** - TSS/FRiP depend on annotation, pipeline, and counting convention; threshold from the data's own distributions.
- **Unified peaks across datasets** - re-quantify against one peak set; peak boundaries are not portable.
- **One genome build everywhere** - peaks, fragments, EnsDb, and BSgenome must match; a build mismatch runs without error but returns coordinate-mismatched garbage.

## Related Skills

single-cell/multimodal-integration - joining the ATAC modality with RNA (Multiome WNN/MultiVI)
single-cell/preprocessing - shared QC and filtering concepts from scRNA-seq
single-cell/clustering - clustering and UMAP shared with scRNA-seq
single-cell/doublet-detection - doublet concepts and rate expectations
atac-seq/atac-peak-calling - bulk ATAC peak-calling background (MACS shift/extend)
atac-seq/motif-deviation - chromVAR deviation scoring in depth
chip-seq/motif-analysis - motif databases (JASPAR/cisBP) and enrichment testing
