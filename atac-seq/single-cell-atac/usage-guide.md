# Single-Cell ATAC-seq - Usage Guide

## Overview

End-to-end single-cell ATAC-seq processing and analysis. Covers ecosystem selection (Signac, ArchR, SnapATAC2), per-cell QC thresholds (looser than bulk because of per-cell sparsity), doublet detection (AMULET, ArchR, scDblFinder), TF-IDF + LSI dimensionality reduction (skipping the depth-correlated first component), per-cluster pseudobulk peak calling, and 10X Multiome RNA+ATAC integration via WNN. Includes per-tool failure modes and reconciliation between ecosystems.

## Prerequisites

```r
# Signac ecosystem
BiocManager::install(c('GenomicRanges', 'EnsDb.Hsapiens.v86', 'BSgenome.Hsapiens.UCSC.hg38'))
install.packages(c('Signac', 'Seurat'))

# ArchR ecosystem
remotes::install_github('GreenleafLab/ArchR', ref='master', repos=BiocManager::repositories())
ArchR::installExtraPackages()
```

```bash
# SnapATAC2 (Python)
pip install snapatac2

# Doublet detection
pip install amulet-py
```

```r
BiocManager::install('scDblFinder')
```

Inputs (from 10X Cell Ranger ATAC or Multiome):
- `outs/fragments.tsv.gz` (and `.tbi` index)
- `outs/filtered_peak_bc_matrix.h5` (Signac path)
- `outs/singlecell.csv` (per-barcode metadata)

## Quick Start

Tell your AI agent what you want to do:
- "Process 10X scATAC with Signac: TF-IDF + LSI (dims 2:30, skip depth) + UMAP + Leiden"
- "Use ArchR for a 200K-cell dataset because Signac will run out of memory"
- "Run AMULET for ATAC-specific doublet detection; cross-validate with ArchR's doublet score"
- "Build per-cluster pseudobulk peaks with MACS3, then iterative-overlap consensus"
- "Integrate 10X Multiome via WNN with PCA on RNA and LSI(2:30) on ATAC"
- "Re-filter cellranger output at fragment count >= 1000 AND TSS enrichment >= 4 because cellranger's cell calling is too lenient"

## Example Prompts

### Standard Signac Pipeline
> "Load 10X scATAC h5 + fragments file into Signac. Compute NucleosomeSignal and TSSEnrichment. Filter cells: peak fragments 1000-20000, FRiP >= 15%, blacklist ratio < 0.05, nucleosome signal < 4, TSS enrichment > 4. Then RunTFIDF -> RunSVD -> RunUMAP(dims=2:30) -> FindNeighbors(dims=2:30) -> FindClusters(algorithm=4, resolution=0.5)."

### Large Dataset with ArchR
> "150K-cell scATAC dataset. Use ArchR: createArrowFiles with minTSS=4 and minFrags=1000, addDoubletScores, filterDoublets, addIterativeLSI on TileMatrix, addClusters Leiden, addUMAP. Then addReproduciblePeakSet per cluster."

### SnapATAC2 Python Workflow
> "Use SnapATAC2: import fragments.tsv.gz via `snap.pp.import_fragments` with hg38 chrom_sizes, compute tsse, filter min_counts=1000 and min_tsse=4, add_tile_matrix bin_size=500, select_features 250k, spectral, leiden, umap."

### Doublet Reconciliation
> "Run AMULET (collision-based) and ArchR's doublet score in parallel; flag cells called doublet by both as high-confidence; report flagging rate per cluster."

### Per-Cluster Pseudobulk DA
> "Aggregate cells per cluster into pseudobulk BAMs, run MACS3 callpeak per cluster, build iterative-overlap consensus across clusters, then DESeq2 differential between cluster A and B (atac-seq/differential-accessibility)."

### Multiome WNN
> "I have Multiome data. Run RNA workflow (LogNormalize -> PCA -> dims 1:30) and ATAC workflow (TF-IDF -> SVD -> dims 2:30) separately, then FindMultiModalNeighbors and RunUMAP nn.name='weighted.nn' for the joint embedding."

### Cell Type Annotation
> "Compute gene activity scores with Signac::GeneActivity (or ArchR getGeneScore); transfer cell type labels from a reference scRNA-seq atlas via Seurat::FindTransferAnchors."

## What the Agent Will Do

1. Verify input is Cell Ranger ATAC or Multiome output (chemistry version matters)
2. Re-filter cellranger cell calls at fragment >= 1000 AND TSS enrichment >= 4
3. Compute per-cell QC: NucleosomeSignal, TSS enrichment, FRiP, mt fraction, blacklist ratio
4. Run doublet detection (AMULET primary; ArchR or scDblFinder confirmatory)
5. Choose ecosystem (Signac for Multiome / Seurat user; ArchR for large; SnapATAC2 for Python)
6. TF-IDF + LSI/spectral; SKIP depth-correlated first component (dims 2:30)
7. UMAP + Leiden clustering
8. Per-cluster pseudobulk peak calling (MACS3 with ATAC parameters)
9. Iterative-overlap consensus across clusters
10. Gene activity scores for cell type annotation
11. Integration with paired scRNA-seq if Multiome (WNN)
12. Optional: trajectory analysis (ArchR getTrajectory or Cicero)

## Per-Cell QC Quick Reference

| Metric | Pass | Reject |
|--------|------|--------|
| Fragment count | 3000-50000 | < 1000 or > 80000 |
| TSS enrichment | >= 4 | < 2 |
| FRiP | >= 0.15 | < 0.05 |
| Mt fraction | < 0.05 | > 0.20 |
| Nucleosome signal | < 4 | > 10 |
| Doublet score | < 0.5 (AMULET) | > 0.7 |

Per-cell thresholds are looser than bulk; cells are sparse by design.

## Ecosystem Quick Reference

| Use case | Ecosystem |
|----------|----------|
| 10X Multiome (RNA + ATAC) | Signac (Seurat integration) |
| Large dataset (>100K cells) | ArchR (Arrow files; memory-efficient) |
| Python ecosystem | SnapATAC2 |
| Trajectory analysis | ArchR getTrajectory; Cicero (Signac) |
| ATAC-only standard analysis | Either Signac or ArchR |

## Doublet Detection Quick Reference

| Tool | Method | When to use |
|------|--------|-------------|
| AMULET | Collision-based | Primary for ATAC; orthogonal to clustering |
| ArchR addDoubletScores | Synthetic + LSI projection | Built into ArchR; auto-thresholds |
| scDblFinder | Synthetic + classifier | RNA-developed; works on Signac/SCE |

Run two; intersection of flagged cells is high-confidence doublet set.

## Tips

- The single most common bug: including LSI component 1 in clustering. It correlates with depth, not biology. Always use `dims=2:30`.
- Cell Ranger ATAC's cell calling is lenient. Always re-filter at fragment >= 1000 AND TSS enrichment >= 4.
- Per-cell thresholds are looser than bulk; the population aggregate is what matters statistically.
- AMULET is ATAC-specific (collision detection from 2-allele biology); pair with ArchR or scDblFinder for orthogonal evidence.
- Multiome WNN equal-weights modalities by default; ATAC's sparseness can dominate noise. Inspect weights and adjust if needed.
- Per-cluster pseudobulk peak calling needs >= 200 cells per cluster; aggregate small clusters or use union peaks.
- Gene activity scores (Signac::GeneActivity / ArchR getGeneScore) are approximate but useful for cross-modality cell-type transfer.
- chromVAR / AddMotifs MUST run after peakset finalized; otherwise motif annotations become stale.
- For trajectory analysis, ArchR's getTrajectory is more polished than Cicero; Cicero is better for cis-regulatory inference.
- Plant / non-model organisms need custom EnsDb / TxDb / BSgenome objects; build with `txdbmaker` or `AnnotationDbi`.
- 10X Multiome ATAC peakset MUST be called from the Multiome ATAC fragments, not transferred from a separate scATAC dataset.

## Related Skills

- atac-seq/atac-qc - Bulk QC adapted for per-cell
- atac-seq/atac-peak-calling - Pseudobulk peak calling
- atac-seq/consensus-peakset - Per-cluster consensus
- atac-seq/differential-accessibility - Pseudobulk DA
- atac-seq/motif-deviation - chromVAR for per-cell TF activity
- atac-seq/footprinting - scprinter for sc footprinting
- atac-seq/co-accessibility - Cicero / SCENIC+
- single-cell/preprocessing - General sc patterns
- single-cell/clustering - Cluster definition
- single-cell/cell-annotation - Marker-based annotation
- single-cell/multimodal-integration - Multiome integration
- single-cell/scatac-analysis - Cross-reference single-cell ATAC-specific patterns
