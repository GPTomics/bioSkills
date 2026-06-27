---
name: bio-single-cell-scatac-analysis
description: Analyze single-cell ATAC-seq with Signac/ArchR (R) and SnapATAC2 (Python alternative). Use when processing scATAC fragments, choosing a framework, calling consensus peaks, running TF-IDF/LSI while diagnosing the depth component, scoring chromVAR motif deviations against GC-matched backgrounds, detecting homotypic vs heterotypic doublets, or deciding whether to binarize the count matrix.
tool_type: r
primary_tool: Signac
---

## Version Compatibility

Reference examples tested with: Signac 1.13+, Seurat 5.0+, ArchR 1.0+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python (SnapATAC2 alternative): `pip show snapatac2` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# scATAC-seq Analysis

**"Analyze my single-cell ATAC-seq data"** -> Process fragments, QC on chromatin signal, reduce dimensions with TF-IDF/LSI, cluster, call consensus peaks per cell type, and score TF motif activity.
- R: `Signac::CreateChromatinAssay()` -> `RunTFIDF()` -> `FindTopFeatures()` -> `RunSVD()` -> `RunChromVAR()`
- R (large data, on-disk): `ArchR::createArrowFiles()` -> `addIterativeLSI()` -> `addReproduciblePeakSet()`
- Python (scverse, >1M cells): `snapatac2.pp.add_tile_matrix()` -> `tl.spectral()` -> `tl.macs3()`

## Governing Principle

A zero in the cell-by-peak matrix is epistemically ambiguous: it can mean "closed in this cell" (biology) or "accessible but no Tn5 fragment captured here" (sampling). With ~2 DNA copies per diploid locus and shallow per-cell coverage, sampling dominates the zeros. The matrix is near-binary by sampling statistics, not by biology; underlying accessibility is continuous but observed as a Bernoulli-like draw.

Binarization is now disfavored. Among non-zero entries the count (1 vs 2 vs >2) is informative, and collapsing to 1 discards it (Martens 2024). Model fragment counts with a count likelihood (Paired-Insertion Counting, SnapATAC2; PoissonVI), never read counts (PCR noise). Caveat: the extra information lives in the count=2 tier, so the benefit scales with sequencing depth; binarized analyses of shallow data are leaving little on the table, deep data more.

Per-cell signal is near-binary by sampling, so single-cell single-gene quantitative claims are unreliable; aggregate to cluster/pseudobulk for graded signal.

Gene-activity scores are a weak cluster-level proxy, structurally not just empirically: (1) enhancer-to-promoter assignment is unknown, and any fixed-distance heuristic (Signac gene body + 2 kb, ArchR exponential decay to 100 kb) is wrong for genes whose enhancers sit outside the window or loop differently by cell type; (2) poised/bivalent promoters are accessible while the gene is silent, so accessibility-to-expression is not monotone. Use gene activity for cluster-level annotation and scRNA-integration anchoring only, never as a single-cell transcriptome surrogate.

The peak set depends on which cells called the peaks: peaks are called on cells already grouped, but the grouping used a feature matrix that depends on a peak/tile choice. This is a circularity. Rare populations unresolved in the first pass never get their peaks called, so their defining elements stay invisible, a self-reinforcing blind spot. This is why iterative per-cluster peak calling (ArchR iterative LSI) exists, and why testing differential accessibility on a peak set called from the same clustering is double-dipping.

## Framework Decision Table

Framework choice is an infrastructure decision (language, memory, multimodal needs), not a statistics decision. Scalability numbers describe each tool's most-optimized path, not every operation.

| Framework | Language / storage | Use when | Fails when |
|---|---|---|---|
| Signac | R, in-memory Seurat ChromatinAssay | Seurat-integrated multimodal (WNN), familiar Seurat API | ~10^5+ cells (RAM-bound; `future` workers copy the object) |
| ArchR | R, on-disk HDF5 Arrow files | Large R workflows (~1M cells), built-in iterative LSI/peak/GRN suite | Networked filesystems (HDF5 file-locking); not a portable matrix |
| SnapATAC2 | Python+Rust, AnnData backed | >1M cells (matrix-free spectral), scverse/scvi-tools stack | Less turnkey footprinting/GRN; faster-moving 2.x API |
| muon + scanpy | Python, MuData | Multimodal Python container (RNA+ATAC) | Not ATAC-optimized for the heaviest steps |

R<->Python interop (reticulate, zellkonverter, sceasy) loses information (ChromatinAssay slots, ArchR HDF5 do not round-trip); plan to stay in one ecosystem. Verify the current best-practice default against installed docs before committing.

## Matrix Type: Tile vs Peak vs Gene Activity

| Matrix | What it is | Use when | Caveat |
|---|---|---|---|
| Tile/bin (500 bp) | Genome binned, no prior peaks | Initial LSI/clustering before peaks exist | Not biology-aware; 500 bp tiles vs 501 bp peaks (off-by-one feature bugs) |
| Peak (consensus) | Per-cluster MACS peaks merged to fixed width | Final accessibility quantification, DA testing | Requires peaks first; circular with clustering |
| Gene activity | Accessibility folded to per-gene scalar | Cluster annotation, scRNA-integration anchors | Weak proxy; repressed/bivalent genes fail; distal enhancers misassigned |

## TF-IDF + LSI: Diagnose the Depth Component

**Goal:** Reduce the sparse, near-binary, depth-confounded matrix without letting technical depth dominate.

**Approach:** Reweight peaks with TF-IDF, reduce with truncated SVD, then drop components that correlate with depth, diagnosed by `DepthCor`, not blindly dropping component 1.

```r
obj <- RunTFIDF(obj)                       # method 1 (default) = log(TF x IDF), Stuart & Butler
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)                         # writes the 'lsi' reduction

DepthCor(obj, n = 10)                      # per-component Pearson correlation with nCount
# LSI_1 usually has |corr| > 0.95 with depth, but verify; occasionally it is component 2/3, or none
```

Component 1 captures depth ~90% of the time but the rule is symptom-based: compute each component's depth correlation (`DepthCor`, or ArchR `corCutOff = 0.75`) and drop whichever exceed the threshold. ArchR `addIterativeLSI()` recomputes LSI on variable features across clustering passes to reduce depth/batch artifacts. A reviewer flags blind `dims = 2:30` with no depth-correlation diagnostic.

## Clustering on LSI

**Goal:** Cluster cells from the depth-cleaned LSI embedding.

**Approach:** Build the neighbor graph and UMAP on the retained LSI dimensions, then cluster.

```r
dims_use <- 2:30                            # set from DepthCor, not assumed
obj <- RunUMAP(obj, reduction = 'lsi', dims = dims_use)
obj <- FindNeighbors(obj, reduction = 'lsi', dims = dims_use)
obj <- FindClusters(obj, algorithm = 3, resolution = 0.5)   # algorithm 3 = SLM
```

## Consensus Peak Calling

**Goal:** Call peaks per cell type and merge into a non-overlapping, reusable feature set, avoiding bias toward abundant cell types.

**Approach:** Pooled bulk calling misses rare-population elements; call per cluster on pseudobulk, then merge. ArchR's fixed-width iterative-overlap set is the most reproducible; Signac's `CallPeaks` is simpler but uses a variable-width union that drops significance metadata.

```r
peaks <- CallPeaks(obj, group.by = 'seurat_clusters')       # per-group MACS, then GRanges::reduce() union
peak_counts <- FeatureMatrix(fragments = Fragments(obj), features = peaks, cells = colnames(obj))
obj[['peaks']] <- CreateChromatinAssay(counts = peak_counts, fragments = Fragments(obj), annotation = Annotation(obj))
```

Fixed-width peaks (ArchR's 501 bp) remove per-peak length normalization and give a stable reusable feature space. ArchR ranks fixed-width candidates by significance, keeps the best, removes overlappers, and requires a peak in >=2 pseudobulk replicates (reproducibility, orthogonal to MACS q-value). Wrapper parameters differ (ArchR shift -75/extsize 150 with `--nolambda`; Signac/SnapATAC2 shift -100/extsize 200), which changes which weak peaks survive. Comparing peak sets across datasets requires re-quantifying against a unified set; peak boundaries are not portable.

## Differential Accessibility

**Goal:** Find peaks more accessible in one group, controlling for the depth confounder.

**Approach:** Use a logistic-regression test with total fragments as a latent variable; do not test on a peak set called from the same clustering being compared (double-dipping).

```r
DefaultAssay(obj) <- 'peaks'
da <- FindMarkers(obj, ident.1 = 'cluster1', ident.2 = 'cluster2',
                  test.use = 'LR', latent.vars = 'nCount_peaks')
```

## chromVAR Motif Deviations

**Goal:** Find which TF motifs vary in accessibility across cells, corrected for GC content and depth.

**Approach:** Attach motif matches, then compute deviations against a GC- and accessibility-matched background; rank with z-scores, never raw deviations.

```r
library(JASPAR2020); library(TFBSTools); library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

pfm <- getMatrixSet(JASPAR2020, opts = list(collection = 'CORE', tax_group = 'vertebrates', all_versions = FALSE))
obj <- AddMotifs(obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
obj <- RunChromVAR(obj, genome = BSgenome.Hsapiens.UCSC.hg38)   # GC-matched background internally

DefaultAssay(obj) <- 'chromvar'
diff_motifs <- FindMarkers(obj, ident.1 = 'cluster1', ident.2 = 'cluster2',
                           mean.fxn = rowMeans, fc.name = 'avg_diff')
```

chromVAR's deviation is meaningful only against a GC- and accessibility-matched background; an unmatched background manufactures apparent enrichment for GC-rich motifs (most TF motifs are GC-rich). Use z-scores (background-normalized) for cross-motif ranking, raw deviations are not comparable across motifs. Motif != TF: paralogous TFs share near-identical motifs, so an enriched motif implicates a family, not a factor; motif presence != occupancy; and a footprint (TOBIAS, needs pseudobulk) is stronger occupancy evidence than motif-in-peak. Disambiguate with TF expression (Multiome) before claiming "TF X drives this program".

## Gene Activity (Cluster-Level Only)

**Goal:** Approximate per-gene accessibility for marker-based annotation and scRNA anchoring.

**Approach:** Sum fragments over the gene body plus a promoter window; treat the output as a cluster-level aid, not measured RNA.

```r
gene_act <- GeneActivity(obj)              # gene body + 2 kb upstream, flat count, no distance weighting
obj[['ACT']] <- CreateAssayObject(counts = gene_act)
obj <- NormalizeData(obj, assay = 'ACT', scale.factor = median(obj$nCount_ACT))
```

## Doublet Detection: Homotypic vs Heterotypic

Two strategies catch different doublet classes; run both and combine. Doublet callers are separate from QC metrics (TSS/nucleosome gate debris, not doublets).

| Tool | Principle | Catches | Key dependency |
|---|---|---|---|
| AMULET | >2 fragments overlapping a diploid locus -> Poisson + BH | Homotypic (same-type) | ~25k read pairs/cell for full recall |
| ArchR `addDoubletScores` | Simulate doublets -> LSI/UMAP -> kNN; use `DoubletEnrichment` | Heterotypic (different-type) | LSI/UMAP quality; structurally blind to homotypic |
| scDblFinder ATAC | Simulate on `nfeatures=25` aggregated meta-features | Heterotypic | Embedding quality |

AMULET silently under-calls below ~25k coverage; CNV/aneuploidy breaks its diploid null (amplified loci exceed 2 copies in true singlet cancer cells -> false positives); multinucleate/S-G2-M cells violate the <=2-copies assumption. ArchR prefers `DoubletEnrichment` over `DoubletScore`. scDblFinder uses `nfeatures=25` (not 1000) and its authors recommend against `clamulet`.

## QC Thresholds

| Metric | Signac column | Threshold | Basis |
|---|---|---|---|
| TSS enrichment | `TSS.enrichment` | >2-3 (Signac); >4 (ArchR human) | ENCODE signal/noise; threshold is annotation-dependent, not portable |
| Total fragments | `nCount_peaks` / nFrags | >1000 (often >3000) | removes empties/debris |
| Nucleosome signal | `nucleosome_signal` | <4 | banding quality; very low can mean over-transposition |
| FRiP | `FRiP` | >0.15-0.40 (study-dependent) | signal in peaks; depends on peak set and counting convention |

TSS scores are not comparable across pipelines/annotations; never port thresholds. `TSSEnrichment(fast=TRUE)` blocks later `TSSPlot()`. FRiP depends on the peak set (circular if the same cells) and counting convention (Signac counts fragments, CellRanger-ATAC counts insertions). Threshold from the joint distributions of the actual data, not copied defaults.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| UMAP separates by depth, not biology | Did not drop the depth-correlated LSI component | Run `DepthCor`; drop components above threshold (often #1, verify) |
| Long flat run of zeros read as "closed" | Zeros are sampling-dominated, ambiguous | Interpret at cluster/pseudobulk level; check effective coverage before structural claims |
| Gene activity disagrees with RNA for a marker | Repressed/bivalent promoter is open but silent; distal enhancer outside window | Use gene activity for cluster annotation only; validate with multiome RNA |
| "Everything is GC-rich enriched" in chromVAR | Unmatched background | Use `getBackgroundPeaks`/RunChromVAR GC+accessibility-matched background; report z-scores |
| Rare cell type never appears | Peaks called from a coarse single-pass clustering missed its elements | Iterative per-cluster peak calling + re-clustering (ArchR iterative LSI) |
| DA peaks look inflated | Tested on a peak set called from the same clustering (double-dipping) | Call peaks independently of the comparison; treat as ranking |
| Doublets pass QC | TSS/nucleosome gate debris, not doublets | Run AMULET (homotypic) and ArchR/scDblFinder (heterotypic) and combine |
| AMULET finds few doublets in cancer | CNV breaks the diploid null; or coverage <25k pairs/cell | Use heterotypic callers in aneuploid samples; check per-cell coverage |
| "TF X drives this" from a motif | Motif implicates a family; presence != occupancy | Confirm with TF expression (multiome) and/or footprint (TOBIAS, pseudobulk) |
| QC, gene activity, and motifs all run but look wrong | Peaks, fragments, EnsDb annotation, and BSgenome are on different genome builds; coordinate mismatch is silently wrong (no crash) | Pin every reference to one build (e.g. all hg38); verify TSS enrichment and a known marker before trusting downstream |

## Related Skills

- single-cell/multimodal-integration - joining the ATAC modality with RNA (Multiome WNN/MultiVI)
- single-cell/preprocessing - shared QC and filtering concepts from scRNA-seq
- single-cell/clustering - clustering and UMAP shared with scRNA-seq
- single-cell/doublet-detection - doublet concepts and rate expectations
- atac-seq/atac-peak-calling - bulk ATAC peak-calling background (MACS shift/extend)
- atac-seq/motif-deviation - chromVAR deviation scoring in depth
- chip-seq/motif-analysis - motif databases (JASPAR/cisBP) and enrichment testing

## References

Buenrostro JD, Giresi PG, Zaba LC, et al. Transposition of native chromatin for fast and sensitive epigenomic profiling (ATAC-seq). Nat Methods 10(12):1213-1218 (2013).
Cusanovich DA, Daza R, Adey A, et al. Multiplex single-cell profiling of chromatin accessibility (TF-IDF/LSI). Science 348(6237):910-914 (2015).
Stuart T, Srivastava A, Madad S, Lareau CA, Satija R. Single-cell chromatin state analysis with Signac. Nat Methods 18:1333-1341 (2021).
Granja JM, Corces MR, Pierce SE, et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nat Genet 53:403-411 (2021).
Zhang K, Zemke NR, Armand EJ, Ren B. A fast, scalable and versatile tool for analysis of single-cell omics data (SnapATAC2). Nat Methods 21(2):217-227 (2024).
Schep AN, Wu B, Buenrostro JD, Greenleaf WJ. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. Nat Methods 14(10):975-978 (2017).
Martens LD, Fischer DS, Theis FJ, Buettner F. Modeling fragment counts improves single-cell ATAC-seq analysis. Nat Methods 21(1):28-31 (2024).
Miao Z, Kim J. Uniform quantification of single-nucleus ATAC-seq data with Paired-Insertion Counting (PIC) and a model-based insertion rate estimator. Nat Methods 21:32-36 (2024).
Thibodeau A, Eroglu A, McGinnis CS, et al. AMULET: a novel read count-based method for effective multiplet detection from single-nucleus ATAC-seq data. Genome Biol 22:252 (2021).
Germain P-L, Lun A, Garcia Meixide C, Macnair W, Robinson MD. Doublet identification in single-cell sequencing data using scDblFinder. F1000Research 10:979 (2022).
Bentsen M, Goymann P, Schultheis H, et al. ATAC-seq footprinting unravels kinetics of transcription factor binding during zygotic genome activation (TOBIAS). Nat Commun 11:4267 (2020).
