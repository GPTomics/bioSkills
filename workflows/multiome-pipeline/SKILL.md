---
name: bio-workflows-multiome-pipeline
description: Orchestrates the end-to-end 10x Multiome (paired scRNA + scATAC) pipeline from Cell Ranger ARC output to a jointly-embedded, annotated object, chaining per-modality QC, AMULET fragment-based ATAC doublet detection, per-modality normalization (RNA SCT/PCA; ATAC TF-IDF/LSI), WNN (or MultiVI) integration, joint clustering, RNA-based annotation, and LinkPeaks peak-to-gene linking. Use when enforcing the shared cell-barcode intersection between modalities (cellranger-ARC not -atac), keeping per-modality QC/doublets before the joint embedding, dropping the depth-correlated LSI component, annotating identity from RNA (ATAC is regulatory state), treating peak-to-gene links as correlational hypotheses, or aggregating to pseudobulk for cross-condition DE. Hands mechanism to the single-cell and atac-seq component skills; not a re-teach of any single step.
tool_type: r
primary_tool: Seurat
goal_approach_exempt: true
workflow: true
depends_on:
  - single-cell/data-io
  - single-cell/preprocessing
  - single-cell/clustering
  - single-cell/multimodal-integration
  - single-cell/scatac-analysis
  - atac-seq/single-cell-atac
  - atac-seq/co-accessibility
  - atac-seq/motif-deviation
qc_checkpoints:
  - after_loading: "Both modalities detected per cell"
  - after_rna_qc: "RNA quality filters passed"
  - after_atac_qc: "TSS enrichment >2, nucleosome signal <4"
  - after_wnn: "Joint embedding separates cell types"
---

## Version Compatibility

Reference examples tested with: Cell Ranger ARC 2.2+, Seurat 5.1+, Signac 1.14+, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38 1.4+, ggplot2 3.5+ (AMULET via the standalone java/python tool or scDblFinder's amulet() in R -- ArchR and snapATAC2 ship their OWN simulation-based doublet callers, not AMULET; MultiVI via scvi-tools if using the Python path)

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: `cellranger-arc count` (NOT `cellranger-atac`) emits the paired RNA+ATAC per-nucleus barcodes. Seurat `FindClusters(algorithm=3)` is SLM, not Leiden (1=Louvain, 2=Louvain-multilevel, 3=SLM, 4=Leiden). The `atac_fragments.tsv.gz` must be block-gzipped + tabix-indexed; the Tn5 +4/-5 offset is already applied by 10x -- do not re-shift. Confirm in-tool before quoting.

# Multiome Pipeline

**"Analyze my 10X Multiome data jointly"** -> Orchestrate Cell Ranger ARC processing, Seurat/Signac scRNA+scATAC integration via WNN, chromatin accessibility peak calling, motif enrichment, and gene regulatory network inference.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| `cellranger-arc` run (NOT `cellranger-atac`) | Only ARC emits the joined RNA+ATAC per-nucleus barcodes; -atac gives ATAC-only barcodes and the join is impossible |
| Shared cell-barcode join | RNA and ATAC QC pass DIFFERENT barcodes; the analyzable set is their INTERSECTION. A namespace mismatch (`-1` suffix, RNA vs ATAC whitelist) silently empties the join |
| Same genome build for GEX + ATAC (EnsDb + BSgenome) | Gene activity, LinkPeaks, and motif coordinates require identical build, else peak-to-gene linking is garbage |
| Consensus peak set (multi-sample) | Peaks are dataset-specific; merging samples on discordant peaks fabricates batch structure -- re-quantify against a unified peak set |
| Intron inclusion (GEX half) | Multiome is nuclei (mostly unspliced) -- introns are essential for the RNA UMI totals |

## Pipeline orchestration: the joint-modality decisions that make or break the result

Multiome's defining feature is that RNA and ATAC are measured in the SAME nucleus, so the two assays share one barcode universe and must be reconciled, not analyzed independently. The orchestration decisions:
- The paired-cell anchor is the whole point: keep only barcodes that pass QC in BOTH modalities. RNA and ATAC QC use different metrics (RNA: gene count, mito %; ATAC: TSS enrichment, nucleosome signal, fragment count) and are computed per modality, but the surviving cell set is their intersection. cellranger-arc (not cellranger-atac) produces the paired barcodes; their universes differ. See single-cell/preprocessing and single-cell/scatac-analysis.
- Per-modality QC and doublet detection run BEFORE the joint embedding. ATAC doublets are missed by RNA-based callers (scDblFinder/Scrublet, see single-cell/doublet-detection) and need a fragment-based detector such as AMULET (see single-cell/scatac-analysis); resolving them after WNN lets fake intermediate states drive the joint clustering.
- Drop the depth-correlated LSI component before joint analysis. The first ATAC LSI/SVD component usually (not always) captures sequencing depth rather than biology; check with DepthCor and exclude whichever component correlates with depth (commonly #1, hence dims 2:30 in WNN). Forgetting this lets depth dominate the ATAC contribution to the joint graph.
- WNN vs a generative joint embedding is a real choice. Seurat/Signac WNN learns a per-cell modality weight on precomputed PCA + LSI and is the default when both modalities are well processed; MultiVI (scvi-tools) jointly models RNA + ATAC counts end-to-end and handles batch and mosaic (RNA-only or ATAC-only) cells better. WNN cannot integrate cells missing a modality. See single-cell/multimodal-integration; verify current best practice against installed docs.
- Embed (WNN) then cluster then annotate, and annotate from RNA primarily. Gene-activity scores derived from ATAC are an approximation of expression, so cell-type labels come from RNA markers; ATAC informs the regulatory state, not the identity call. See single-cell/clustering, single-cell/markers-annotation, single-cell/cell-annotation.
- Peak-to-gene linking is correlational, not causal. LinkPeaks correlates peak accessibility with gene expression across cells within a window; a link is a hypothesis to validate, not a proven enhancer-target pair. For genome-wide enhancer-gene mapping use ABC/ENCODE-rE2G. See atac-seq/co-accessibility and atac-seq/enhancer-gene-linking.
- Cross-condition questions still need pseudobulk and separate composition testing. Condition DE on either modality aggregates RAW counts per sample x cell type (cells-as-replicates is pseudoreplication, Squair 2021); proportion shifts between conditions are tested separately and can masquerade as DE. See single-cell/differential-abundance and differential-expression/deseq2-basics.

## Workflow Overview

```
10X Multiome data
    |
    v
[1. Load Data] ---------> Read RNA + ATAC
    |
    v
[2. RNA Processing] ----> Standard scRNA workflow
    |
    v
[3. ATAC Processing] ---> Peak calling, LSI
    |
    v
[4. WNN Integration] ---> Weighted nearest neighbors
    |
    v
[5. Joint Analysis] ----> Clustering, markers
    |
    v
[6. Linked Features] ---> Gene-peak links
    |
    v
Integrated multiome object
```

## Step 1: Load Multiome Data

```r
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)

# Load RNA
rna_counts <- Read10X_h5('filtered_feature_bc_matrix.h5')
# For multiome, this returns a list with 'Gene Expression' and 'Peaks'

# Create Seurat object with RNA
seurat_obj <- CreateSeuratObject(
    counts = rna_counts$`Gene Expression`,
    assay = 'RNA'
)

# Load ATAC
atac_counts <- rna_counts$Peaks
# Or from fragments file
frags <- CreateFragmentObject('atac_fragments.tsv.gz', cells = colnames(seurat_obj))

# EnsDb returns Ensembl seqnames (1,2,X); cellranger-arc peaks/fragments are UCSC (chr1,chr2).
# Convert, or TSSEnrichment and LinkPeaks silently fail on zero seqname overlap.
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'

# Create ChromatinAssay
atac_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(':', '-'),
    fragments = frags,
    annotation = annotations
)

seurat_obj[['ATAC']] <- atac_assay
```

## Step 2: RNA Quality Control and Processing

```r
# QC metrics
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')

# Filter
seurat_obj <- subset(seurat_obj,
    nCount_RNA > 1000 &
    nCount_RNA < 25000 &
    percent.mt < 20
)

# Normalize RNA
seurat_obj <- SCTransform(seurat_obj, assay = 'RNA', verbose = FALSE)

# PCA
seurat_obj <- RunPCA(seurat_obj, assay = 'SCT', verbose = FALSE)
```

## Step 3: ATAC Quality Control and Processing

```r
# ATAC QC metrics
DefaultAssay(seurat_obj) <- 'ATAC'

seurat_obj <- NucleosomeSignal(seurat_obj)
seurat_obj <- TSSEnrichment(seurat_obj)

# Visualize
VlnPlot(seurat_obj, features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
        pt.size = 0, ncol = 3)

# Filter ATAC
seurat_obj <- subset(seurat_obj,
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    TSS.enrichment > 2 &
    nucleosome_signal < 4
)

# Normalize ATAC (TF-IDF + SVD = LSI)
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q0')
seurat_obj <- RunSVD(seurat_obj)

# Check LSI components (first often correlates with depth)
DepthCor(seurat_obj)
```

## Step 3b: Doublet detection (per modality, BEFORE WNN)

Remove doublets before the joint embedding, or fake intermediate states drive the joint clustering. RNA-based callers (scDblFinder/Scrublet) MISS ATAC doublets -- ATAC needs a fragment-based caller (AMULET), run on the same nuclei. Detect per modality, drop the union of doublets, then build WNN. Mechanism: single-cell/doublet-detection (RNA) and single-cell/scatac-analysis (AMULET).

## Step 4: Weighted Nearest Neighbors (WNN)

```r
# Build WNN graph using both modalities
seurat_obj <- FindMultiModalNeighbors(
    seurat_obj,
    reduction.list = list('pca', 'lsi'),
    dims.list = list(1:30, 2:30),  # Skip LSI component 1 if depth-correlated
    modality.weight.name = 'RNA.weight'
)

# UMAP on WNN graph
seurat_obj <- RunUMAP(seurat_obj, nn.name = 'weighted.nn',
                       reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_')

# Cluster on WNN
seurat_obj <- FindClusters(seurat_obj, graph.name = 'wsnn',
                            algorithm = 3, resolution = 0.5, verbose = FALSE)
```

## Step 5: Visualization and Markers

```r
# Compare modality-specific and joint embeddings
p1 <- DimPlot(seurat_obj, reduction = 'pca', label = TRUE) + ggtitle('RNA PCA')
p2 <- DimPlot(seurat_obj, reduction = 'lsi', label = TRUE) + ggtitle('ATAC LSI')
p3 <- DimPlot(seurat_obj, reduction = 'wnn.umap', label = TRUE) + ggtitle('WNN UMAP')
p1 + p2 + p3

# Modality weights per cell
VlnPlot(seurat_obj, features = 'RNA.weight', group.by = 'seurat_clusters', pt.size = 0)

# Find markers (RNA)
DefaultAssay(seurat_obj) <- 'SCT'
rna_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25)

# Find markers (ATAC - differentially accessible peaks)
DefaultAssay(seurat_obj) <- 'ATAC'
atac_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.05,
                                test.use = 'LR', latent.vars = 'nCount_ATAC')
```

## Step 6: Gene-Peak Linkage

```r
# Link peaks to genes
DefaultAssay(seurat_obj) <- 'ATAC'
seurat_obj <- RegionStats(seurat_obj, genome = BSgenome.Hsapiens.UCSC.hg38)

seurat_obj <- LinkPeaks(
    seurat_obj,
    peak.assay = 'ATAC',
    expression.assay = 'SCT',
    genes.use = c('CD8A', 'CD4', 'MS4A1', 'CD14')  # Example genes
)

# Visualize links
CoveragePlot(seurat_obj, region = 'CD8A', features = 'CD8A',
             expression.assay = 'SCT', extend.upstream = 10000, extend.downstream = 10000)
```

## Complete Workflow Script

```r
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

# Configuration
data_dir <- 'multiome_output'
output_dir <- 'multiome_results'
dir.create(output_dir, showWarnings = FALSE)

# === Load Data ===
cat('Loading data...\n')
counts <- Read10X_h5(file.path(data_dir, 'filtered_feature_bc_matrix.h5'))
frags <- file.path(data_dir, 'atac_fragments.tsv.gz')

seurat_obj <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = 'RNA')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'   # match cellranger-arc UCSC seqnames or TSS/LinkPeaks fail
seurat_obj[['ATAC']] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(':', '-'),
    fragments = frags,
    annotation = annotations
)
cat('Cells:', ncol(seurat_obj), '\n')

# === RNA QC ===
cat('RNA QC...\n')
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')
seurat_obj <- subset(seurat_obj, nCount_RNA > 1000 & nCount_RNA < 25000 & percent.mt < 20)

# === ATAC QC ===
cat('ATAC QC...\n')
DefaultAssay(seurat_obj) <- 'ATAC'
seurat_obj <- NucleosomeSignal(seurat_obj)
seurat_obj <- TSSEnrichment(seurat_obj)
seurat_obj <- subset(seurat_obj, nCount_ATAC > 1000 & TSS.enrichment > 2 & nucleosome_signal < 4)
cat('After QC:', ncol(seurat_obj), 'cells\n')

# === Process RNA ===
cat('Processing RNA...\n')
DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# === Process ATAC ===
cat('Processing ATAC...\n')
DefaultAssay(seurat_obj) <- 'ATAC'
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q0')
seurat_obj <- RunSVD(seurat_obj)

# === WNN Integration ===
cat('WNN integration...\n')
seurat_obj <- FindMultiModalNeighbors(seurat_obj,
    reduction.list = list('pca', 'lsi'),
    dims.list = list(1:30, 2:30),
    modality.weight.name = 'RNA.weight'
)
seurat_obj <- RunUMAP(seurat_obj, nn.name = 'weighted.nn',
    reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_')
seurat_obj <- FindClusters(seurat_obj, graph.name = 'wsnn', algorithm = 3, resolution = 0.5, verbose = FALSE)

# === Save ===
saveRDS(seurat_obj, file.path(output_dir, 'multiome_analyzed.rds'))

# === Plots ===
pdf(file.path(output_dir, 'wnn_umap.pdf'), width = 10, height = 8)
DimPlot(seurat_obj, reduction = 'wnn.umap', label = TRUE)
dev.off()

cat('Results saved to:', output_dir, '\n')
cat('Clusters:', length(unique(seurat_obj$seurat_clusters)), '\n')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Near-empty joint object / no cells after join | RNA vs ATAC barcode namespace mismatch (`-1` suffix, different whitelist) | Reconcile barcodes; intersect on identical strings; confirm cellranger-ARC (not -atac) |
| ATAC depth dominates the joint graph | Kept the depth-correlated LSI component | `DepthCor` -> drop it (WNN `dims.list` `2:30` for ATAC) |
| Fake intermediate joint clusters | ATAC doublets not removed (RNA caller is blind to them) | AMULET fragment-based doublet call per modality before WNN |
| Cell types mislabeled | Annotated from ATAC gene-activity | Annotate identity from RNA markers; activity is a cluster-level proxy |
| Spurious batch across multi-sample multiome | Merged on discordant peak sets | Unify peaks and re-quantify |
| "Enhancer regulates gene" overclaim | Read LinkPeaks correlation as causal | Treat as a composition-confounded hypothesis; validate |
| Inflated cross-condition DE | Tested cells as replicates on either modality | Pseudobulk RAW per sample x cell-type (Squair 2021) |

## References

- Hao Y, Hao S, Andersen-Nissen E, et al (2021) Integrated analysis of multimodal single-cell data. *Cell* 184:3573-3587.e29. DOI 10.1016/j.cell.2021.04.048. (WNN.)
- Ashuach T, Gabitto MI, Koodli RV, et al (2023) MultiVI: deep generative model for the integration of multimodal data. *Nature Methods* 20:1222-1231. DOI 10.1038/s41592-023-01909-9. (mosaic-capable joint RNA+ATAC alternative to WNN.)
- Squair JW, Gautier M, Kathe C, et al (2021) Confronting false discoveries in single-cell differential expression. *Nature Communications* 12:5692. DOI 10.1038/s41467-021-25960-2. (pseudobulk for cross-condition DE.)

## Related Skills

- single-cell/data-io - Loading 10X, h5ad, RDS, and h5mu formats
- single-cell/preprocessing - Per-modality QC and normalization choice
- single-cell/doublet-detection - RNA-based and hashing doublet removal
- single-cell/clustering - Resolution sweep and cluster validation on the joint graph
- single-cell/markers-annotation - Marker discovery, manual labeling, and pseudobulk condition DE
- single-cell/cell-annotation - Automated reference-based label transfer from the RNA modality
- single-cell/differential-abundance - Test whether cell-type proportions shifted between conditions
- single-cell/multimodal-integration - WNN, totalVI/MultiVI, and MOFA joint-embedding details
- single-cell/scatac-analysis - ATAC-specific processing, LSI, gene activity, and AMULET fragment-based doublet detection
- differential-expression/deseq2-basics - Pseudobulk condition DE engine for aggregated counts
- atac-seq/single-cell-atac - Signac / ArchR / SnapATAC2 ecosystem decision; AMULET; cellranger-arc
- atac-seq/co-accessibility - Cicero / ArchR getCoAccessibility for cis-regulatory inference
- atac-seq/enhancer-gene-linking - ABC / ENCODE-rE2G for enhancer-gene mapping
- atac-seq/motif-deviation - chromVAR for per-cell TF motif activity
- atac-seq/footprinting - scprinter for sc footprinting
- workflows/grn-pipeline - Downstream: the paired object feeds SCENIC+ enhancer-GRN inference (Path B)
