---
name: bio-atac-seq-single-cell-atac
description: Process and analyze single-cell ATAC-seq data with Signac, ArchR, SnapATAC2, or Cell Ranger ATAC. Use when handling 10X scATAC or 10X Multiome (paired RNA+ATAC) data, performing per-cell QC, choosing between ArchR/Signac/SnapATAC2 ecosystems, building per-cluster consensus peaksets, integrating with paired scRNA-seq, doublet detection (AMULET vs ArchR vs scDblFinder), or running pseudobulk differential accessibility per cluster.
tool_type: mixed
primary_tool: Signac
---

## Version Compatibility

Reference examples tested with: Cell Ranger ATAC 2.1+, Signac 1.13+, Seurat 5.0+, ArchR 1.0.2+, SnapATAC2 2.8+, AMULET 1.1+, scDblFinder 1.16+, scater 1.30+, scvi-tools 1.1+, GenomicRanges 1.54+, JASPAR2024 0.99+, BSgenome.Hsapiens.UCSC.hg38 1.4+, EnsDb.Hsapiens.v86 2.99+, MACS3 3.0+. SnapATAC2 2.8+ uses `pp.import_fragments`; older 2.5-2.7 used `pp.import_data` (renamed/removed in 2.9).

Verify before use:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# Single-Cell ATAC-seq

**"Process my 10X scATAC data from cellranger output"** -> Build a per-cell fragment matrix, compute per-cell QC, dimensionality reduction (TF-IDF + LSI / spectral / autoencoder), cluster, call cluster-level pseudobulk peaks, annotate cell types via gene-activity scores, and integrate with paired scRNA-seq if Multiome.

- R: `Signac::CreateChromatinAssay()` -> `Seurat` workflow (TF-IDF + SVD + UMAP + Leiden)
- R: `ArchR::createArrowFiles()` -> ArchR project (TileMatrix + LSI + UMAP)
- Python: `snapatac2.pp.import_fragments()` -> SnapATAC2 (spectral / diffusion-map clustering)
- CLI (preprocessing): `cellranger-atac count` (10X) or `chromap` (alignment-only fragment files)

## Ecosystem Choice (The Most Important Decision)

| Ecosystem | Language | Strength | Fails when | Best for |
|-----------|---------|---------|------------|----------|
| Signac (Stuart 2021) | R, Seurat-based | Tightest scRNA-seq integration; Seurat ecosystem mature | Memory hungry on >100K cells; slower than ArchR | Multiome RNA+ATAC; small-to-medium datasets; Seurat user |
| ArchR (Granja 2021) | R, Arrow/HDF5 | Memory-efficient (Arrow files); fast on 100K-1M cells; built-in trajectory + doublet | Less RNA-seq integration; ArchR-specific format | Large bulk-cohort scATAC; trajectory analysis; ATAC-only |
| SnapATAC2 (Zhang 2024) | Python, AnnData | Memory-efficient; modern Python ecosystem; spectral clustering performant | Newer; benchmarks evolving; ecosystem smaller than R | Python-first labs; very large datasets (>1M cells) |
| Cell Ranger ATAC | CLI (10X-specific) | 10X official preprocessing | Closed; fixed pipeline | Only as preprocessing step; analysis happens elsewhere |
| scATAC-pro | CLI-based pipeline | Alternative preprocessing | Less maintained | Legacy; not recommended for new projects |

Methodology evolves; verify against Granja 2021 (ArchR), Stuart 2021 (Signac), Zhang 2024 (SnapATAC2), Heumos 2023 (best practices) before locking pipelines.

## 10X Multiome Caveat (Paired RNA + ATAC)

10X Multiome chemistry profiles RNA AND ATAC from the same cell. Outputs are joined by a shared barcode. Multiome workflows use Signac for ATAC and Seurat for RNA, integrated through the same Seurat object via WNN (Weighted Nearest Neighbor; Hao 2021).

Single-modality 10X scATAC (chemistry v1, v2) does NOT produce paired RNA. Verify the chemistry on the cellranger summary before assuming Multiome.

## Per-Cell QC Thresholds

| Metric | Definition | Pass | Caution | Reject | Source |
|--------|-----------|------|---------|--------|--------|
| Fragment count per cell | n_fragments after dedup | 3000-50000 | 1000-3000 | < 1000 or > 80000 | 10X recommendation; high = doublet |
| TSS enrichment per cell | Signal at TSS / flanks (per cell) | >= 4 | 2-4 | < 2 | ArchR / Signac default; lower than bulk because per-cell |
| Nucleosome signal | Mono / NFR fragment ratio | <= 4 | 4-10 | > 10 | Signac default; high = poor library |
| % reads in peaks (per cell) | FRiP per cell at consensus | >= 0.15 | 0.10-0.15 | < 0.05 | ArchR / Signac defaults |
| Mitochondrial fraction (per cell) | chrM / total per cell | < 0.05 | 0.05-0.15 | > 0.20 | Lower than bulk; per-cell more sensitive |
| Doublet score | AMULET / ArchR doublet | < 0.5 | 0.5-0.7 | > 0.7 | Tool-dependent threshold |
| Blacklist ratio | Reads in ENCODE blacklist / total | < 0.05 | 0.05-0.10 | > 0.10 | Standard |

Per-cell thresholds are looser than bulk because individual cells have orders of magnitude less signal; the population aggregate is what matters.

## Doublet Detection: Three Approaches

| Tool | Method | Strength | Fails when |
|------|--------|---------|------------|
| AMULET (Thibodeau 2021) | Collision-based: detects cells with too many fragments at same position (mathematically impossible from single cell because of 2-allele limit) | Specific to ATAC biology; orthogonal to clustering | Requires high depth; sensitive at >5000 fragments/cell |
| ArchR addDoubletScores | Synthetic doublet simulation + projection into LSI | Built into ArchR; auto-thresholds | Tied to ArchR's LSI; not portable |
| scDblFinder (Germain 2021) | Synthetic doublets + classifier | Works on Signac and SCE objects; well-benchmarked | RNA-developed; ATAC adaptation requires careful settings |

**Operational rule:** Run AMULET as the primary ATAC-specific check; verify with ArchR or scDblFinder as orthogonal evidence. Double-flagged cells are high-confidence doublets.

## Per-Tool Failure Modes

### Signac TF-IDF + SVD -- First component is depth

**Trigger:** Running `RunTFIDF` -> `RunSVD` and using all components for clustering.

**Mechanism:** Component 1 of LSI is highly correlated with sequencing depth per cell, not biology. Including it pulls clusters along depth axis.

**Symptom:** UMAP shows linear cell-density gradient that tracks fragment count.

**Fix:** Skip component 1 in downstream UMAP and clustering: `RunUMAP(object, dims=2:30)`, `FindNeighbors(object, dims=2:30)`. ArchR and SnapATAC2 do this automatically.

### ArchR -- TileMatrix vs PeakMatrix confusion

**Trigger:** Using TileMatrix for differential testing or motif analysis.

**Mechanism:** TileMatrix is regular fixed bins (default 500 bp) covering the whole genome; useful for clustering but not for biology because ATAC signal is at peaks, not arbitrary bins.

**Fix:** Use TileMatrix for embedding/clustering (it's faster); use PeakMatrix (after `addReproduciblePeakSet`) for differential, motif, and gene-activity analysis.

### SnapATAC2 -- Memory layout assumes integer counts

**Trigger:** Loading non-integer or negative-valued matrices.

**Fix:** SnapATAC2 expects raw fragment counts (Int32). Convert from float matrices before loading.

### Cell Ranger ATAC -- Empty droplet detection

**Trigger:** Default cellranger-atac calls cells based on UMI count but ATAC has no UMIs; uses fragment-based heuristic. Sometimes calls < 1000-fragment cells as real.

**Fix:** Re-filter cellranger output to require fragment count >= 1000 AND TSS enrichment >= 4 BEFORE downstream analysis.

### Multiome WNN -- ATAC weighting

**Trigger:** Default WNN equally weights RNA and ATAC modalities.

**Mechanism:** ATAC's per-cell signal is much sparser than RNA; equal weighting can swamp the joint embedding with ATAC noise.

**Fix:** Inspect per-modality weights with `IntegrateLayers`; if ATAC weights dominate noise, manually adjust or use `FindMultiModalNeighbors` carefully.

### Per-cluster pseudobulk peak calling -- Empty clusters

**Trigger:** Calling MACS3 per cluster when one cluster has < 200 cells.

**Mechanism:** MACS3 needs >= 1M reads in pseudobulk to call peaks reliably; small clusters do not produce enough reads.

**Fix:** Aggregate small clusters into a "rare" group OR drop them from per-cluster calling. Use the union of larger-cluster peaks as a fallback for rare-cell-type analysis.

## Decision Tree by Goal

| Goal | Recommended pipeline |
|------|---------------------|
| Standard 10X scATAC analysis (R user) | Signac: CreateChromatinAssay -> RunTFIDF -> RunSVD -> RunUMAP (dims 2:30) -> FindClusters |
| Standard 10X scATAC analysis (Python user) | SnapATAC2: pp.import_fragments -> add_tile_matrix -> spectral -> UMAP -> leiden |
| Large dataset (> 100K cells) | ArchR (memory-efficient Arrow files) |
| 10X Multiome (paired RNA + ATAC) | Signac + Seurat: per-modality embedding then WNN integration |
| Trajectory / pseudotime | ArchR getTrajectory; or Signac + Cicero |
| Differential accessibility per cluster | Pseudobulk per cluster -> consensus peakset -> DESeq2 |
| Cell-type annotation | Gene-activity scores via ArchR or Signac; then run `single-cell/markers-annotation` |
| Multimodal trajectories (RNA + ATAC) | MOFA+, ArchR + scRNA integration, or SCENIC+ |
| Plant / non-model | Signac with custom EnsDb / TxDb; ArchR with custom annotations |

## Standard Signac Workflow

```r
library(Signac); library(Seurat); library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# 1. Load 10X output
counts <- Read10X_h5('outs/filtered_peak_bc_matrix.h5')
metadata <- read.csv('outs/singlecell.csv', row.names=1)

chrom_assay <- CreateChromatinAssay(
    counts=counts,
    sep=c(':', '-'),
    fragments='outs/fragments.tsv.gz',
    annotation=GetGRangesFromEnsDb(EnsDb.Hsapiens.v86),
    genome='hg38')

obj <- CreateSeuratObject(counts=chrom_assay, assay='ATAC', meta.data=metadata)

# 2. Per-cell QC
obj <- NucleosomeSignal(obj)              # Mono/NFR ratio
obj <- TSSEnrichment(obj, fast=FALSE)
obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments

# 3. QC filter (looser per-cell thresholds)
obj <- subset(obj,
    subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 &
             pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
             nucleosome_signal < 4 & TSS.enrichment > 4)

# 4. Dimensionality reduction (skip component 1 - it's depth)
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff='q0')
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction='lsi', dims=2:30)
obj <- FindNeighbors(obj, reduction='lsi', dims=2:30)
obj <- FindClusters(obj, algorithm=4, resolution=0.5)    # Leiden

# 5. Gene-activity for annotation
gene.activities <- GeneActivity(obj)
obj[['ACT']] <- CreateAssayObject(counts=gene.activities)
DefaultAssay(obj) <- 'ACT'
obj <- NormalizeData(obj, normalization.method='LogNormalize',
                     scale.factor=median(obj$nCount_ACT))
```

## Per-Cluster Pseudobulk Peak Calling

```r
# Signac wrapper around MACS3
peaks <- CallPeaks(obj, group.by='seurat_clusters',
                   macs2.path='/path/to/macs3', cleanup=FALSE,
                   format='BAMPE', extra.args='--shift -75 --extsize 150 -p 0.01')
# Then iterative-overlap consensus across clusters (atac-seq/consensus-peakset)
```

## ArchR Workflow

```r
library(ArchR)
addArchRGenome('hg38')

# 1. Create Arrow files from fragment files
ArrowFiles <- createArrowFiles(
    inputFiles=c('rep1=fragments_rep1.tsv.gz', 'rep2=fragments_rep2.tsv.gz'),
    sampleNames=c('rep1', 'rep2'),
    minTSS=4, minFrags=1000,
    addTileMat=TRUE, addGeneScoreMat=TRUE)

# 2. Doublet detection (built-in)
doubletScores <- addDoubletScores(input=ArrowFiles, k=10, knnMethod='UMAP')

# 3. Project + filter
proj <- ArchRProject(ArrowFiles=ArrowFiles, outputDirectory='ArchR_out')
proj <- filterDoublets(proj)

# 4. LSI + UMAP + clustering
proj <- addIterativeLSI(proj, useMatrix='TileMatrix', name='IterativeLSI')
proj <- addClusters(proj, reducedDims='IterativeLSI', method='Seurat', resolution=0.5)
proj <- addUMAP(proj, reducedDims='IterativeLSI')

# 5. Reproducible peakset (per cluster)
proj <- addGroupCoverages(proj, groupBy='Clusters')
proj <- addReproduciblePeakSet(proj, groupBy='Clusters', pathToMacs2='/path/macs3')
proj <- addPeakMatrix(proj)
```

## SnapATAC2 Workflow (Python)

```python
import snapatac2 as snap

# 1. Load 10X fragments. SnapATAC2 uses snap.pp.import_fragments (NOT snap.read_10x, which doesn't exist).
data = snap.pp.import_fragments(
    fragment_file='outs/fragments.tsv.gz',
    chrom_sizes=snap.genome.hg38,
    file='out.h5ad',                        # backed AnnData; backend handled by file path
    sorted_by_barcode=False)
data.obs['sample_id'] = 'rep1'

# 2. Per-cell QC
snap.metrics.tsse(data, gene_anno=snap.genome.hg38)
snap.metrics.frag_size_distr(data)
snap.pp.filter_cells(data, min_counts=1000, min_tsse=4)

# 3. Tile matrix + spectral
snap.pp.add_tile_matrix(data, bin_size=500)
snap.pp.select_features(data, n_features=250000)
snap.tl.spectral(data)
snap.tl.umap(data)
snap.tl.leiden(data)

# 4. Per-cluster peak calling (uses MACS3)
snap.tl.macs3(data, groupby='leiden')

# 5. Gene activity (gene score) for annotation
gene_mat = snap.pp.make_gene_matrix(data, gene_anno=snap.genome.hg38)
```

## Reconciliation Across Tools

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Signac UMAP shows tight clusters; ArchR UMAP shows loose | Different LSI implementation; ArchR uses iterative LSI by default | Both valid; biology should match in cluster labels |
| Different tools call different doublet rates | Different algorithms (collision vs simulation) | Use intersection (cells flagged by 2+ tools) as high-confidence doublets |
| Cluster boundaries differ | Different clustering algorithm or resolution | Standardize on Leiden algorithm 4 with same resolution |
| Per-cluster peak count differs | Different peak callers or pseudobulk depth | Ensure same MACS3 parameters; pool small clusters |

**Operational rule:** For high-confidence cell-type annotation, agree across two ecosystems (e.g., Signac + ArchR clusters) and report agreement metrics.

## cellranger-atac vs cellranger-arc

| Pipeline | Use for | Output |
|----------|---------|--------|
| cellranger-atac | Single-modality 10X scATAC (chemistry v1, v2) | per-cell ATAC barcodes; fragments.tsv.gz; per-barcode metadata |
| cellranger-arc | 10X Multiome (paired RNA + ATAC same cell, "Multiome chemistry") | joint barcodes for RNA + ATAC; separate fragment / count files; one barcode whitelist |

**Trigger:** Loading 10X output without checking which pipeline produced it.

**Mechanism:** cellranger-atac and cellranger-arc produce different output structures. cellranger-arc fragments.tsv has barcodes paired with the RNA matrix; cellranger-atac fragments are ATAC-only with their own barcode universe.

**Fix:** Verify the chemistry on the cellranger summary (look for "Multiome" in the run config). Use `Read10X_h5` for cellranger-arc Multiome RNA output; use `CreateChromatinAssay` with the matched fragments for the ATAC. Mixing barcodes across pipelines fails silently.

## Cell-Cycle Correction for scATAC

**Trigger:** Proliferating cell types in the dataset; cells distributed across cell cycle phases.

**Mechanism:** Replication-associated chromatin opening adds 5-10% global accessibility per cell as it moves G1 -> S -> G2/M; this confounds clustering and DA.

**Detection:** Score cells with a chromatin-adapted S-phase signature (regulated origin loci, replication-stress-response gene accessibility) analogous to Seurat's CellCycleScoring on RNA. Or compute a Repli-seq peak overlap score.

**Fix:** Regress S-phase score on the TF-IDF residuals before downstream LSI: `ScaleData(obj, vars.to.regress='S.score')` analogous to RNA workflow. For DA between cell-cycle-mismatched conditions, add S-phase as covariate in pseudobulk DESeq2.

## Sex-Chromosome QC for scATAC

**Trigger:** Mixed-sex donors in the dataset; sample-swap detection.

**Mechanism:** XIST locus accessibility is high in female cells (X-inactivation); chrY peak count is essentially zero in females. Per-cell or per-sample, the XIST/chrY ratio identifies sex unambiguously.

**Detection:** In a per-cell counts matrix, compute fraction of fragments at XIST locus (chrX:73820651-73852723 hg38) and chrY peak count; classify cells; flag cells/samples where assignment disagrees with metadata.

**XCI escapees:** Genes that escape X-inactivation (KDM6A, DDX3X, EIF1AX) are biallelically accessible in female cells but not male; can be used as additional sex confirmation if XIST is ambiguous.

## scArches Reference Mapping

**Trigger:** Projecting query scATAC onto a reference atlas; cross-study integration without batch effects.

**Mechanism:** scArches (Lotfollahi 2022) provides transfer-learning to project a new dataset onto an existing reference's latent space without retraining the reference. For ATAC, the relevant model is PEAKVI (Ashuach 2022). PEAKVI lives in `scvi-tools` (`scvi.model.PEAKVI`); the `load_query_data` classmethod implements the scArches algorithm directly, so an explicit `import scarches` is not required.

```python
import scvi

# Pre-trained reference model (e.g. PBMC scATAC atlas)
query_model = scvi.model.PEAKVI.load_query_data(adata_query, reference_path)
query_model.train(max_epochs=200)
adata_query.obsm['X_emb'] = query_model.get_latent_representation()
```

Reference atlases for ATAC are still emerging; the most-developed are PBMC (Granja 2021) and brain (BRAIN Initiative).

## chromBPNet Per-Cluster Pseudobulk

**Trigger:** Cell-type-specific variant effect prediction; per-cluster bias-corrected calling.

**Mechanism:** Train chromBPNet (atac-seq/deep-learning-atac) per pseudobulk cluster; outputs are bias-corrected per-base profiles + variant effect predictions specific to that cell type.

**Workflow:** Aggregate fragments per cluster into pseudobulk BAMs; run chromBPNet pipeline per cluster (~24h GPU per cluster); use the resulting model for in silico variant scoring at GWAS / rare-variant SNPs in that cell type.

For most studies, this is reserved for top 5-10 priority clusters; running chromBPNet per cluster on >20 cell types is computationally heavy.

## AMULET Depth Threshold

AMULET is most reliable at >= 5,000 fragments per cell. Below 5,000 fragments, the collision-based detection has insufficient power: the expected number of collisions per cell is too small for the binomial test to distinguish doublet from singleton. For < 5,000 fragment libraries, use synthetic-doublet methods (ArchR addDoubletScores, scDblFinder) instead.

## Multiome WNN Integration (Signac)

```r
# Assume `obj` has both 'RNA' and 'ATAC' assays from same Multiome experiment
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(obj) <- 'ATAC'
obj <- RunTFIDF(obj) %>% FindTopFeatures(min.cutoff='q0') %>% RunSVD()

# Joint embedding
obj <- FindMultiModalNeighbors(obj, reduction.list=list('pca', 'lsi'),
                               dims.list=list(1:30, 2:30))
obj <- RunUMAP(obj, nn.name='weighted.nn', reduction.name='wnn.umap')
```

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| UMAP shows depth gradient | LSI component 1 included in clustering | `dims=2:30` instead of `1:30` |
| Cell Ranger output has many low-quality cells | cellranger ATAC uses lenient cell calling | Re-filter at fragment count >= 1000 + TSS enrichment >= 4 |
| ArchR "TileMatrix not found" | Forgot `addTileMat=TRUE` in createArrowFiles | Re-create Arrow files with the flag |
| Signac CreateChromatinAssay fails on missing fragments file | Path to fragments.tsv.gz incorrect; or missing tabix index | Provide full path; run `tabix -p bed fragments.tsv.gz` |
| MACS3 fails on tiny pseudobulk | Cluster too small | Use cluster aggregation; require >= 200 cells per cluster |
| AMULET reports 100% doublets | Threshold mis-set or input is technical replicates | Check fragment-count distribution; AMULET works at >= 5000 fragments/cell |
| Multiome WNN clusters dominated by ATAC noise | Equal modality weighting | Inspect modality weights; manually adjust if needed |
| chromVAR / motif assay all NA | Run before peakset finalized | Re-run AddMotifs / RunChromVAR after peaks stable |
| EnsDb / BSgenome version mismatch | hg38 BSgenome with wrong-build EnsDb | Match builds; `EnsDb.Hsapiens.v86` is GRCh38 (Ensembl v86); use `EnsDb.Hsapiens.v75` for hg19. Newer hg38 EnsDb releases (v98+) reflect newer GENCODE annotations |

## References

- Stuart T et al 2021 Nat Methods 18:1333 (Signac)
- Granja JM et al 2021 Nat Genet 53:403 (ArchR)
- Zhang K et al 2024 Nat Methods 21:217 (SnapATAC2)
- Hao Y et al 2021 Cell 184:3573 (Seurat WNN)
- Thibodeau A et al 2021 Genome Biol 22:252 (AMULET doublet detection)
- Germain PL et al 2021 F1000Res 10:979 (scDblFinder)
- Cusanovich DA et al 2015 Science 348:910 (sciATAC; LSI for sc data)
- Chen H et al 2019 Genome Biol 20:241 (scATAC analysis benchmark)
- Heumos L et al 2023 Nat Rev Genet 24:550 (Best practices for single-cell)
- 10X Genomics Cell Ranger ATAC documentation

## Related Skills

- atac-seq/atac-qc - Bulk QC patterns adapted for per-cell
- atac-seq/atac-peak-calling - Pseudobulk peak calling per cluster
- atac-seq/consensus-peakset - Across-cluster consensus
- atac-seq/differential-accessibility - Pseudobulk DA per cluster
- atac-seq/motif-deviation - chromVAR for per-cell TF activity
- atac-seq/footprinting - scprinter for sc footprinting
- atac-seq/co-accessibility - Cicero for cis-regulatory connections
- atac-seq/deep-learning-atac - chromBPNet / scBasset for per-cluster bias correction and variant effects
- atac-seq/enhancer-gene-linking - Per-cell-type enhancer-gene maps
- atac-seq/allele-specific-accessibility - sc allelic imbalance for cis-effects
- single-cell/preprocessing - General sc QC patterns
- single-cell/clustering - Cluster definition
- single-cell/cell-annotation - Marker-based annotation
- single-cell/multimodal-integration - Multiome RNA+ATAC integration
- single-cell/scatac-analysis - Cross-reference single-cell ATAC-specific patterns
- single-cell/batch-integration - scArches reference mapping; Harmony
