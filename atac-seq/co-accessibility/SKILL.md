---
name: bio-atac-seq-co-accessibility
description: Infer cis-regulatory connections (peak-to-peak co-accessibility) from scATAC-seq using Cicero, ArchR getCoAccessibility, or SCENIC+. Use when linking enhancer accessibility to promoter accessibility, identifying enhancer-gene pairs from chromatin alone (without paired RNA), running gene-regulatory inference combining ATAC + RNA, or comparing predicted regulatory contacts against Hi-C/Micro-C ground truth.
tool_type: r
primary_tool: cicero
---

## Version Compatibility

Reference examples tested with: Cicero 1.20+, monocle3 1.3+, ArchR 1.0.2+, SCENIC+ 1.0+, pycisTopic 1.0+, Signac 1.13+, GenomicRanges 1.54+, GenomicInteractions 1.36+, BSgenome.Hsapiens.UCSC.hg38 1.4+.

Verify before use:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# Co-accessibility (cis-Regulatory Linkage)

**"Which enhancers connect to which promoters in my scATAC data?"** -> Use cell-to-cell variability in joint accessibility of nearby peaks to infer cis-regulatory connections without explicit RNA expression. Output is a peak-pair graph with co-accessibility scores; thresholding produces enhancer-gene candidate pairs.

- R: `cicero::run_cicero(input_cds, genomic_coords)` -> peak-pair connection scores
- R: `ArchR::addCoAccessibility(proj)` -> ArchR-internal Cicero wrapper
- Python: `pycisTopic` + `SCENIC+` for network-level inference combining ATAC + RNA + motifs

Co-accessibility is NOT 3D contact; it's a statistical association based on cell-to-cell co-variation. Strong co-accessibility correlates with Hi-C/Micro-C contacts (~30-50% concordance) but is not equivalent.

## What Co-accessibility Captures vs What It Doesn't

| Captures | Misses |
|----------|--------|
| Peak pairs that vary together across cell states | 3D physical contacts that don't vary in accessibility |
| Cis-regulatory grammar within a cell type | Trans-chromosomal interactions |
| Active enhancer-promoter pairs | Constitutive structural contacts |
| Lineage-specific regulation | Developmental contacts that opened before scATAC sample |
| Distance-decay biology of enhancer-promoter | Hub enhancers that contact many distal targets |

For physical contact, use Hi-C, Micro-C, or PCHi-C. Co-accessibility is the chromatin-only proxy.

## Algorithmic Taxonomy

| Tool | Method | Input | Output | Strength | Fails when |
|------|--------|-------|--------|----------|------------|
| Cicero (Pliner 2018) | Graphical lasso on aggregated cell metacells | scATAC peak-cell matrix + cell trajectory | Peak-pair connection score (0-1) | Original, well-validated; integrates with Monocle3 | Slow on >50K cells; sensitive to alpha tuning |
| ArchR getCoAccessibility | Cicero-based; uses ArchR's metacell aggregation | ArchR project | Same as Cicero | Built-in to ArchR pipeline; faster on large datasets | Tied to ArchR; same biology as Cicero |
| SCENIC+ (Bravo 2023) | Multi-step: co-accessibility + motif scoring + RNA correlation | Multiome (ATAC + RNA) or paired | TF-driven enhancer-gene networks | Most comprehensive; multi-modal | Multiome data required; computationally heavy |
| LinkPeaks (Signac) | Pearson correlation of accessibility with paired gene expression | Multiome | Peak-gene linkage score | Direct enhancer-gene from RNA correlation | Multiome-only; not pure ATAC |
| GeneHancer / FANTOM5 / EpiMap | Bulk-derived enhancer-gene reference | None (database lookup) | Pre-computed enhancer-gene pairs | Comprehensive; published references | Cell-type-agnostic; may not match the biology of interest |

Methodology evolves; verify against Pliner 2018 (Cicero), Bravo 2023 (SCENIC+), Nasser 2021 (ABC model alternative for enhancer-gene), and current Hi-C concordance benchmarks.

## How Cicero Works (Conceptually)

Cell-to-cell variability is too sparse for direct correlation. Cicero solves this via metacells:

1. Reduce dimensionality (UMAP from input).
2. Build k-NN graph of cells.
3. Aggregate k cells into metacells (default k = 50).
4. Compute correlation in accessibility across metacells, restricted to peak pairs within `genomic_distance_max` (default 500 kb cis).
5. Apply graphical lasso with regularization `alpha` to sparsify the correlation matrix.
6. Output: per-pair connection score; positive = co-variation, negative = anti-co-variation.

Connection thresholds typically 0.05-0.5; > 0.25 is high-confidence.

## Per-Tool Failure Modes

### Cicero -- alpha tuning shifts results

**Trigger:** Default alpha (sometimes computed automatically from data); custom alpha < 0.5 or > 5.

**Mechanism:** Alpha controls graphical lasso regularization. Too low: dense graph with many spurious connections; too high: sparse with biology missing.

**Symptom:** Connection count varies 10-100x across alpha sweeps.

**Fix:** Use Cicero's `estimate_distance_parameter()` to get data-driven alpha; verify connection count is biologically plausible (~10-50% of peaks have at least one strong connection).

### Cicero -- metacell aggregation hides cell-type-specific connections

**Trigger:** Running Cicero on heterogeneous dataset spanning multiple cell types.

**Mechanism:** Metacells aggregate across cell types; connections that exist only in one cell type get diluted.

**Fix:** Run Cicero per-cluster separately; combine results with cluster annotations. Cell-type-specific connections often differ.

### Cicero -- distance assumption

**Trigger:** Default `genomic_distance_max=500000` (500 kb cis only).

**Mechanism:** Distal connections beyond 500 kb cis are excluded; trans-chromosomal entirely missed.

**Fix:** For specific use cases (e.g., gene desertless TADs), increase `genomic_distance_max` to 1 Mb or more. Trans connections require Hi-C, not co-accessibility.

### SCENIC+ -- RNA scaling

**Trigger:** RNA-side dropouts in Multiome data.

**Mechanism:** SCENIC+ requires reasonable RNA quantification per cell. Sparse Multiome RNA with many zero genes causes correlation degradation.

**Fix:** Filter cells with insufficient RNA; aggregate cells if necessary. Multiome RNA should look comparable to standalone scRNA-seq.

### LinkPeaks (Signac) -- Distance default

**Trigger:** Default `LinkPeaks(..., distance=5e+05)`.

**Mechanism:** Same as Cicero; 500 kb cis only by default.

**Fix:** Same; widen if needed but trans not supported.

## Decision Tree by Goal

| Goal | Tool |
|------|------|
| ATAC-only enhancer-promoter inference | Cicero |
| ATAC-only inside ArchR ecosystem | ArchR getCoAccessibility |
| Multiome (RNA + ATAC) enhancer-gene inference | LinkPeaks (Signac) for direct correlation; SCENIC+ for TF network |
| TF-driven regulatory networks | SCENIC+ (requires Multiome) |
| Comparison against Hi-C / Micro-C | Cicero output -> overlap with HiCCUPS loops |
| Published reference enhancer-gene pairs | GeneHancer, FANTOM5, EpiMap (pre-computed lookup) |
| Gene desertless distal regulation | Cicero with widened distance; or H3K27ac HiChIP |

## Cicero Standard Workflow

```r
library(cicero); library(monocle3); library(GenomicRanges)

# Input: peak-cell binary matrix from Signac/ArchR (rows = peaks, cols = cells)
# Convert peaks to "chrN_start_end" format
peak_names <- paste0(seqnames(peaks), '_', start(peaks), '_', end(peaks))
input_cds <- new_cell_data_set(peak_matrix, cell_metadata=metadata,
                               gene_metadata=peak_metadata)

# Reduce dimensionality (UMAP from input)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method='LSI')
input_cds <- reduce_dimension(input_cds, reduction_method='UMAP',
                              preprocess_method='LSI')

# Build metacell-aggregated CDS
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords, k=50)

# Run Cicero with hg38 chrom sizes
genome_df <- data.frame(chr=seqnames(seqinfo(BSgenome.Hsapiens.UCSC.hg38)),
                        length=seqlengths(seqinfo(BSgenome.Hsapiens.UCSC.hg38)))
conns <- run_cicero(cicero_cds, genomic_coords=genome_df,
                    window=500000, sample_num=100)

# Filter to high-confidence connections
strong <- conns[conns$coaccess > 0.25, ]
cat(sprintf('Total conns: %d; strong (>0.25): %d\n', nrow(conns), nrow(strong)))
```

## ArchR getCoAccessibility

```r
library(ArchR)
proj <- loadArchRProject('ArchR_out')
proj <- addCoAccessibility(proj, reducedDims='IterativeLSI',
                          k=100, knnIteration=500,
                          maxDist=250000)               # 250 kb cis (tighter than default)
co_acc <- getCoAccessibility(proj, corCutOff=0.5,       # Default 0.5; lower for more
                             returnLoops=FALSE)
```

ArchR returns connections as a GRanges object compatible with `GenomicInteractions` for direct overlap with Hi-C loops.

## Visualizing Connections

```r
# As arc plot at a locus of interest
library(Gviz); library(GenomicInteractions)
gi <- GenomicInteractions(anchorOne=GRanges(strong[, 1:3]),
                          anchorTwo=GRanges(strong[, 4:6]),
                          counts=as.integer(strong$coaccess * 100))
plotInteractions(gi, view='arc', interactions.color='red')
```

For genome-browser visualization with ArchR: `plotPeak2GeneHeatmap()` shows the peak-gene linkage matrix; `plotBrowserTrack()` overlays connections on tracks.

## SCENIC+ TF-Driven Networks

```python
import scenicplus as sp

# Multi-step: requires preprocessed pycisTopic models + paired RNA AnnData
scplus_obj = sp.create_SCENICPLUS_object(
    GEX_anndata=rna_adata,
    cisTopic_obj=topic_obj,
    menr=motif_enr_dict)

# Calculate eRegulons (TF + target genes + linked enhancers)
sp.run_scenicplus(scplus_obj, ...)
scplus_obj.uns['eRegulon_metadata']
```

SCENIC+ is significantly more complex than Cicero; budget 1-2 days for setup. The benefit is that outputs are TF -> enhancer -> gene triples, not just peak-peak co-accessibility.

## Cicero Alpha Mathematics

**Trigger:** Tuning Cicero's regularization parameter for the graphical lasso step.

**Mechanism:** `estimate_distance_parameter()` fits a regression of pairwise correlation versus genomic distance for nearby peaks; the slope gives the expected co-accessibility decay with distance. Cicero's alpha penalty is set proportional to this slope so the graph sparsifies at biologically appropriate distance scales.

**Implementation:** Cicero internally calls `estimate_distance_parameter(cicero_cds, window=window, maxit=100, sample_num=100)` -- the function bootstraps 100 windows, fits per-window distance-correlation regression, and averages the alpha estimates. Default value is robust at typical peak densities; manual override is rarely needed.

**When manual tuning helps:** Very dense peaksets (>200k peaks) may need higher alpha to control false positives; very sparse (<10k peaks) may need lower alpha to recover signal. Verify by running on technical replicates -- the expected outcome is ~0 strong connections.

## ABC Model Cross-Reference

For enhancer-to-gene linking with paired Hi-C/Micro-C, the canonical method is the ABC model (Fulco 2019, Nasser 2021), not Cicero. ABC computes ABC = (Activity_E * Contact_E,G) / sum_e(Activity_e * Contact_e,G); standardizes on combined ATAC + H3K27ac activity and Hi-C contact frequencies. ENCODE-rE2G (2024) is the modern logistic-regression replacement.

See atac-seq/enhancer-gene-linking for full ABC and ENCODE-rE2G coverage. Cicero is the ATAC-only fallback when no Hi-C is available.

## HiChIP H3K27ac as Orthogonal Anchor

| Decision | Action |
|----------|--------|
| Have Hi-C / Micro-C | Use ABC (atac-seq/enhancer-gene-linking) primary; Cicero as ATAC-only sanity check |
| Have HiChIP H3K27ac | FitHiChIP loops (FDR < 0.05, count >= 5) primary; ABC + HiChIP intersection is high-confidence |
| Have ATAC + H3K27ac, no 3D | ABC with average HiC fallback (Fulco 2019); document degraded performance |
| Have only ATAC | Cicero (this skill); known concordance with Hi-C ~30-50% |

Cicero is appropriate when no 3D data exists; do not use Cicero in lieu of ABC when Hi-C/Micro-C are available.

## Hi-C / Micro-C Concordance

| Hi-C concordance | Action |
|-----------------|--------|
| > 50% of strong Cicero connections overlap Hi-C loops | High-confidence; Cicero captures real 3D structure |
| 30-50% | Standard; some 3D contacts don't vary in accessibility |
| < 20% | Co-accessibility may not reflect contacts; lineage-specific contacts may be missing |

```r
# Compare Cicero against published Hi-C loops
library(GenomicInteractions)
hic_loops <- import('hiccups_loops.bedpe', format='bedpe')
ci <- GenomicInteractions(anchorOne=GRanges(strong$Peak1),
                          anchorTwo=GRanges(strong$Peak2))
overlap <- countOverlaps(ci, hic_loops, type='equal') > 0
cat(sprintf('Cicero connections overlapping HiCCUPS loops: %.1f%%\n',
            100 * mean(overlap)))
```

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Cicero many weak connections; ArchR few strong | Different alpha or aggregation | Standardize parameters |
| LinkPeaks (Multiome) finds connections Cicero misses | LinkPeaks uses RNA expression as the anchor; Cicero is ATAC-only | Both valid; report intersection as high-confidence |
| Co-accessibility doesn't match Hi-C in heterochromatin | Heterochromatic contacts are constitutive; co-accessibility needs variation | Expected; co-accessibility complements Hi-C |
| SCENIC+ network has ENCODE-validated TFs but missing some | Motif database limited or RNA imputation missed | Expand motif database; integrate paired ChIP-seq if available |

**Operational rule:** Co-accessibility is a hypothesis generator. Validate with Hi-C, ChIP-seq, or experimental enhancer-promoter interaction (CRISPRi-FlowFISH).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Cicero `make_cicero_cds` slow / crashes | k too high or cell count too large | Reduce k or subsample cells |
| All connections near zero | alpha set too high | Use `estimate_distance_parameter()` |
| Connection score > 1 reported | Bug in older Cicero versions | Update; check `as.numeric(coaccess)` for outliers |
| ArchR getCoAccessibility "TileMatrix" error | Need PeakMatrix not TileMatrix | `addPeakMatrix()` first |
| SCENIC+ install fails | Many heavy dependencies | Use the published Docker image |
| Connection count varies wildly per run | Stochastic metacell aggregation | Set seed; or aggregate at higher k for stability |
| LinkPeaks all NaN | RNA expression has too many zeros | Re-filter cells with sufficient RNA |
| Peak names not matching | format mismatch (chr_start_end vs chr:start-end) | Standardize naming convention |

## References

- Pliner HA et al 2018 Mol Cell 71:858 (Cicero)
- Granja JM et al 2021 Nat Genet 53:403 (ArchR getCoAccessibility)
- Bravo Gonzalez-Blas C et al 2023 Nat Methods 20:1355 (SCENIC+)
- Stuart T et al 2021 Nat Methods 18:1333 (Signac LinkPeaks)
- Nasser J et al 2021 Nature 593:238 (ABC model; alternative enhancer-gene)
- Fulco CP et al 2019 Nat Genet 51:1664 (CRISPRi-FlowFISH; gold-standard validation)
- Mumbach MR et al 2017 Nat Genet 49:1602 (HiChIP H3K27ac for enhancer-promoter)
- Boix CA et al 2021 Nature 590:300 (EpiMap; bulk enhancer-gene reference)

## Related Skills

- atac-seq/single-cell-atac - scATAC preprocessing (input)
- atac-seq/consensus-peakset - Peak set used for connection inference
- atac-seq/motif-deviation - chromVAR for TF activity (complement)
- atac-seq/enhancer-gene-linking - ABC, ENCODE-rE2G, CRISPRi-FlowFISH validation when Hi-C is available
- atac-seq/deep-learning-atac - chromBPNet variant effect at predicted enhancers
- gene-regulatory-networks/scenic-regulons - Standalone SCENIC for TF networks
- hi-c-analysis/loop-calling - Physical contacts from Hi-C
- hi-c-analysis/contact-pairs - Hi-C / Micro-C contact pairs
- single-cell/multimodal-integration - Multiome integration
- chip-seq/peak-annotation - Cross-validate with TF ChIP
- pathway-analysis/gsea - Downstream gene-level enrichment
