---
name: bio-crispr-screens-perturb-seq-analysis
description: Analyzes single-cell pooled CRISPR screens (Perturb-seq, CROP-seq, Perturb-CITE-seq, ECCITE-seq, multiome) where each cell carries an sgRNA and a scRNA-seq / surface-protein / chromatin readout. Covers experimental design (direct-capture Perturb-seq Dixit 2016 vs CROP-seq 3'UTR-barcoded Datlinger 2017 vs ECCITE-seq vs Multiome), MOI for sgRNA assignment, escaper-cell filtering (Mixscape, Papalexi 2021), SCEPTRE NB GLM + permutation for low-MOI (Barry 2024 Genome Biol 25:124), the Pertpy framework, factor decomposition, genome-scale Perturb-seq (Replogle 2022 Cell, 2.5M cells), and per-perturbation single-cell DE. Use when running a single-cell CRISPR screen, choosing direct-capture vs CROP-seq architecture, filtering escaper cells, performing single-cell DE, integrating Perturb-seq with pathway analysis, scaling to GW CRISPRi via Replogle protocol, or analyzing multi-omics screens.
tool_type: python
primary_tool: Pertpy
---

## Version Compatibility

Reference examples tested with: Pertpy 0.6+, SCEPTRE 0.10+ (R / katsevich-lab/sceptre), Mixscape via Seurat 4.3+ or Pertpy, scanpy 1.10+, anndata 0.10+, pandas 2.2+, numpy 1.26+, scipy 1.12+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show pertpy scanpy anndata`
- R: `packageVersion('sceptre')`; `?sceptre`; `?Seurat::PrepLDA`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Single-Cell Perturb-Seq Analysis

**"Analyze a single-cell pooled CRISPR perturbation screen"** -> Assign sgRNAs to cells, filter unperturbed escapers, normalize counts, fit per-gene differential expression conditioned on perturbation, and rank perturbations by their molecular effect.

- Python: `pertpy` unified framework for Mixscape + SCEPTRE-via-R + differential expression
- R: `sceptre` for low-MOI NB GLM + permutation testing
- Python/R: `Seurat::MixscapeLDA` and downstream

## Experimental Architecture Comparison

| Method | Year | Architecture | Readout | MOI | Single-cell sgRNA detection |
|--------|------|--------------|---------|-----|------------------------------|
| Perturb-seq (Dixit 2016, *Cell*) | 2016 | sgRNA expressed in cassette; direct PCR capture | scRNA-seq | High (4-50 sgRNAs/cell) | Yes via amplicon-PCR pre-sequencing |
| CROP-seq (Datlinger 2017, *Nat Methods*) | 2017 | sgRNA expressed via tRNA-spacer-sgRNA-pgk-puro; barcoded in 3'UTR of fluorescent marker; captured in 3' droplet | scRNA-seq | Low (1-2 sgRNAs/cell) | Native via 10X 3' chemistry |
| Perturb-CITE-seq (Frangieh 2021, *Nat Genet*) | 2021 | Adds surface-protein hashtag oligos to CROP-seq | scRNA-seq + ADT (protein) | Low | CROP-seq architecture |
| ECCITE-seq (Mimitou 2019, *Nat Methods*) | 2019 | Surface-protein hashtag with sgRNA-marked cells | scRNA-seq + ADT | Low | Hash + sgRNA |
| Perturb-ATAC (Rubin 2019, *Cell*) | 2019 | scATAC-seq readout | scATAC | Low | sgRNA capture via separate library prep |
| Perturb-multiome (10X) | 2021+ | scRNA + scATAC simultaneously | scRNA + ATAC | Low | Direct capture from sgRNA cassette |
| Replogle GW Perturb-seq (2022, *Cell*) | 2022 | Multiplexed CRISPRi with sgRNA barcoding | scRNA-seq | 1 sgRNA/cell | Direct capture |

**Decision rule:** For genome-wide CRISPRi screens, Replogle's CRISPRi + 10X 3' direct-capture protocol is the gold standard (2.5M cells across 2,058 perturbations in Replogle 2022). For protein readout, Perturb-CITE-seq. For chromatin, Perturb-multiome. For low-throughput pilot, original Dixit Perturb-seq.

## MOI and sgRNA Assignment

**The central technical challenge:** Each cell must receive exactly one sgRNA (otherwise the perturbation is undefined). At MOI 0.3, ~26% of cells get ≥1 sgRNA, but 4% get ≥2; the cells with multiple sgRNAs must be filtered or analyzed as combinatorial perturbations.

**Assignment workflow:**

1. **Detect sgRNA reads per cell:** From the sgRNA library prep (direct capture or 3'UTR barcode), count reads per sgRNA per cell.
2. **Threshold:** Most pipelines use 10+ reads of one sgRNA to assign that perturbation.
3. **Multiplets:** Cells with 2+ sgRNAs at >10 reads each are either multi-perturbed (analyzable as combinatorial) or doublets.
4. **Doublet detection:** Use scDblFinder, Scrublet, or AMULET (multiome) to identify doublets independently from sgRNA assignment.

**Goal:** Assign a single perturbation identity (or 'multiplet'/'none') to every cell from the sgRNA counts matrix.

**Approach:** Threshold per-cell sgRNA reads at ≥10 (Pertpy convention); cells exceeding the threshold for exactly one sgRNA are assigned that perturbation; cells with multiple sgRNAs above threshold are flagged as multiplets for filtering or combinatorial analysis.

```python
# sgRNA assignment via threshold counting
def assign_sgrna(adata, sgrna_counts_layer='sgrna_counts', threshold=10):
    '''Per-cell sgRNA assignment. Returns single assignment or 'multiplet'/'none'.'''
    import numpy as np
    counts = adata.layers[sgrna_counts_layer]  # cells x sgRNAs
    above_thresh = counts >= threshold
    n_sgrna_per_cell = above_thresh.sum(axis=1)
    assignments = np.where(
        n_sgrna_per_cell == 0, 'none',
        np.where(n_sgrna_per_cell == 1,
                  [adata.var_names[i] for i in counts.argmax(axis=1)],
                  'multiplet'))
    adata.obs['sgrna_assignment'] = assignments
    return adata
```

## Escaper Cell Filtering (Mixscape)

**Why this matters:** Not all sgRNA-positive cells actually edit. 30-60% of cells may receive sgRNA but fail to undergo gene knockdown ("escapers"). Including escapers dilutes the perturbation effect; Mixscape (Papalexi 2021) identifies and filters them.

**Mixscape algorithm:** For each perturbed cell, compute a "perturbation signature" = (its expression) - (mean of K nearest non-targeting-control cells). This signature isolates the perturbation effect from cell-state variation. Cells with perturbation signature similar to NTC distribution are escapers.

```python
import pertpy as pt
import scanpy as sc

# adata is a scRNA-seq AnnData with 'sgrna_assignment' column
# Pertpy 0.6+ Mixscape API (verify against installed pertpy with help(pt.tl.Mixscape))
mixscape = pt.tl.Mixscape()
mixscape.perturbation_signature(
    adata=adata,
    pert_key='sgrna_assignment',     # .obs column with sgRNA target per cell
    control='NTC',                   # name of non-targeting control in pert_key column
    n_neighbors=20,                  # K neighbors for KNN-NTC subtraction
)
# Writes .layers['X_pert'] with perturbation-signature-corrected expression

# Filter escapers: classify perturbed cells as KO (true perturbation) or NP (non-perturbed/escaper)
mixscape.mixscape(
    adata=adata,
    pert_key='sgrna_assignment',     # pert_key (not 'labels' in modern pertpy)
    control='NTC',
    new_class_name='mixscape_class', # .obs column to write
)
# Defaults to layer='X_pert' (output of perturbation_signature)

# Keep only KO cells for downstream analysis
adata_ko = adata[adata.obs['mixscape_class'].isin(['KO'])].copy()
print(f'KO cells: {adata_ko.n_obs} ({adata_ko.n_obs/adata.n_obs:.1%} of perturbed)')
```

**Critical:** Mixscape can fail when the perturbation has weak phenotype; empirically Mixscape detects perturbations with log-fold-change <-0.5 (depletion) reliably, but weaker effects collapse into the NTC distribution. For genome-wide screens, run Mixscape per perturbation; for low-effect perturbations, trust the assignment without filtering.

## SCEPTRE for Low-MOI Differential Expression

**Why this matters:** Standard differential-expression tools (DESeq2, MAST) assume Gaussian-mixture distribution and fail at single-cell scale with sparse, zero-inflated data. SCEPTRE (Katsevich Lab, 2021; low-MOI variant Barry 2024 Genome Biol) uses a negative-binomial GLM with conditional resampling:

1. Per gene, fit NB GLM: `log(expr_g) ~ pert_indicator + technical_factors`
2. Compute z-score for the perturbation coefficient
3. Resample the pert_indicator (conditional on counts) 500-1000 times; compute permutation null
4. Get FDR via permutation; not parametric

```r
library(sceptre)

# Input: sce object or sparse matrix + metadata
# Required: gene_expression_matrix, perturbation_indicator (binary per cell per pert),
#           technical_factors (batch, n_genes, etc.)

# For each gene + perturbation pair:
results <- run_sceptre_low_moi(
    expression_matrix    = gene_expr,            # genes x cells
    grouping_var         = pert_indicator,       # length(cells); 0 = NTC, 1 = perturbed
    technical_factors    = covariates_df,        # n_genes, n_umi, batch
    response_grouping_var = 'gene_id',
    n_permutations       = 1000
)
# Output: per-gene-per-pert p-value, log-fold-change, FDR
```

**Advantage over MAST:** SCEPTRE's permutation NB GLM is the only method that maintains calibrated FDR in pooled-screen scRNA-seq (Barry 2024 benchmark). MAST and Wilcoxon are over-confident due to data sparsity.

## Pertpy Unified Framework

**Pertpy** (https://pertpy.readthedocs.io) integrates Mixscape, distance-based perturbation comparison, EdgeR/PyDESeq2/WilcoxonTest DE, and factor models in a single AnnData-based interface. For SCEPTRE specifically, invoke the R sceptre package separately (Pertpy does not wrap it).

```python
import pertpy as pt
import scanpy as sc

# Load data
adata = pt.dt.papalexi_2021()  # built-in example from Mixscape paper

# Standard scRNA-seq preprocessing (scanpy)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Mixscape escaper filtering (writes .layers['X_pert'] and .obs['mixscape_class'])
ms = pt.tl.Mixscape()
ms.perturbation_signature(adata, pert_key='perturbation', control='NT', n_neighbors=20)
ms.mixscape(adata, pert_key='perturbation', control='NT')

# Filter to KO cells
adata_ko = adata[adata.obs['mixscape_class'].isin(['KO', 'NT'])].copy()

# Pseudobulk differential expression via pertpy (PyDESeq2 backend)
de = pt.tl.PyDESeq2(adata_ko, design='~perturbation')
de.fit()
results_df = de.test_contrasts(contrast=('perturbation', 'GENE_X', 'NT'))

# For calibrated SCEPTRE on low-MOI single-cell data, use R sceptre directly
# (Barry 2024 Genome Biol; not bundled in pertpy)
```

## Genome-Wide Perturb-Seq (Replogle 2022)

**Replogle 2022 *Cell* 185:2559** demonstrated genome-wide Perturb-seq:
- 2.5M cells across 2,058 perturbations (essentially every expressed protein-coding gene)
- CRISPRi via dCas9-KRAB
- Native 10X 3' direct-capture for sgRNA
- ~1,000 cells per perturbation (enough for per-perturbation DE)
- Cluster-based analysis of perturbed cells reveals gene-program organization

**Scaling principles:**
- Cells per perturbation: 500-1,000 minimum for stable DE
- 10X channels: 10-30 channels at 5,000-10,000 cells each
- Cost: ~$50-100K for genome-scale

```python
# Replogle-style genome-wide design
# Each cell -> 1 sgRNA (low MOI)
# Each gene -> 5 sgRNAs (CRISPRi)
# Each pert -> 1000 cells
# Total: 19,000 genes x 5 sgRNAs = 95,000 sgRNAs
# Cells needed: 2 million (covering all 19k genes at ~100 cells/gene/sgRNA)
```

## Factor-Based Analysis

For complex perturbation responses, decompose the per-cell perturbation effect into shared latent factors:

```python
import pertpy as pt

# FR-Perturb (Factor-Regularized Perturb-seq) decomposes perturbations into shared factors
fr = pt.tl.FRPerturb()
factors = fr.fit(adata, pert_key='perturbation', n_factors=20)
# Output: per-perturbation factor loadings; similar perturbations share factors
```

## Multiomic Perturb-seq (RNA + ATAC)

**For chromatin readout:** Use 10X Multiome with CRISPRi/a; sgRNA assignment via the same scATAC-seq library.

```python
# Multiome: scRNA + scATAC + sgRNA
# Use ArchR or Signac for ATAC integration
# Use Pertpy for RNA-side DE
import muon as mu
mdata = mu.MuData({'rna': adata_rna, 'atac': adata_atac})
# Joint differential analysis across modalities
```

## Failure Modes

### Low sgRNA detection per cell

**Trigger:** Direct-capture method on CROP-seq library, or 3'UTR barcoding on direct-capture library.
**Mechanism:** Architecture mismatch -- the sgRNA can't be detected by the wrong library prep.
**Symptom:** sgRNA assignment rate <50% of cells.
**Fix:** Match library prep to architecture; for CROP-seq, use 10X 3' chemistry; for direct-capture Perturb-seq, use the Dixit amplicon-PCR pre-sequencing.

### Mixscape filters too many cells as escapers

**Trigger:** Weak perturbation phenotype; Mixscape's NTC-subtracted signature is similar to NTC null.
**Mechanism:** Mixscape assumes a detectable signal; weak knockdown is misclassified as escaper.
**Symptom:** >50% of perturbed cells classified as "NP" (non-perturbed); known essentials show no effect.
**Fix:** Lower Mixscape stringency; skip Mixscape for low-effect perturbations; verify Cas9 expression first.

### Doublet contamination drives apparent multi-perturbation cells

**Trigger:** High cell density loading on 10X channels.
**Mechanism:** Two cells in one droplet appear to carry two sgRNAs.
**Symptom:** "Multiplet" rate >5% after sgRNA assignment.
**Fix:** Reduce cell loading per channel (5,000-7,000 instead of 10,000); run Scrublet or scDblFinder; remove doublets before sgRNA assignment.

### MAST or Wilcoxon over-call hits

**Trigger:** Using parametric DE tools on sparse, zero-inflated scRNA-seq.
**Mechanism:** These tools assume Gaussian or simpler null; single-cell data has zero-inflation that makes them over-confident.
**Symptom:** Thousands of significant DE genes per perturbation; FDR uncalibrated.
**Fix:** Use SCEPTRE (permutation-based NB GLM); Barry 2024 benchmark shows this is the only method with calibrated FDR.

### Genome-scale Perturb-seq with insufficient cells per perturbation

**Trigger:** <500 cells per perturbation in genome-scale experiment.
**Mechanism:** DE estimation requires sufficient cells per condition; <500 lacks power for moderate effects.
**Symptom:** Inconsistent hit calls across replicates; pathway analysis non-specific.
**Fix:** Scale up cell numbers; or run focused (sub-genome) Perturb-seq with more cells per pert.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| MOI for single sgRNA per cell | 0.3 | Poisson math; ~26% infected, 4% multi-infected |
| sgRNA assignment threshold | ≥10 reads of one sgRNA | Pertpy / direct-capture convention |
| Multiplet rate (post-doublet filter) | <5% | Typical 10X 3' chemistry |
| Mixscape escaper-filter (true KO retention) | 30-60% of perturbed cells | Papalexi 2021 |
| Cells per perturbation (DE power) | 500-1,000 minimum | Replogle 2022 |
| SCEPTRE permutations | 1,000+ | Barry 2024 |
| Genes per cell (QC) | ≥500-1,000 | Standard scRNA QC |
| Mt% threshold | <15-20% | Standard scRNA QC |
| Doublet detection threshold | scDblFinder, Scrublet defaults | Methods agree |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Low sgRNA detection | Architecture mismatch | Match library prep |
| Too many escapers in Mixscape | Weak phenotype | Skip Mixscape; verify Cas9 |
| Inflated DE hits | MAST / Wilcoxon used | Switch to SCEPTRE |
| Inconsistent gene effects between channels | Channel batch effect | Add channel as covariate in SCEPTRE |
| Multiplet rate >10% | Over-loading cells | Reduce loading; doublet filter |
| Per-pert DE with <100 cells | Insufficient power | Increase cell numbers; or accept low resolution |

## References

- Dixit A et al. 2016. *Cell* 167:1853. Original Perturb-seq.
- Datlinger P et al. 2017. *Nat Methods* 14:297. CROP-seq.
- Frangieh CJ et al. 2021. *Nat Genet* 53:332. Perturb-CITE-seq.
- Mimitou EP et al. 2019. *Nat Methods* 16:409. ECCITE-seq.
- Rubin AJ et al. 2019. *Cell* 176:361. Perturb-ATAC.
- Papalexi E et al. 2021. *Nat Genet* 53:322. Mixscape.
- Barry T, Mason K, Roeder K, Katsevich E. 2024. *Genome Biol* 25:124. SCEPTRE for low-MOI Perturb-seq.
- Replogle JM et al. 2022. *Cell* 185:2559. Genome-wide Perturb-seq.
- Heumos L et al. 2023. *Nat Methods* 20:1349. Pertpy framework.
- Cao J et al. 2023. *Nat Genet* 55:1894. Mixscale (dose-response Perturb-seq).

## Related Skills

- crispr-screens/library-design - Direct-capture vs CROP-seq library design
- crispr-screens/screen-qc - sgRNA assignment rates as QC
- crispr-screens/mageck-analysis - Pseudobulk analysis as alternative
- crispr-screens/hit-calling - Pseudo-bulk hit calling alternative
- single-cell/preprocessing - scRNA-seq preprocessing
- single-cell/clustering - Post-DE clustering
- single-cell/multimodal-integration - Multiome Perturb-seq
- single-cell/perturb-seq - General single-cell screen analysis
- pathway-analysis/go-enrichment - Pathway enrichment of perturbation hits
