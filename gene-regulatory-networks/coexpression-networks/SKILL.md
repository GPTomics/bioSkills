---
name: bio-gene-regulatory-networks-coexpression-networks
description: Build weighted gene co-expression networks to identify modules of co-regulated genes, relate them to phenotypes, and find hub genes using WGCNA, hdWGCNA, MEGENA, CEMiTool, and Gaussian graphical models. Covers signed-network choice, soft-threshold selection, module preservation, and the marginal-vs-partial-correlation distinction. Use when finding co-expression modules, identifying hub genes, relating gene networks to clinical or experimental traits, or building single-cell co-expression networks. For directed TF-target inference see scenic-regulons and grn-inference; for condition rewiring see differential-networks.
tool_type: r
primary_tool: WGCNA
---

## Version Compatibility

Reference examples tested with: WGCNA 1.72+, hdWGCNA 0.3+, CEMiTool 1.26+, GENIE3-adjacent GGM via GeneNet 1.2.16+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

WGCNA argument defaults differ by entry point: `pickSoftThreshold()` and `blockwiseModules()` default to `networkType='unsigned'`. The signed-network choice (below) must be set identically at every step or the soft power and modules silently mismatch.

# Co-expression Networks

**"Find co-expression modules and hub genes from my expression data"** -> Build a weighted gene co-expression network, detect modules of co-regulated genes by clustering a topological-overlap dissimilarity, summarize each module by its eigengene, and relate modules to sample traits.
- R: `WGCNA::blockwiseModules()` for network construction + module detection (bulk)
- R: `hdWGCNA` metacell workflow for single-cell expression
- R: `GeneNet`/graphical lasso when direct (not indirect) edges are required

## The Single Most Important Modern Insight -- Co-expression Is Marginal Correlation; Regulation Is Conditional Dependence

A WGCNA edge between genes B and C exists whenever a common driver A correlates with both -- even if B and C never interact. If A regulates B and C, then B and C are conditionally independent given A, yet a co-expression network still draws a strong B-C edge. **A co-expression module is a descriptive object ("these genes vary together"), not a regulatory network, even when its hub is a transcription factor.** The dividing line is marginal vs partial correlation: WGCNA, MEGENA, and CEMiTool measure **marginal** correlation (direct + indirect edges mixed together), while Gaussian graphical models (GeneNet, graphical lasso) estimate **partial** correlation (conditional dependence, i.e. direct edges only). Decide which object the biological question needs before choosing a tool. Causal/regulatory language ("module X is driven by hub TF Y") from a marginal co-expression edge is the single most common over-claim in this domain.

A secondary, load-bearing caveat: the **scale-free topology assumption underpinning soft-threshold selection is empirically weak**. Broido & Clauset 2019 (*Nat Commun* 10:1017) found scale-free structure in only ~4% of ~1000 real networks; Khanin & Wit 2006 (*J Comput Biol* 13:810) found none of 10 biological networks fit a pure power law. The scale-free R^2 >= 0.8 rule is a heuristic for picking the power where the connectivity curve flattens, **not** a hypothesis test confirming the biology is scale-free -- so it does not license "the network is scale-free, therefore hubs are master regulators."

## Co-expression Method Taxonomy

| Method | Edge type | Strength | Use / fails when |
|--------|-----------|----------|------------------|
| WGCNA (Langfelder & Horvath 2008) | marginal, soft-threshold + TOM | de facto standard; eigengenes; preservation | bulk, n>=15-20; fails on raw scRNA-seq |
| hdWGCNA (Morabito 2023) | marginal, on metacells | flattens dropout/zero-inflation | single-cell; fails if metacells overlap (pseudo-replication) |
| MEGENA (Song & Zhang 2015) | marginal, planar filtered network | principled sparsification; multiscale/nested modules | hierarchical biology; nested output breaks WGCNA-style one-module-per-gene code |
| CEMiTool (Russo 2018) | marginal, auto-beta | fast first pass + GSEA/ORA + HTML report | exploratory; S4 object (accessors, not `$`); silently variance-filters genes |
| GeneNet / graphical lasso | **partial** (direct edges) | distinguishes direct from indirect | when "who regulates whom directly" matters; collapses at small n / large p |

Rule of thumb: module discovery and trait association -> WGCNA (signed); single-cell -> hdWGCNA; "is this edge direct or indirect?" -> a Gaussian graphical model.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bulk RNA-seq, >=15-20 samples, want modules + trait links | WGCNA signed, bicor | stable marginal modules; eigengene-trait correlation |
| Single-cell / sparse counts | hdWGCNA on metacells | naive WGCNA fails on zero-inflation; metacells flatten dropout |
| Need direct vs indirect edges | GeneNet / graphical lasso | partial correlation removes confounded indirect edges |
| Hierarchical / multiscale structure expected | MEGENA | nested modules at multiple resolutions |
| Want directed TF -> target regulons | -> scenic-regulons (single-cell) / grn-inference (bulk) | co-expression is undirected; motif/regression priors add direction |
| Compare networks across conditions | -> differential-networks | rewiring is a different question from module discovery |
| <15 samples | none reliable | correlation estimates too noisy; report this, do not force WGCNA |

## WGCNA: Build a Signed Network

**Goal:** Detect robust co-expression modules from a bulk expression matrix and choose network parameters defensibly.

**Approach:** Use a **signed** network with **biweight midcorrelation (bicor)** so activators and repressors are not merged into one module and outliers are down-weighted; pick the soft power on the SAME network type used for construction; set `maxBlockSize` above the gene count to avoid the block artifact.

```r
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# WGCNA convention: genes as columns, samples as rows
expr <- t(as.matrix(read.csv('normalized_counts.csv', row.names = 1)))
gsg <- goodSamplesGenes(expr, verbose = 0)
expr <- expr[gsg$goodSamples, gsg$goodGenes]

# Soft power on the SIGNED fit (must match construction below)
powers <- 1:20
sft <- pickSoftThreshold(expr, powerVector = powers, networkType = 'signed', verbose = 0)
soft_power <- sft$powerEstimate            # signed networks usually land near 12

# Signed network, bicor with maxPOutliers to avoid spurious outlier flags at modest n.
# maxBlockSize >= n_genes keeps everything in one block (no cross-block blindness).
net <- blockwiseModules(
    expr, power = soft_power,
    networkType = 'signed', TOMType = 'signed',
    corType = 'bicor', maxPOutliers = 0.05,
    minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2,
    maxBlockSize = ncol(expr) + 1,
    numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 0
)
module_colors <- labels2colors(net$colors)
table(module_colors)                       # large 'grey' fraction = poor fit, not a module
```

## WGCNA: Eigengenes, Hubs, and Trait Relationships

**Goal:** Relate modules to sample traits and identify defensible hub genes.

**Approach:** Summarize each module by its eigengene (first PC), correlate eigengenes with traits, and define hubs by **module membership kME** (signed, bounded, comparable across modules) rather than raw connectivity.

```r
# Recompute eigengenes from the COLOR labels so ME/kME columns are MEturquoise/kMEturquoise.
# (net$MEs uses numeric names ME0/ME1... under numericLabels=TRUE and won't match colors.)
MEs <- orderMEs(moduleEigengenes(expr, module_colors)$eigengenes)
traits <- read.csv('sample_traits.csv', row.names = 1)

module_trait_cor  <- cor(MEs, traits, use = 'p')
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(expr))

# Hub = high module membership (kME), the WGCNA-preferred definition.
# kME is signed and bounded [-1,1] with a p-value; prefer it over kWithin connectivity.
kME <- signedKME(expr, MEs)
moi <- 'turquoise'
moi_genes <- colnames(expr)[module_colors == moi]
hubs <- sort(kME[moi_genes, paste0('kME', moi)], decreasing = TRUE)
head(hubs, 20)
```

## WGCNA: Module Preservation (the step that makes a module real)

**Goal:** Test whether discovered modules reproduce in an independent dataset -- the scientific claim, not module detection itself.

**Approach:** Run `modulePreservation()` with the discovery network as reference and an independent cohort as test; read Zsummary and medianRank.

```r
# multiData/multiColor hold reference + test expression and the reference module labels
mp <- modulePreservation(
    multiData, multiColor,
    referenceNetworks = 1, nPermutations = 200,
    randomSeed = 1, verbose = 0
)
stats <- mp$preservation$Z$ref.reference$inColumnsAlsoPresentIn.test
stats[, c('moduleSize', 'Zsummary.pres')]
# Zsummary > 10 strong preservation; 2-10 weak/moderate; < 2 no evidence (Langfelder 2011).
# medianRank (mp$preservation$observed) compares modules to each other, size-independent.
```

## hdWGCNA: Single-Cell Co-expression

**Goal:** Build co-expression networks from scRNA-seq without dropout-driven spurious correlation.

**Approach:** Aggregate transcriptionally similar cells into metacells (pseudobulk) to flatten zero-inflation; cap metacell sharing to avoid pseudo-replication; then run standard WGCNA on the metacells.

```r
library(hdWGCNA); library(Seurat)
seurat_obj <- readRDS('clustered.rds')
seurat_obj <- SetupForWGCNA(seurat_obj, gene_select = 'fraction', fraction = 0.05,
                            wgcna_name = 'hdwgcna')

# Build metacells WITHIN cell types; max_shared caps how many cells two metacells share.
# Heavy overlap inflates downstream correlations (pseudo-replication) -- keep it low.
seurat_obj <- MetacellsByGroups(seurat_obj, group.by = 'cell_type',
                                k = 25, max_shared = 10, ident.group = 'cell_type')
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(seurat_obj, group.by = 'cell_type', group_name = 'all')
seurat_obj <- TestSoftPowers(seurat_obj, networkType = 'signed')
seurat_obj <- ConstructNetwork(seurat_obj, setDatExpr = FALSE, tom_name = 'hdwgcna')
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)   # kME per module
```

## Direct vs Indirect Edges: Gaussian Graphical Model

**Goal:** Recover edges that survive conditioning on all other genes (direct dependence), not marginal co-expression.

**Approach:** Estimate a shrinkage partial-correlation matrix (n << p safe) and test edges by local FDR.

```r
library(GeneNet)
# expr: samples x genes
pcor   <- ggm.estimate.pcor(expr)              # shrinkage partial correlations
edges  <- network.test.edges(pcor, direct = FALSE, plot = FALSE)
net_gg <- extract.network(edges, method.ggm = 'prob', cutoff.ggm = 0.9)
# Far sparser than WGCNA -- that is the point: indirect edges have been removed.
```

CEMiTool (`cemitool(expr, annot)`) returns an **S4 object** -- use accessors (`module_genes()`, `get_hubs()`, `generate_report()`), never `$` -- and silently variance-filters genes with `filter=TRUE`. MEGENA returns nested modules; do not feed them to one-module-per-gene WGCNA code.

## Per-Method Failure Modes

### Unsigned-by-default network
**Trigger:** leaving `networkType='unsigned'` (the default) then interpreting modules as co-regulated programs. **Mechanism:** unsigned uses |cor|, so a gene and its strong anti-correlate land in the same module. **Symptom:** modules mixing clearly opposing programs (e.g. proliferation and quiescence). **Fix:** use `networkType='signed'` (and `TOMType='signed'`) at pickSoftThreshold AND blockwiseModules.

### Soft power chosen on the wrong network type
**Trigger:** `pickSoftThreshold()` run unsigned, then `blockwiseModules(networkType='signed')`. **Mechanism:** the scale-free fit and chosen power differ by network type. **Symptom:** near-disconnected network or one giant module. **Fix:** set `networkType` identically in both calls (signed power is roughly double unsigned).

### blockwiseModules block artifact
**Trigger:** n_genes > `maxBlockSize` (historical default 5000) with no mention of blocking. **Mechanism:** TOM is computed within each block only; cross-block co-expression is invisible and assignments shift with block size/RAM. **Symptom:** module membership changes when rerun on a different machine. **Fix:** set `maxBlockSize >= n_genes`, or report the block size.

### Naive WGCNA on raw scRNA-seq
**Trigger:** running WGCNA on sparse single-cell counts. **Mechanism:** zero-inflation creates a spike of zero correlations and dropout-driven spurious correlation. **Symptom:** a pathological correlation distribution and unstable modules. **Fix:** use hdWGCNA metacells (and cap `max_shared`).

### Detection without preservation
**Trigger:** reporting modules from one cohort with no replication test. **Mechanism:** any dendrogram can be cut into modules. **Symptom:** modules that do not reproduce in independent data. **Fix:** run `modulePreservation()`; report Zsummary/medianRank.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Samples >= 15 (ideally >= 20) | Horvath/Langfelder WGCNA FAQ | correlation estimates too noisy below 15; modules become random |
| Scale-free R^2 >= 0.8 (0.8-0.9) | WGCNA convention | heuristic for power where connectivity flattens -- NOT proof of scale-free biology |
| Signed soft power ~12 (vs ~6 unsigned) | WGCNA FAQ table | signed adjacency = ((1+cor)/2)^power needs higher power for the same connectivity |
| maxPOutliers = 0.05-0.10 | Langfelder & Horvath 2012 | prevents bicor flagging legitimate observations as outliers at modest n |
| minModuleSize = 30 (default); mergeCutHeight = 0.25 (chosen; default is 0.15) | WGCNA | smaller modules are usually noise; 0.25 merges eigengenes correlating > 0.75 |
| Zsummary > 10 / 2-10 / < 2 | Langfelder 2011 | strong / weak-moderate / no preservation evidence |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| One giant module + large grey | power too low or wrong network type | re-pick power on the correct `networkType` |
| `cemitool` results via `$` return NULL | CEMiTool returns an S4 object | use accessors `module_genes()`, `get_hubs()` |
| hub list dominated by housekeeping genes | hubness tracking mean expression | define hubs by kME, control for expression level |
| top module is the sequencing batch | batch correlates with many genes | correct batch before network construction (-> differential-expression/batch-correction) |
| modules irreproducible across reruns | block artifact | set `maxBlockSize >= n_genes` |

## References

- Langfelder P, Horvath S. 2008. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9:559.
- Langfelder P, Horvath S. 2012. Fast R functions for robust correlations and hierarchical clustering. *J Stat Softw* 46(11):1-17. -- bicor / maxPOutliers.
- Langfelder P, Luo R, Oldham MC, Horvath S. 2011. Is my network module preserved and reproducible? *PLoS Comput Biol* 7(1):e1001057.
- Broido AD, Clauset A. 2019. Scale-free networks are rare. *Nat Commun* 10:1017. -- soft-threshold caveat.
- Khanin R, Wit E. 2006. How scale-free are biological networks? *J Comput Biol* 13(3):810-818.
- Song WM, Zhang B. 2015. Multiscale embedded gene co-expression network analysis (MEGENA). *PLoS Comput Biol* 11(11):e1004574.
- Russo PST, et al. 2018. CEMiTool. *BMC Bioinformatics* 19:56.
- Morabito S, Reese F, Rahimzadeh N, Miyoshi E, Swarup V. 2023. hdWGCNA. *Cell Rep Methods* 3(6):100498.

## Related Skills

- differential-networks - compare co-expression structure between conditions (rewiring)
- scenic-regulons - directed TF regulons from single-cell data (motif-pruned)
- grn-inference - bulk directed GRN inference and TF protein-activity (VIPER)
- differential-expression/batch-correction - remove batch effects before network construction
- single-cell/preprocessing - QC and normalization for single-cell inputs to hdWGCNA
- pathway-analysis/go-enrichment - functional enrichment of co-expression modules
