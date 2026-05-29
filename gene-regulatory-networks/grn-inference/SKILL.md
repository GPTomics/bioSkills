---
name: bio-gene-regulatory-networks-grn-inference
description: Infer gene regulatory networks from bulk or general expression data with mutual-information (ARACNe) and tree-ensemble (GENIE3, GRNBoost2) methods, and infer transcription-factor protein activity from regulons with VIPER and msVIPER. Covers the activity-not-edges paradigm, the undirected-association caveat, the DREAM5 wisdom-of-crowds and method-complementarity result, AUPRC-over-AUROC evaluation, and gold-standard incompleteness. Use when inferring a regulatory network from a bulk expression matrix, finding master regulators, or scoring TF activity from a signature. For single-cell motif-pruned regulons see scenic-regulons; for co-expression modules see coexpression-networks.
tool_type: mixed
primary_tool: VIPER
---

## Version Compatibility

Reference examples tested with: VIPER 1.36+ (Bioconductor), GENIE3 1.24+ (Bioconductor), ARACNe-AP (Java, build from source), arboreto 0.1.6+ (Python GRNBoost2).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

GENIE3 expects the expression matrix as genes-in-rows, samples-in-columns (the transpose of the WGCNA convention); a transposed matrix silently produces a meaningless network.

# GRN Inference and TF Activity

**"Infer a gene regulatory network from my bulk expression data and find the master regulators"** -> Reverse-engineer TF -> target edges from an expression matrix, assemble them into regulons, then score the protein activity of each TF from a gene-expression signature.
- R: `GENIE3()` (tree-ensemble) or ARACNe-AP (mutual information) for edges
- R: `viper::aracne2regulon()` -> `msviper()` / `viper()` for TF activity

## The Single Most Important Modern Insight -- Activity, Not Edges: A Regulon Is a Multiplexed Reporter

A GRN inferred from observational expression is, by default, an **undirected statistical association graph**: correlation and mutual information are symmetric and cannot distinguish TF -> target from target -> TF or from a shared upstream driver. Tree-ensemble methods (GENIE3/GRNBoost2) appear directed only because they restrict predictors to a TF list -- that direction is an input assumption, not an inference. Individual inferred edges are therefore unreliable, and benchmarks confirm it (below).

The paradigm that survives this (the Califano-lab lineage: ARACNe -> VIPER) is to **stop trusting individual edges and instead read TF protein activity from the regulon as a whole**. VIPER treats a regulon (a TF and its inferred targets, each carrying a Mode-of-Regulation sign) as a **multiplexed reporter assay**: even if many edges are wrong, the *coordinated* up/down shift of the targets in a signature is a robust estimate of the regulator's activity (Alvarez 2016 *Nat Genet* 48:838). This is why VIPER can identify an active master regulator whose **own mRNA is unchanged** -- because the protein is regulated post-transcriptionally. Master-regulator analysis (MARINa/VIPER) is a fundamentally different computation from "the most-connected node," and edge-level precision matters less than activity inference. The Mode-of-Regulation signs are load-bearing: without them VIPER collapses to a plain enrichment test.

## Method Taxonomy

| Family | Tool | Citation | Mechanism | Structural bias |
|--------|------|----------|-----------|-----------------|
| Mutual information | ARACNe-AP | Lachmann 2016 *Bioinformatics* | MI + data-processing-inequality pruning of indirect edges | good on feed-forward loops; deletes the direct leg of true FFLs |
| Tree ensemble (RF) | GENIE3 | Huynh-Thu 2010 *PLoS ONE* | per-target random-forest variable importance | good on cascades; trades away FFLs |
| Tree ensemble (GBM) | GRNBoost2 | Moerman 2019 *Bioinformatics* | per-target gradient boosting; fast/scalable | as GENIE3; stochastic without a seed |
| Info-theoretic | CLR / MRNET | Faith 2007; Meyer 2007 | MI z-scored against per-gene background / mRMR | suppress hub artifacts |
| TF activity | VIPER / msVIPER | Alvarez 2016 *Nat Genet* | regulon-enrichment (aREA) on a signature | needs a regulon + Mode of Regulation |

DREAM5 (Marbach 2012 *Nat Methods* 9:796): no single method dominates; families make complementary errors, so an **ensemble ("wisdom of crowds") is the most robust**. And methods that excel on synthetic data collapse on real eukaryotic data (yeast was near-random) because TF and target mRNA decorrelate -- synthetic AUPRC does not transfer.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bulk RNA-seq, want a TF -> target network | GENIE3 or ARACNe-AP | tree-ensemble or MI edge inference |
| Find master regulators of a phenotype | ARACNe regulon -> msVIPER | activity inference is robust to edge errors |
| Score per-sample TF activity for stratification | `viper()` per-sample matrix | turns expression into an activity readout |
| Want robustness / no single best method | ensemble multiple inferences | DREAM5 wisdom-of-crowds |
| Single-cell data with motif resources | -> scenic-regulons | motif pruning adds directness SCENIC-style |
| Just co-expression modules (no direction) | -> coexpression-networks | WGCNA modules, no TF privileging |
| Compare TF activity between conditions | msViper on a 2-group signature | differential activity, not differential edges |

## Edge Inference with GENIE3 (R)

**Goal:** Reverse-engineer a ranked TF -> target network from a bulk expression matrix.

**Approach:** Fit a per-target random forest predicting each gene from candidate regulators (TFs); the regulator's variable importance is the edge weight. Restrict predictors to a TF list to orient edges.

```r
library(GENIE3)

# GENIE3 convention: genes in ROWS, samples in COLUMNS (transpose of WGCNA).
expr <- as.matrix(read.csv('normalized_counts.csv', row.names = 1))
regulators <- readLines('tf_list.txt')                  # candidate TFs only

set.seed(42)                                            # tree ensembles are stochastic
weight_matrix <- GENIE3(expr, regulators = regulators, treeMethod = 'RF',
                        K = 'sqrt', nTrees = 1000, nCores = 8)
link_list <- getLinkList(weight_matrix)                 # ranked edge list (NOT thresholded)
head(link_list)
```

## Edge Inference with ARACNe-AP (Java CLI)

**Goal:** Build a mutual-information network with indirect edges pruned.

**Approach:** ARACNe-AP is a two-phase Java pipeline: compute the MI threshold, run many bootstrap reconstructions, then consolidate them (with data-processing-inequality pruning) into a final network. Running a single bootstrap or skipping consolidation is the classic misuse.

```bash
# Phase 1: MI threshold at a chosen p-value (needs the expression matrix + TF list).
java -Xmx32G -jar aracne.jar -e expr.txt -o out/ --tfs tf_list.txt \
    --pvalue 1E-8 --seed 1 --calculateThreshold

# Phase 2: many bootstraps (vary --seed) -- 100 is conventional.
for s in $(seq 1 100); do
    java -Xmx32G -jar aracne.jar -e expr.txt -o out/ --tfs tf_list.txt \
        --pvalue 1E-8 --seed $s
done

# Phase 3: consolidate bootstraps into the final network (DPI + a Poisson edge-significance
# test with Bonferroni correction across bootstraps).
java -Xmx32G -jar aracne.jar -o out/ --consolidate
```

## TF Activity with VIPER / msVIPER (R)

**Goal:** Infer transcription-factor protein activity from a regulon and an expression signature, and rank master regulators.

**Approach:** Convert an ARACNe network into a regulon object (assigning each target a Mode-of-Regulation sign and likelihood), build a null model by sample permutation, then run msVIPER on a two-group signature (master regulators) or VIPER per sample (activity matrix).

```r
library(viper)

# Build the regulon from the ARACNe network + matched expression (assigns Mode of Regulation).
# ARACNe-AP network.txt has a header + 4 columns (Regulator, Target, MI, p-value); viper's
# '3col' reader wants Regulator/Target/MI with no header, so strip them first (shell):
#   tail -n +2 out/network.txt | cut -f1-3 > net_3col.txt
# ('adj' is the legacy ARACNE adjacency-matrix format, not ARACNe-AP.)
regulon <- aracne2regulon('net_3col.txt', eset, format = '3col')

# msVIPER: master regulators of a two-group contrast.
signature <- rowTtest(eset, pheno = 'group', group1 = 'tumor', group2 = 'normal')
sig_z <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
nullmodel <- ttestNull(eset, pheno = 'group', group1 = 'tumor', group2 = 'normal', per = 1000)
mra <- msviper(sig_z, regulon, nullmodel)
summary(mra)                                            # top master regulators by NES

# VIPER: a per-sample TF-activity matrix for clustering/stratification.
activity <- viper(eset, regulon, method = 'scale')
```

For single cells or tissues lacking a matched network, `metaVIPER` integrates multiple interactomes; DIGGIT then intersects master regulators with genetic alterations to nominate causal drivers.

## Per-Method Failure Modes

### Fabricated directionality
**Trigger:** presenting an MI/correlation network as a directed causal GRN. **Mechanism:** symmetric measures carry no direction; the TF-list restriction is an assumption. **Symptom:** arrowheads with no perturbation/time/sequence support. **Fix:** state edges are associations; reserve causal claims for perturbation-validated edges.

### Master regulators from edge counts
**Trigger:** calling the most-connected node a master regulator. **Mechanism:** MRA (VIPER) is regulon enrichment in a signature, not node degree. **Symptom:** "hub = driver" with no activity computation. **Fix:** run msVIPER; report NES.

### Skipping ARACNe consolidation
**Trigger:** a single bootstrap, or no `--consolidate`. **Mechanism:** the network is unstabilized and DPI/Bonferroni unapplied. **Symptom:** noisy, non-reproducible edges. **Fix:** run ~100 bootstraps then consolidate.

### Missing Mode of Regulation in VIPER
**Trigger:** a regulon without target signs. **Mechanism:** aREA needs activating/repressing signs so a repressed-target down-shift counts toward activation. **Symptom:** VIPER behaves like a plain enrichment test. **Fix:** build the regulon with `aracne2regulon` (which assigns MoR).

### AUROC-only / synthetic-only validation
**Trigger:** reporting AUROC near 1, or validating only on simulated data. **Mechanism:** with ~0.1-1% true edges AUROC hides near-random AUPRC; synthetic success does not transfer (DREAM5). **Symptom:** no AUPRC, no real-data gold standard. **Fix:** report AUPRC + early precision against an independent gold standard; acknowledge gold-standard incompleteness.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| ARACNe bootstraps ~100 then consolidate | ARACNe-AP workflow | stabilizes edges; DPI + Poisson edge test, Bonferroni-corrected |
| ARACNe MI p-value 1E-8 | ARACNe-AP default-scale | controls edge false positives genome-wide |
| GENIE3 nTrees = 1000, K = 'sqrt' | GENIE3 defaults | variance/runtime trade-off for importances |
| VIPER/msVIPER null permutations ~1000 | VIPER convention | calibrates the NES null distribution |
| Report AUPRC + early precision (not AUROC) | Marbach 2012 / Pratapa 2020 | AUROC misleads under sparse positives |
| Set a seed for GENIE3/GRNBoost2 | reproducibility | tree ensembles are stochastic |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| meaningless GENIE3 network | matrix transposed (samples in rows) | genes in rows, samples in columns |
| ARACNe network unstable across runs | single bootstrap / no consolidate | run ~100 bootstraps then `--consolidate` |
| VIPER acts like plain enrichment | regulon lacks Mode of Regulation | build via `aracne2regulon` |
| top regulator is just highly expressed | using degree/expression as "activity" | use msVIPER NES |
| great synthetic accuracy, fails on real data | over-fit to in-silico benchmark (DREAM5) | validate on real gold standards; report AUPRC |

## References

- Margolin AA, et al. 2006. ARACNE: reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics* 7(Suppl 1):S7.
- Lachmann A, Giorgi FM, Lopez G, Califano A. 2016. ARACNe-AP. *Bioinformatics* 32(14):2233-2235.
- Huynh-Thu VA, et al. 2010. Inferring regulatory networks using tree-based methods (GENIE3). *PLoS ONE* 5(9):e12776.
- Moerman T, et al. 2019. GRNBoost2 and Arboreto. *Bioinformatics* 35(12):2159-2161.
- Faith JJ, et al. 2007. Large-scale mapping and validation of E. coli transcriptional regulation (CLR). *PLoS Biol* 5(1):e8.
- Alvarez MJ, et al. 2016. Network-based inference of protein activity (VIPER). *Nat Genet* 48(8):838-847.
- Marbach D, et al. 2012. Wisdom of crowds for robust gene network inference (DREAM5). *Nat Methods* 9(8):796-804.
- Pratapa A, et al. 2020. Benchmarking single-cell GRN inference (BEELINE). *Nat Methods* 17(2):147-154.

## Related Skills

- scenic-regulons - single-cell regulons with motif-pruning directness
- coexpression-networks - undirected co-expression modules (no TF privileging)
- differential-networks - VIPER differential activity / rewiring between conditions
- multiomics-grn - enhancer-driven directed GRNs from accessibility
- differential-expression/de-results - signatures that feed msVIPER
- single-cell/perturb-seq - interventional validation of inferred regulation
