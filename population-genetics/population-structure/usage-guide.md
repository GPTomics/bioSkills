# Population Structure - Usage Guide

## Overview

Every population-structure method returns a model-conditioned description of genotype variance, not the truth: the output is a deterministic function of the panel, the surviving SNP set, and the imposed model (continuous PCs vs K discrete clusters vs a tree-with-admixture). PCA conflates ancestry with LD blocks, inversions, relatedness, and batch, so the real work is pruning those artifacts before reading axes; ADMIXTURE Q-values are panel- and K-dependent weights on abstract clusters, not ancestry percentages, and the cross-validation minimum K is a guide, not the true number of populations; FST must be combined across SNPs as a ratio of averages (never an average of per-SNP FST), negative per-SNP values are kept not clamped, and every f3/f4/D statistic needs a block-jackknife standard error or its significance is fabricated.

## Prerequisites

- PLINK 2.0 installed: `conda install -c bioconda plink2`
- ADMIXTURE (standalone CLI): `conda install -c bioconda admixture`
- EIGENSOFT for smartpca/Tracy-Widom, FlashPCA2 for biobank-scale PCA: `conda install -c bioconda eigensoft flashpca`
- Python plotting and FST: `pip install scikit-allel numpy pandas matplotlib`; f-statistics via R `admixr` over AdmixTools
- LD-pruned, QC'd, relatedness-pruned PLINK genotypes (see plink-basics and linkage-disequilibrium)
- Conceptual prerequisites and big notes:
  - PCs are an unsupervised projection of covariance, not ancestry labels; LD pruning and outlier removal before PCA decide the axes.
  - Remove relatives before computing PCs, not after, or a relative cluster grabs a spurious component.
  - Exclude long-range-LD and inversion regions (MHC, 8p23, 17q21.31, LCT) by coordinate; they survive LD pruning and create karyotype PCs.
  - ADMIXTURE Q is a model artifact conditioned on K and panel; present a span of K and check Q stability across seeds, do not read the CV argmin as the true population count.
  - Combine per-SNP FST as a ratio of averages, prefer Hudson under sample-size asymmetry, and never clamp negative per-SNP FST.
  - Every f3/f4/D needs a block-jackknife SE; report Z with the conventional `|Z|` > 3 bar.

## Quick Start

Tell your AI agent what you want to do:
- "Run PCA on my LD-pruned genotypes and give me GWAS stratification covariates"
- "Exclude inversion regions and relatives before PCA, then plot PC1 vs PC2 by population"
- "Run ADMIXTURE for K=2 to 8 and show the CV curve as a guide, not a verdict"
- "Compute pairwise Hudson FST as a ratio of averages"
- "Test whether population C is admixed with f3 and a block jackknife"
- "Project ancient samples onto a reference PCA without shrinkage bias"

## Example Prompts

### PCA and stratification
> "LD-prune my data, exclude the MHC and 17q21.31 inversions, remove second-degree relatives, then compute 20 PCs and export them as covariates."

> "Run smartpca with Tracy-Widom significance and tell me which PCs are real versus noise."

> "Project my low-coverage ancient samples onto a modern reference PCA using least-squares projection."

### Model-based clustering
> "Run ADMIXTURE from K=2 to K=8 with cross-validation, plot the CV error, and check Q stability across three seeds before interpreting any K."

> "Make a stacked ancestry barplot for K=4 with cluster labels aligned across replicates."

### Differentiation and admixture
> "Compute Hudson FST between two populations as a ratio of averages and keep the negative per-SNP values."

> "Test whether my target population is admixed using f3, and whether there is gene flow using a D-statistic, both with block-jackknife Z scores."

## What the Agent Will Do

1. Confirm the genotypes are QC'd, LD-pruned, and relatedness-pruned, and exclude long-range-LD/inversion regions by coordinate before any PCA.
2. Compute PCs with plink2 (or smartpca for Tracy-Widom significance and projection), reporting which PCs are interpretable rather than treating every axis as ancestry.
3. Run ADMIXTURE across a span of K with cross-validation, present the CV curve as a guide, and check Q stability across seeds with label alignment before plotting.
4. Estimate FST as a ratio of averages with the Hudson estimator under sample-size asymmetry, keeping negative per-SNP values.
5. Run any f3/f4/D test through AdmixTools/admixr with a block-jackknife SE and report Z, never a bare point estimate.
6. State the panel, SNP set, model, and reference/outgroup choices so each result is interpretable as model output, not truth.

## Tips

- LD-prune and remove inversion/long-range-LD regions before PCA; a PC loading on one chromosome arm is an inversion, not a deme.
- Remove related individuals with KING (`--king-cutoff 0.0884`, second-degree) before computing PCs, then project relatives back.
- Present a span of K and check Q stability across seeds; the lowest cross-validation error is a prediction-accuracy guide, not the true number of populations.
- Combine per-SNP FST as a ratio of averages (sum numerators / sum denominators), keep negative per-SNP values, and prefer Hudson over Weir-Cockerham when sample sizes differ.
- Attach a block-jackknife standard error to every f3/f4/D; a non-negative f3 is inconclusive, and a nonzero D can be ancient structure rather than introgression.
- Use FlashPCA2 or plink2 `--pca approx` for biobank-scale data; exact PCA does not scale past ~50k samples.

## Related Skills

- plink-basics - QC, KING relatedness pruning, and fileset preparation before structure analysis
- linkage-disequilibrium - LD pruning the SNP set that PCA and ADMIXTURE require
- scikit-allel-analysis - array-scale FST and diversity windows in Python
- phasing-imputation/haplotype-phasing - phased haplotypes for haplotype-based structure methods
- comparative-genomics/introgression-detection - D/f4 introgression scans beyond pairwise tests
