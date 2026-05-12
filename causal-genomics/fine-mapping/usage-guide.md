# Fine-Mapping - Usage Guide

## Overview

Resolve GWAS lead SNPs to credible sets of likely causal variants by fitting sparse Bayesian regressions that propagate linkage disequilibrium (LD) into posterior inclusion probabilities (PIPs). Modern fine-mapping is dominated by SuSiE (Wang 2020) and its variants: susie_rss for summary statistics, SuSiE-inf for non-sparse loci, SuSiEx for cross-ancestry, and coloc.susie for downstream colocalization. FINEMAP, CAVIAR, DAP-G, PAINTOR, PolyFun, MultiSuSiE, and FOCUS cover specialized scenarios. The hardest practical problem is LD reference mismatch; this skill bakes in the `estimate_s_rss` and `kriging_rss` diagnostics as a mandatory step.

## Prerequisites

```r
install.packages('susieR')                # CRAN; pin >= 0.12.27 for stable susie_rss API
install.packages('coloc')                 # CRAN; >= 5.2.3 for coloc.susie
install.packages(c('ggplot2', 'patchwork', 'dplyr', 'readr'))
```

```bash
# FINEMAP (CLI binary; not an R package)
# Download from http://www.christianbenner.com/

# PolyFun (Python)
git clone https://github.com/omerwe/polyfun
pip install -r polyfun/requirements.txt
# Pre-baked baseline-LF priors:
#   https://data.broadinstitute.org/alkesgroup/UKBB_LD/baselineLF2.2.UKB.tar.gz

# PAINTOR (C++ CLI)
git clone https://github.com/bogdanlab/PAINTOR_V3.0
make

# SuSiEx (C++ CLI)
git clone https://github.com/getian107/SuSiEx
make -C src

# DAP-G (CLI)
git clone https://github.com/xqwen/dap

# FOCUS (Python; for TWAS fine-mapping)
pip install pyfocus

# PLINK 1.9 / 2.0 for LD matrix generation
conda install -c bioconda plink plink2
```

## Quick Start

Tell the AI agent what to do in natural language:
- "Fine-map this GWAS locus to a 95 percent credible set using susie_rss"
- "Run the LD diagnostic estimate_s_rss before reporting credible sets"
- "Compare SuSiE and FINEMAP at this locus and reconcile disagreements"
- "Apply PolyFun functional priors to sharpen PIPs"
- "Cross-ancestry fine-map using SuSiEx with EUR, EAS, and AFR summary statistics"
- "Fine-map the HLA region with L=30 and explain why credible sets stay wide"
- "Feed susie_rss credible sets into coloc.susie for two-trait colocalization"

## Example Prompts

### Single-Locus EUR GWAS
> "Fine-map a 1 Mb window around rs12345 on chr6 using susie_rss with 1000G EUR LD. Run estimate_s_rss and report lambda. Extract 95 percent credible sets, purity, and top PIP variants."

> "Same locus, but try L=10 and L=20 and report whether the larger L changes the credible-set count."

### Polygenic / Non-Sparse Locus
> "This UK Biobank locus has 200 SNPs with -log10(p) > 4. Standard SuSiE produces 8 small credible sets that do not replicate. Refit with SuSiE-inf and compare."

### Cross-Ancestry
> "Run SuSiEx on EUR (N=500k), EAS (N=200k), and AFR (N=80k) summary statistics at chr1:50-51Mb. Compare credible set size to single-ancestry susie_rss in EUR."

### Functional Priors
> "Compute PolyFun per-SNP priors genome-wide using the baseline-LF reference, then fine-map this locus with susie_rss prior_weights set from PolyFun output. Compare to uniform-prior PIPs."

### TWAS Fine-Mapping
> "I have FUSION TWAS Z-scores for 20 co-regulated genes at chr19:45Mb. Run FOCUS to identify the likely causal gene."

### LD Diagnostic
> "Run estimate_s_rss and kriging_rss on this locus; flag any SNPs with |z_obs - z_exp| > 3 and explain whether the LD reference is suitable."

### HLA Caveat
> "I tried to fine-map a chr6:30-33 Mb autoimmune locus with susie_rss; the credible set has 60 SNPs at low purity. Explain why and recommend an HLA-specific workflow."

### Coloc Integration
> "Fine-map trait1 and trait2 separately at the same locus, then run coloc.susie. Report PP.H4 per credible-set pair."

### Reconciliation
> "SuSiE finds 3 credible sets but FINEMAP finds 1 at the same locus. Diagnose: is it convergence, LD mismatch, or non-sparse architecture?"

## What the Agent Will Do

1. **Extract locus** - 1-3 Mb window around the lead SNP from harmonized GWAS summary statistics
2. **Build / load LD matrix** - In-sample LD when available; otherwise ancestry-stratified reference panel (PLINK `--r square`)
3. **Run LD diagnostic** - `estimate_s_rss()` to compute lambda; `kriging_rss()` for per-SNP flags
4. **Choose method** - susie_rss for sparse loci, SuSiE-inf for polygenic shoulders, SuSiEx for cross-ancestry, FINEMAP for independent confirmation
5. **Fit model** - L=10 default; L=20-30 for HLA or complex loci; optional `prior_weights` from PolyFun
6. **Extract credible sets** - `fit$sets$cs`, `fit$sets$purity`, `fit$pip`; filter to purity >= 0.5
7. **Report** - Number of credible sets, size of each, purity, top PIP variant per set, and the operational caveat that the credible set is the unit of inference
8. **Optionally** - Feed into coloc.susie for colocalization; annotate variants with VEP / Ensembl

## Tips

- **In-sample LD beats reference LD** - When cohort genotypes are accessible, compute LD on the GWAS samples directly. External LD is the dominant failure mode.
- **Always run estimate_s_rss** - Lambda < 0.05 is acceptable; > 0.10 means the LD reference is wrong for this cohort.
- **Credible set is the unit of inference** - Report the set, its size, and its purity. The top PIP variant within a set is a candidate, not a conclusion.
- **L is cheap to increase** - SuSiE prunes unused effects; raise L until `sum(!fit$sets$pruned) < L`.
- **Purity filter** - Sets with `min_abs_corr < 0.5` (r2 < 0.25) are LD-confounded; drop or flag.
- **PolyFun argument** - Per-SNP causal priors go to `prior_weights`, NOT `prior_variance`. Verify with `?susie_rss`.
- **Non-sparse loci** - At biobank scale, SuSiE-inf is the default; vanilla SuSiE over-states credible-set count.
- **HLA / chr8 inversion** - Document the caveat; standard methods are unreliable. Use HLA-specific imputation for HLA.
- **PSD violations** - Add `diag(1e-4)` to the LD matrix if eigenvalues are slightly negative from finite-precision storage; or use `Matrix::nearPD`.
- **Cross-ancestry** - SuSiEx joint fine-mapping shrinks credible sets when AFR or East Asian populations contribute (shorter LD blocks resolve EUR-tagged regions).

## Related Skills

- causal-genomics/colocalization-analysis - coloc.susie operates on credible sets
- causal-genomics/mendelian-randomization - Fine-mapped cis-instruments
- causal-genomics/pleiotropy-detection - Per-credible-set pleiotropy
- population-genetics/linkage-disequilibrium - LD matrix construction
- population-genetics/association-testing - Upstream GWAS summary statistics
- workflows/gwas-pipeline - End-to-end GWAS to fine-mapping pipeline
- variant-calling/variant-annotation - Annotate credible-set variants
