# population-genetics

## Overview

Decision-grade population and statistical genetics from genotypes: allele-aware data management and QC, linkage disequilibrium, population structure, single-variant and rare-variant association, and selection scans. The category emphasizes the choices that silently change results - effect-allele bookkeeping, FST estimator and aggregation, PC-covariates vs mixed models, the mask that defines a gene-based test, and the demography confound in selection statistics.

**Tool type:** mixed | **Primary tools:** PLINK 1.9/2.0, BOLT-LMM, SAIGE, regenie, ADMIXTURE, smartpca, SKAT, scikit-allel

## Skills

| Skill | Description |
|-------|-------------|
| plink-basics | PLINK file formats and conversion, allele-order discipline (A1/A2 vs REF/ALT), and ordered sample/variant QC (missingness, MAF, HWE, sex check, KING relatedness) |
| linkage-disequilibrium | r2 vs D' (tagging vs recombination history), LD pruning vs clumping, ancestry-matched LD references, long-range-LD region handling, LD score regression |
| population-structure | PCA (smartpca/plink2) and admixture (ADMIXTURE) as model-conditioned descriptions, FST as ratio-of-averages (Weir-Cockerham vs Hudson), f-statistics with block jackknife |
| association-testing | Single-variant GWAS: linear/logistic GLM, PC-covariates vs linear mixed models (BOLT-LMM/SAIGE/regenie), SPA/Firth for case-control imbalance, LOCO, lambda vs LDSC intercept |
| rare-variant-association | Gene/region-based aggregation: burden, SKAT, SKAT-O, ACAT, STAAR via SAIGE-GENE+/regenie/SKAT R, with variant masks (functional class x MAF) as the hypothesis |
| selection-statistics | SFS-based (Tajima's D, Fay-Wu H), haplotype-based (iHS, XP-EHH, nSL, H12), and differentiation-based (FST, PBS) scans, with the demography confound and polarization requirements |
| scikit-allel-analysis | Python array-based population genetics (GenotypeArray/AlleleCountsArray), accessibility-masked diversity, ratio-of-sums FST and f-statistics, dask/zarr scaling, sgkit successor |

## Method Selection

| Task | Decision axis | Route |
|------|---------------|-------|
| Confounding control in GWAS | unrelated + only continuous ancestry -> PC-covariate GLM; related/fine-scale/cryptic structure -> LMM | association-testing |
| Case-control imbalance or low MAC | imbalance worse than ~1:10 or low minor-allele count -> SPA (SAIGE) or Firth | association-testing |
| Rare-variant signal | single carriers underpowered per-variant -> aggregate by gene/region; burden (one direction) vs SKAT (mixed) vs SKAT-O (unknown) | rare-variant-association |
| Differentiation estimate | combine per-SNP FST as ratio-of-averages; Hudson under sample-size asymmetry | population-structure, scikit-allel-analysis |
| Selection scan | incomplete/ongoing sweep -> iHS/nSL; fixed/near-fixed -> XP-EHH/Rsb; always use empirical outliers + multiple signals | selection-statistics |
| Variant set for PCA/relatedness | LD-prune and exclude long-range-LD regions by coordinate first | linkage-disequilibrium, plink-basics |

## Example Prompts

- "Convert my VCF to PLINK and run QC keeping the reference allele pinned"
- "Run QC in the correct order with controls-only HWE and KING relatedness pruning"
- "LD-prune my genotypes and clump GWAS summary statistics against an ancestry-matched reference"
- "Run PCA after removing long-range-LD regions and relatives, and test which PCs are significant"
- "Estimate FST between two populations as a ratio of averages with a block-jackknife standard error"
- "Run admixture analysis and tell me why the lowest-CV K is not the true number of populations"
- "Run a GWAS with a linear mixed model and decide whether lambda indicates confounding or polygenicity"
- "Run a gene-based rare-variant burden and SKAT-O test on my exome data"
- "Build LoF and LoF+missense masks and run SAIGE-GENE+ across genes"
- "Compute iHS and XP-EHH and explain which sweeps each detects"
- "Calculate accessibility-masked nucleotide diversity in windows with scikit-allel"

## Requirements

```bash
# PLINK 1.9 and 2.0, ADMIXTURE, EIGENSOFT (smartpca)
conda install -c bioconda plink plink2 admixture eigensoft

# Mixed-model and rare-variant association (install as needed)
conda install -c bioconda bolt-lmm saige regenie

# scikit-allel for Python population genetics; SKAT for rare-variant tests in R
pip install scikit-allel
# R: install.packages('SKAT')   # SKAT/SKAT-O/STAAR via Bioconductor/CRAN
```

## Related Skills

- **variant-calling** - VCF generation, filtering, and annotation before conversion
- **phasing-imputation** - haplotype phasing and dosage imputation feeding association and selection scans
- **causal-genomics** - fine-mapping, Mendelian randomization, and TWAS downstream of GWAS summary statistics
- **clinical-databases** - polygenic risk scores and variant prioritization consuming association results
