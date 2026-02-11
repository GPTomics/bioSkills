# Conservation Genetics Usage Guide

## Overview

Assesses the genetic health of populations for conservation management. Covers F-statistics and genetic diversity metrics (hierfstat), allelic richness, pairwise population differentiation, runs of homozygosity for inbreeding detection (detectRUNS), and effective population size estimation across multiple time scales: contemporary Ne (NeEstimator LD method), recent Ne trajectory over ~200 generations (GONE2), and deep demographic history over thousands of generations (Stairway Plot 2, PSMC). Designed for microsatellite, RADseq, or WGS data from threatened or managed species.

## Prerequisites

- R with hierfstat (`install.packages('hierfstat')`)
- adegenet (`install.packages('adegenet')`)
- detectRUNS (`install.packages('detectRUNS')`)
- GONE2 (`remotes::install_github('esrud/GONE2')`)
- poppr for private alleles (`install.packages('poppr')`)
- NeEstimator v2 (CLI: https://github.com/bunop/NeEstimator2.X)
- Stairway Plot 2 (Java: https://github.com/xiaoming-liu/stairway-plot-v2)
- PSMC (CLI: https://github.com/lh3/psmc)

## Quick Start

Tell your AI agent what you want to do:

- "Calculate genetic diversity and Fst for my endangered species populations"
- "Estimate effective population size from my RADseq data"
- "Detect runs of homozygosity and inbreeding in my conservation dataset"
- "Reconstruct demographic history for my species using PSMC"
- "Compare allelic richness across populations with different sample sizes"

## Example Prompts

### Genetic Diversity

> "I have genepop files from 5 populations of an endangered bird. Calculate Ho, He, Fis, pairwise Fst with bootstrap confidence intervals, and rarefied allelic richness."

> "Compare genetic diversity metrics between island and mainland populations of my study species. Test whether island populations have significantly lower heterozygosity."

### Inbreeding Detection

> "I have PLINK ped/map files from a captive breeding program. Detect runs of homozygosity, calculate F_ROH for each individual, and classify ROH by length classes."

> "Identify the most inbred individuals in my dataset for breeding management recommendations."

### Effective Population Size

> "Estimate the recent Ne trajectory for my species using GONE2 from PLINK bed/bim/fam files. Has the population been declining?"

> "Run PSMC on my whole-genome BAM file to reconstruct demographic history over the last million years. Use a mutation rate of 1.4e-8 and generation time of 5 years."

> "Estimate contemporary Ne using the LD method in NeEstimator. Is it above the critical threshold of 50?"

### Bottleneck Detection

> "Test whether my population shows signatures of a recent bottleneck using heterozygosity excess."

## What the Agent Will Do

1. Load genotype data and convert to the appropriate format (genind, hierfstat, PLINK)
2. Calculate per-population diversity metrics: Ho, He, Fis, allelic richness, private alleles
3. Compute pairwise Fst with bootstrap confidence intervals
4. Detect runs of homozygosity and calculate F_ROH inbreeding coefficients
5. Classify ROH into length classes corresponding to the timing of inbreeding events
6. Estimate effective population size using the appropriate method for the data and time scale
7. Generate diagnostic plots: Ne trajectories, ROH distributions, diversity comparisons
8. Interpret results against conservation thresholds (Ne < 50 critical, Ne < 500 vulnerable)

## Tips

- The 50/500 rule is a guideline: Ne > 50 to avoid inbreeding depression, Ne > 500 for adaptive potential
- GONE2 requires at least 10,000 SNPs and 50 diploid individuals for reliable recent Ne estimates
- PSMC requires a single high-coverage (>20x) whole genome; not suitable for RADseq data
- Stairway Plot 2 requires only the site frequency spectrum and works with RADseq or WGS
- F_ROH is more informative than Fis for inbreeding detection because it captures genomic-level homozygosity
- ROH > 16 Mb indicate very recent inbreeding (parents or grandparents); 1-4 Mb indicates historical bottlenecks
- Rarefied allelic richness is essential for comparing populations with different sample sizes
- Ne/N ratio is typically 0.1-0.3; very low ratios suggest high variance in reproductive success
- For NeEstimator LD method, use pcrit = 0.02 to exclude rare alleles that inflate Ne estimates
- Always pair Ne estimates with census size (N) and trend data for conservation assessments

## Related Skills

- landscape-genomics - Adaptive variation and genotype-environment associations
- species-delimitation - Taxonomic unit definition for conservation
- population-genetics/population-structure - Population stratification
- population-genetics/selection-statistics - Selection signatures
- variant-calling/vcf-basics - VCF input from RADseq/WGS
