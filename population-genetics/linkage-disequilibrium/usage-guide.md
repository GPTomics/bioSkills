# Linkage Disequilibrium - Usage Guide

## Overview

Linkage disequilibrium is the non-random association of alleles at different loci, but the practical point is that r2 and D' answer different questions: r2 (= chi2/N) is the tagging and GWAS-power currency, while D' marks whether recombination has been observed and is upward-biased for rare variants; pruning (genotype-blind) and clumping (p-value-aware) are distinct operations, and any clumping or fine-mapping LD reference must be ancestry-matched or it fails silently.

## Prerequisites

- PLINK 1.9 and PLINK 2.0 installed: `conda install -c bioconda plink plink2 vcftools`
- scikit-allel for Python LD: `pip install scikit-allel scipy numpy`
- Input genotypes as BED/BIM/FAM or pgen (PLINK) or VCF (scikit-allel/vcftools); GWAS summary statistics for clumping
- The genome build (hg19 vs hg38) known, since long-range-LD region coordinates are build-specific
- Conceptual prerequisites and big notes:
  - Use r2 for tagging, pruning, imputation, and power; use D' only for block boundaries and recombination history.
  - D'=1 does not imply high r2 - allele-frequency asymmetry caps r2 (r2/r2max = D'^2), so a common SNP cannot tag a rare causal variant.
  - PLINK 2.0 has no bare `--r2`: choose `--r2-unphased` (composite, robust default) or `--r2-phased` (EM, assumes HWE).
  - Pruning (`--indep-pairwise`) eats genotypes and emits an independent marker set; clumping (`--clump`) eats summary statistics plus an LD reference and emits lead SNPs - they are not interchangeable.
  - Exclude long-range-LD regions (MHC, 8p23.1, 17q21.31) by coordinate before pruning or PCA; their internal r2 is high and real.
  - The clumping/fine-mapping/LDSC LD reference must match the study ancestry, ideally the study sample itself.

## Quick Start

Tell your AI agent what you want to do:
- "Prune my genotypes for PCA, excluding long-range-LD regions first"
- "Compute composite r2 between SNPs without phasing"
- "Clump my GWAS summary stats to one lead SNP per locus"
- "Find tag SNPs at r2 >= 0.8 for my variants of interest"
- "Define Gabriel haplotype blocks for a candidate region"
- "Plot LD decay with distance, stratified by population"

## Example Prompts

### LD Pruning
> "Exclude the MHC and known inversions by coordinate, then LD-prune my PLINK fileset at r2 0.1 for PCA."

> "Create an independent SNP set for ADMIXTURE with a 200kb window."

### LD Calculation
> "Compute composite (unphased) r2 between all pairs within 500kb and report D' for the HLA region."

> "Build an LD matrix for my candidate locus and tell me which SNPs tag a variant at r2 >= 0.8."

### GWAS Clumping
> "Clump my GWAS summary statistics with p1 5e-8 and r2 0.1 using an ancestry-matched LD reference."

> "Find lead SNPs for each associated locus, and flag where over-clumping might hide a secondary signal."

### Blocks and Decay
> "Identify Gabriel haplotype blocks across chromosome 6p21."

> "Plot LD decay with physical distance, stratifying by population so admixture LD does not flatten the curve."

## What the Agent Will Do

1. Identify the operation - measure (r2 vs D'), prune, clump, block, or decay - from the stated goal, not from a default.
2. For measurement, default to composite `--r2-unphased` / Rogers-Huff unless phase is trustworthy, and reserve D' for block/recombination questions.
3. For pruning, exclude long-range-LD regions by coordinate first, then slide an r2 window with a step of 1 for kb windows.
4. For clumping, override PLINK's permissive defaults (p1 5e-8, r2 0.1) and use an ancestry-matched LD reference.
5. For blocks, run the Gabriel confidence-interval method; for decay, bin r2 by distance after stratifying by population.
6. Report variant counts before/after pruning, lead SNPs per locus, and any rare-variant or wrong-ancestry caveats.

## Tips

- D'=1 is not "high LD" for prediction; a frequency mismatch makes the proxy useless even at D'=1.
- Do not interpret D' for variants below minor-allele count ~10-20; it is upward-biased and the blocks are sampling artifacts.
- Prefer composite r2 (`--r2-unphased`, Rogers-Huff) when the sample is small, structured, inbred, or has missingness.
- A flat or high-floored LD decay curve is telling about demography (bottleneck, admixture), not distance - stratify first.
- LDSC intercept-minus-1 estimates confounding but absorbs sample overlap; use the attenuation ratio.

## Related Skills

- plink-basics - format conversion and QC before any LD operation
- population-structure - PCA and ADMIXTURE on the LD-pruned marker set
- association-testing - GWAS whose summary statistics feed clumping
- selection-statistics - haplotype statistics that depend on LD structure
- causal-genomics/fine-mapping - resolving independent causal variants beyond clumping
- phasing-imputation/haplotype-phasing - phased haplotypes for EM/haplotype-based r2
