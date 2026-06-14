# PLINK Basics - Usage Guide

## Overview

PLINK is the workhorse for genotype data management and quality control, but it is a stateful rewriter of allele bookkeeping rather than a calculator: PLINK 1.x re-derives the A1 (effect) allele as the minor allele on every load, so effect estimates and PRS weights silently lose their allele anchor unless the alleles are pinned, while PLINK 2.0 tracks genuine REF/ALT. The QC decisions that matter most are not the thresholds but the order (variant missingness before sample missingness), the controls-only HWE convention (automatic in 1.9, manual in 2.0), and the case/control differential-missingness confounder that injects false association.

## Prerequisites

- PLINK 1.9 and PLINK 2.0 installed: `conda install -c bioconda plink plink2`
- Input genotypes as VCF/BCF, BED/BIM/FAM, PED/MAP, or pgen/pvar/psam
- A genome build (hg19 vs hg38) known for the data - PAR boundaries and liftover depend on it
- Conceptual prerequisites and big notes:
  - A1/A2 (1.x) is the counted/effect-allele pair and carries no reference or strand meaning; REF/ALT (2.0) is genuine. Harmonize alleles to a fixed reference before combining cohorts.
  - Run `--geno` (variant) and `--mind` (sample) in separate invocations, variant first, to control the implicit order.
  - HWE filtering is controls-only and should use `midp`; in plink2 the controls-only behavior must be added explicitly.
  - Relatedness in a structured sample needs KING (`--make-king`), not PI_HAT (`--genome`).
  - Split the pseudoautosomal region (`--split-par`) before any sex check or X-specific analysis.

## Quick Start

Tell your AI agent what you want to do:
- "Convert this VCF to a PLINK fileset, keeping REF/ALT honest"
- "Run standard GWAS QC on my genotypes in the correct order"
- "Apply controls-only HWE filtering with mid-p"
- "Find related individuals with KING and prune to unrelated"
- "Check reported sex against X heterozygosity after splitting PAR"
- "Test for differential missingness between cases and controls"
- "Merge two genotype datasets and resolve strand flips"

## Example Prompts

### Conversion
> "Convert imputed.vcf.gz to a PLINK 2.0 pgen fileset preserving dosages, then export a biallelic BED with the reference allele pinned from the GRCh38 FASTA."

### Quality control
> "Run QC on cohort.bed: variant missingness at 0.02, then sample missingness at 0.02, MAF 0.01, and controls-only HWE at 1e-6 with mid-p."

### Sample-level QC
> "Compute KING kinship and remove one of each pair closer than second-degree, and flag heterozygosity outliers beyond 3 SD on LD-pruned SNPs."

### Confounder checks
> "Run a differential-missingness test between my cases and controls and give me the list of variants to drop before GWAS."

## What the Agent Will Do

1. Identify the input format and target, and choose PLINK 1.9 vs 2.0 from the operation (PED/MAP and `--genome` force 1.9; dosages and KING favor 2.0).
2. Convert format, preserving REF/ALT where possible and applying `--keep-allele-order` on any BED export.
3. Apply QC in the order that avoids dropping good samples: variant missingness, then sample missingness, MAF, and controls-only HWE with mid-p.
4. Run sample-level QC: sex check after PAR split, heterozygosity outliers on LD-pruned SNPs, and KING relatedness pruning.
5. Check the differential-missingness confounder before any association.
6. Report which samples and variants were removed at each step and why, so the filtering is auditable.

## Tips

- Stay in `.pgen` through QC if every downstream tool supports it; export `.bed` only as a last step and always with `--keep-allele-order`.
- The thresholds in the skill are starting conventions. Plot missingness, MAF, het-F, and KING distributions and pick cutoffs at the natural gaps.
- Sex-check and heterozygosity failures usually mean a sample swap or contamination - investigate the sample rather than silently dropping it.
- Palindromic A/T and C/G SNPs cannot be strand-resolved by `--flip`; resolve by allele frequency or drop them before merging cohorts.
- After merging, re-apply allele-order discipline - a merge can re-derive A1/A2.

## Related Skills

- linkage-disequilibrium - LD pruning and clumping on QC'd genotypes
- population-structure - PCA and ADMIXTURE after QC and relatedness pruning
- association-testing - GWAS with `--glm` on the filtered fileset
- variant-calling/vcf-basics - VCF generation and manipulation before conversion
- phasing-imputation/genotype-imputation - imputed dosages that enter as `.pgen`
