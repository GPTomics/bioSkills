# Allele-Specific Accessibility - Usage Guide

## Overview

Detect cis-regulatory genetic effects on chromatin accessibility from heterozygous SNPs. Covers WASP reference-bias correction (mandatory), GATK ASEReadCounter for allele counting, RASQUAL for joint cohort-level + allelic imbalance caQTL mapping, QuASAR for genotype-free ASE, and per-peak aggregation strategies. Builds chromatin QTL (caQTL) maps for fine-mapping and validates GWAS variant function with allelic imbalance evidence.

## Prerequisites

```bash
# WASP for mapping bias correction (mandatory)
git clone https://github.com/bmvdgeijn/WASP

# GATK 4 for allele counting
conda install -c bioconda gatk4

# RASQUAL for joint caQTL
git clone https://github.com/natsuhiko/rasqual

# Phasing (input to WASP)
conda install -c bioconda shapeit5 beagle whatshap

# Other tools
conda install -c bioconda samtools bcftools vcftools plink2 bowtie2 bwa-mem2
```

```r
install.packages('MatrixEQTL')                # CRAN
remotes::install_github('piquelab/QuASAR')    # GitHub-only (not on Bioconductor/CRAN)
```

Inputs:
- ATAC-seq BAM (deduplicated, MAPQ-filtered)
- Genotype VCF for the same individual (filtered to heterozygous SNPs)
- Reference genome FASTA
- Phased haplotypes (SHAPEIT5/BEAGLE/whatshap output)
- Consensus peakset (atac-seq/consensus-peakset)

## Quick Start

Tell your AI agent what you want to do:
- "Run WASP filter to correct reference-allele mapping bias before any ASE analysis"
- "Use GATK ASEReadCounter on WASP-filtered BAM to count REF and ALT reads at heterozygous SNPs"
- "Aggregate ASE within peaks: per-peak binomial test with effect size >= 0.2 and FDR < 0.05"
- "Run RASQUAL joint cis-caQTL combining total counts and allelic counts on a cohort of N=50"
- "Phase genotypes with SHAPEIT5 before WASP since unphased sites can't assign reads to haplotypes"
- "Cross-validate predicted variant effects (chromBPNet) against observed allelic imbalance"

## Example Prompts

### Mandatory WASP Filtering
> "Run WASP on a sample's ATAC BAM: find_intersecting_snps, re-align flipped reads with bowtie2, filter_remapped_reads to keep only consistently-mapped, and merge with non-overlapping reads. The output BAM is safe for ASE counting."

### GATK ASE Counting
> "On a WASP-filtered BAM with heterozygous SNPs in VCF, run GATK ASEReadCounter to produce per-SNP REF/ALT counts. Filter to SNPs with totalCount >= 30 for adequate power."

### Within-Peak ASE
> "Aggregate per-SNP ASE counts within each peak in the consensus peakset; pooled binomial test per peak with effect size |ref_frac - 0.5| >= 0.2 and BH-adjusted p < 0.05."

### Cohort caQTL with RASQUAL
> "Run RASQUAL per feature: stream the cis-window genotype VCF (with allele-specific counts) through tabix into rasqual, which models genotype/allelic correlation internally -- no external LD-precompute step. Joint p-value typically gives 1.5-3x more caQTLs than MatrixEQTL alone at N=50."

### Genotype-Free ASE
> "I don't have genotypes but want to detect ASE. Use QuASAR which infers genotypes from the BAM directly using a HMM."

### Variant Effect Cross-Validation
> "Compare chromBPNet-predicted variant effects (atac-seq/deep-learning-atac) against observed allelic imbalance at the same SNPs; concordance > 70% supports the deep-learning model."

## What the Agent Will Do

1. Verify single-sample setup (BAM + VCF + reference + phased haplotypes)
2. Phase genotypes (SHAPEIT5/BEAGLE for cohort; whatshap for read-based)
3. Run WASP: find_intersecting_snps, re-align, filter_remapped_reads, merge -- mandatory
4. Run GATK ASEReadCounter on WASP-filtered BAM
5. Filter to SNPs with totalCount >= 30 (sufficient power for moderate effects)
6. Within-peak aggregation: pooled binomial test per peak
7. For cohort-scale caQTL: combine MatrixEQTL on peak counts + RASQUAL joint analysis
8. Apply BH-FDR; effect size >= 0.2 (30:70 imbalance) for biologically meaningful
9. Cross-validate against deep-learning-atac variant effect predictions
10. Report cis-regulatory variants with high-confidence triple-evidence rule

## Method Decision Quick Reference

| Setting | Method |
|---------|--------|
| Single individual, has genotypes | WASP + GATK ASEReadCounter |
| Cohort N >= 100 | WASP + GATK + MatrixEQTL on peak counts |
| Cohort N = 20-100 | WASP + RASQUAL (joint power gain) |
| No genotypes | QuASAR |
| Validating GWAS variant | Lookup het samples; aggregate ASE |
| Trios / family designs | Per-trio phasing + per-individual ASE |

## Threshold Quick Reference

| Test | Threshold |
|------|-----------|
| Per-SNP coverage | totalCount >= 30 |
| Per-peak SNP aggregation | snp_count >= 2 |
| Effect size | abs(ref_frac - 0.5) >= 0.2 (30:70 imbalance) |
| Significance | BH-adjusted p < 0.05 |
| Cohort caQTL | p < 1e-5 (genome-wide) |
| Allele frequency filter | 0.1 < AF < 0.9 (avoids rare-allele extremes) |

## Tips

- WASP filtering is MANDATORY. There are no exceptions. Reference-allele mapping bias is universal.
- Phase genotypes BEFORE WASP. Unphased het sites can't assign reads to haplotypes.
- Per-SNP ASE has limited power at standard depth (10-30 reads); aggregate within peaks for biology-level inference.
- RASQUAL gains 1.5-3x power over MatrixEQTL alone at moderate N (20-100); standard for caQTL studies.
- Effect size 0.2 (30:70) is the typical lower bound for biological relevance; 0.4 (10:90) is strong.
- For GWAS variant function, look up het carriers in the cohort and aggregate ASE at the lead variant; meta-analyze across individuals.
- iPSC line genotype validation (single-individual) often uses ASE at known SNPs to confirm parental origin.
- Cross-validate predicted vs observed: chromBPNet variant effects (atac-seq/deep-learning-atac) at same SNPs should correlate with observed allelic imbalance.
- Multi-allelic sites are problematic for GATK ASEReadCounter; pre-filter with bcftools to biallelic only.
- Reference allele bias is ~1-5% at typical aligner stringencies; WASP corrects this systematically. Without WASP, expect 0.51-0.55 reference fraction even at null sites.
- For high-confidence cis-regulatory variant claims, require: WASP-filtered ASE FDR < 0.05 + effect size >= 0.2 + cohort caQTL p < 1e-5 + accessibility peak overlap.

## Related Skills

- atac-seq/atac-peak-calling - Generate peaks for ASE aggregation
- atac-seq/consensus-peakset - Cohort consensus
- atac-seq/differential-accessibility - Cohort-level (cis + trans) DA
- atac-seq/deep-learning-atac - Predicted vs observed variant effects
- atac-seq/enhancer-gene-linking - Map ASE-supported variants to target genes
- variant-calling/vcf-basics - VCF input
- variant-calling/joint-calling - Cohort genotypes
- phasing-imputation/haplotype-phasing - Phase before WASP
- causal-genomics/fine-mapping - caQTL for fine-mapping
- population-genetics/association-testing - GWAS context
