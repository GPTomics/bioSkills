# Allele-Specific Binding - Usage Guide

## Overview

Identify variants with allelic effects on TF or histone modification binding from heterozygous-variant ChIP-seq. Covers WASP (mandatory reference-bias filter), RASQUAL (joint cis-QTL + bias-aware testing), BaalChIP (Bayesian beta-binomial with copy-number-aware overdispersion), and AlleleSeq (personalized diploid genome). Embeds the three universal pitfalls (reference bias, imprinted loci, X-inactivation, copy-number imbalance) and integrates with downstream caQTL / bQTL / fine-mapping.

## Prerequisites

```bash
# WASP (mandatory upstream)
git clone https://github.com/bmvdgeijn/WASP.git
# Requires Python 3, pysam, pytables

# Aligner
conda install -c bioconda bowtie2 samtools bcftools

# Variant calling (for sample-specific hetSNPs)
conda install -c bioconda gatk4
```

```r
BiocManager::install(c('BaalChIP', 'AllelicImbalance'))
```

```bash
# RASQUAL
git clone https://github.com/dg13/rasqual.git
# Build with `make` in repo

# AlleleSeq
# https://github.com/gerstein-lab/AlleleSeq
```

## Quick Start

Tell the agent what to do:
- "Apply WASP reference-bias filter to my ChIP-seq BAM before any ASB analysis (mandatory step)"
- "Run BaalChIP on TNBC cancer samples with copy-number-aware overdispersion"
- "Use RASQUAL for joint cis-QTL + ASB analysis on a population of 50 samples"
- "Build a personalized diploid genome with AlleleSeq for a single-sample ASB analysis"
- "Filter imprinted loci and X-inactivated chrX before reporting ASB"
- "Validate chromBPNet variant effect predictions against ASB measurements"

## Example Prompts

### WASP pipeline
> "Run the WASP mapping pipeline on a bowtie2-aligned ChIP-seq BAM with my sample's heterozygous SNP VCF. Apply the three steps: find intersecting SNPs, re-map swapped reads, filter discordant mappings. Output WASP-filtered BAM (will lose 22-31% of reads)."

### BaalChIP for cancer
> "I have FOXA1 ChIP-seq in HCC1395 (TNBC). Run BaalChIP with copy-number-aware overdispersion using my ASCAT CNV calls. Output ASB sites with posterior > 0.95."

### RASQUAL for population
> "I have 50 individuals' CTCF ChIP-seq with phased genotypes. Run RASQUAL for cis-QTL mapping with built-in bias correction (`phi` parameter). Window 250 kb around each peak."

### AlleleSeq personalized genome
> "Build a personalized diploid genome from my sample's phased VCF. Align ChIP reads to maternal and paternal copies separately, then merge for binomial ASB test."

### Cross-validation with DL
> "I have chromBPNet variant predictions and BaalChIP ASB measurements. Compare: variants with |log2_fc| > 1 in chromBPNet should show ASB; variants with no ASB but predicted strong effect may be in extrapolation regime."

### Filter universal pitfalls
> "Before reporting ASB: (1) verify WASP was applied; (2) filter imprinted loci (H19, IGF2, MEG3, KCNQ1OT1); (3) exclude chrX in female samples; (4) for cancer, exclude CN-altered regions or use BaalChIP with CN file."

## What the Agent Will Do

1. **Verify hetSNP source**: sample-specific genotype VCF (from GATK / DeepVariant) preferred over population SNP panel
2. **Apply WASP filter**: mandatory; produces unbiased BAM with 22-31% read loss
3. **Filter known biology:**
   - Imprinted loci (H19, IGF2, etc.) from imprinted-genes BED
   - chrX in female samples (X-inactivation)
   - Copy-number-altered regions (for cancer; or use BaalChIP)
4. **Choose ASB method:**
   - Single sample, normal: BaalChIP or AlleleSeq
   - Cancer sample with CNV: BaalChIP with CN file
   - Population (>20 samples): RASQUAL for joint cis-QTL
   - Need maximum bias correction: AlleleSeq with phased genome
5. **Run test** with appropriate parameters
6. **Output**: per-hetSNP ASB calls with allelic ratio, posterior, FDR
7. **Cross-validate** with deep-learning predictions (chromBPNet, EnFormer) where applicable
8. **Document**: WASP version + filter loss, ASB method + parameters, imprinted/X-filter applied, CN integration

## Tips

- **WASP is mandatory.** Skipping it produces reference-allele-skewed false positives.
- **Use sample-specific hetSNP VCF, not population panel.** Sample-specific variants may not be in 1KG.
- **WASP drops 22-31% of reads.** Plan sequencing depth accordingly; ASB power requires roughly 2× the depth needed for peak calling.
- **For cancer, BaalChIP with CN file is the only rigorous option.** CN imbalance changes effective allele dose; methods that ignore it produce false ASB calls.
- **Imprinted loci and X-inactivation are the most common artifact sources.** Always filter.
- **AlleleSeq is conceptually cleanest but expensive.** Personalized diploid genome construction + alignment doubles compute cost.
- **RASQUAL `phi` parameter** is the per-feature bias estimate; can be cross-checked against WASP filter loss.
- **Validate DL variant predictions against measured ASB.** chromBPNet predicts in counterfactual; ASB measures in actual sample. Concordance increases confidence.
- **ASB ≠ caQTL.** ASB is single-sample at hetSNPs; caQTL is population-level association of any variant with chromatin accessibility / binding.

## Troubleshooting

### Many ASB calls still REF-biased

WASP not applied or sample-specific hetSNPs missing from WASP SNP file. Build WASP table from sample's own genotype VCF.

### Most ASB calls at imprinted loci

Imprinted loci not filtered. Apply imprinted-genes BED before reporting.

### Most ASB at chrX (female sample)

X-inactivation. Filter chrX or use X-inactivation-aware methods.

### Bimodal allelic ratios in cancer

Copy-number imbalance. Use BaalChIP with CN file OR exclude CN-altered regions.

### BaalChIP "no hetSNPs in peaks"

hetSNP / peak chromosome naming mismatch. Standardize chr prefixes.

### RASQUAL convergence failure

1. Sparse hetSNPs per peak → require more samples or combine peaks
2. Strong CN imbalance → use BaalChIP
3. Imputation quality poor → require r² > 0.8

### AlleleSeq "diploid genome too large"

Many large SVs in genome. Use small-variant-only VCF; exclude SV-rich regions.

### Low statistical power

ASB needs ≥20 reads per allele to detect a 2:1 ratio at p<0.05. Increase sequencing depth (2× normal ChIP for ASB power) or combine replicates.

### Disagreement with chromBPNet predictions

1. ASB sample has low coverage at variant → boost depth
2. chromBPNet model may be in extrapolation regime → check ensemble agreement
3. Variant in context the model didn't see → caution

## Related Skills

- chip-seq/peak-calling - Peak calling upstream
- chip-seq/chipseq-qc - QC before ASB
- chip-seq/chip-deep-learning - Cross-validate ASB with DL variant predictions
- chip-seq/peak-annotation - Annotate ASB variants to genes / cCREs
- atac-seq/allele-specific-accessibility - Parallel ATAC ASB workflow
- causal-genomics/fine-mapping - ASB as fine-mapping orthogonal evidence
- variant-calling/variant-annotation - Annotate hetSNPs
- phasing-imputation/haplotype-phasing - Required for AlleleSeq
