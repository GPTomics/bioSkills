# Enhancer-Gene Linking - Usage Guide

## Overview

Predict which gene each distal enhancer regulates using ABC, ENCODE-rE2G, HiChIP, or Cicero. Choose the right method by available data: ATAC + H3K27ac + Hi-C/Micro-C enables ABC and ENCODE-rE2G; ATAC alone falls back to Cicero. Validate predictions against CRISPRi-FlowFISH gold-standard. Build cell-type-specific regulatory maps for fine-mapping causal SNPs to target genes or for therapeutic target discovery.

## Prerequisites

```bash
# ABC pipeline
git clone https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
pip install -r ABC-Enhancer-Gene-Prediction/requirements.txt

# ENCODE-rE2G
git clone https://github.com/EngreitzLab/ENCODE_rE2G
pip install snakemake

# HiChIP loop calling
conda install -c bioconda fithichip hicpro samtools bedtools
```

```r
# Cicero (atac-seq/co-accessibility skill)
BiocManager::install(c('cicero', 'monocle3', 'GenomicInteractions'))
```

Inputs:
- ATAC-seq bigWig (RPGC-normalized)
- H3K27ac ChIP-seq bigWig (matched to ATAC; RPM-IP normalized)
- Hi-C/Micro-C contact matrix (cooler or hic format) at 5-10 kb resolution
- Gene annotation BED (RefSeq protein-coding TSSs)
- ENCODE blacklist v2

## Quick Start

Tell your AI agent what you want to do:
- "Run ABC on K562 ATAC + H3K27ac + Micro-C; threshold at ABC.Score >= 0.02"
- "Run ENCODE-rE2G with the K562 pre-trained model; report scores >= 0.5"
- "I have ATAC but no Hi-C. Use ABC's average HiC fallback OR fall back to Cicero co-accessibility"
- "Cross-validate ABC predictions against published Fulco 2019 CRISPRi-FlowFISH ground truth"
- "Combine ABC + ENCODE-rE2G + HiChIP -- intersection is high-confidence"
- "Score the variant effect at fine-mapped GWAS SNPs at predicted enhancers (link to deep-learning-atac)"

## Example Prompts

### Standard ABC
> "Run ABC on K562 with cell-type-matched Micro-C at 5 kb resolution. Pre-filter promoter regions from candidate enhancers. Threshold ABC.Score >= 0.02 for high-confidence predictions."

### ENCODE-rE2G
> "Use ENCODE-rE2G with the GM12878-trained model on my B-cell ATAC. Report predictions with re2g_score >= 0.5 and effective HiC; cross-check against ABC."

### No Hi-C Fallback
> "I only have ATAC + H3K27ac, no Hi-C. Run ABC with the average HiC fallback (5 ENCODE cell types pooled) and acknowledge degraded performance in methods. Or fall back to Cicero co-accessibility."

### Method Combination
> "Run both ABC and ENCODE-rE2G; intersect with HiChIP H3K27ac loops at FDR < 0.05. Report the triple-overlapping enhancer-gene pairs as high-confidence and the union as exploratory."

### CRISPRi Validation
> "Compare my ABC predictions in K562 against the Fulco 2019 CRISPRi-FlowFISH catalog. Report sensitivity and specificity at threshold ABC.Score >= 0.02."

### Variant-to-Gene
> "I have a GWAS lead SNP in a non-coding region. Identify the candidate enhancer (ATAC + H3K27ac peak overlap), then use ABC to find the most likely target gene."

## What the Agent Will Do

1. Verify available inputs (ATAC, H3K27ac, Hi-C, gene annotations)
2. Choose method by data availability (decision tree above)
3. For ABC: define candidate enhancers (non-promoter ATAC peaks), compute Activity (ATAC * H3K27ac), Contact (Hi-C), apply ABC formula
4. For ENCODE-rE2G: choose pre-trained model by cell-type proximity; run Snakemake pipeline
5. For HiChIP: call loops at FDR < 0.05 + contact count >= 5; intersect with ABC
6. Cross-validate against CRISPRi-FlowFISH ground truth where available
7. Threshold at recommended cutoffs (ABC >= 0.02; rE2G >= 0.5)
8. Report intersection of methods as high-confidence; union as exploratory
9. Document cell-type proxies if exact match not available

## Decision Quick Reference

| Data | Method |
|------|--------|
| ATAC + H3K27ac + Hi-C/Micro-C (matched) | ABC OR ENCODE-rE2G |
| ATAC + H3K27ac, no Hi-C | ABC with avg HiC fallback OR ENCODE-rE2G no-hic model |
| ATAC only | Cicero (atac-seq/co-accessibility) |
| ATAC + H3K27ac HiChIP | FitHiChIP + ABC, intersect |
| ENCODE cell type | Pre-computed ENCODE-rE2G predictions (download) |
| Need experimental validation | CRISPRi-FlowFISH design |

## Threshold Quick Reference

| Method | Standard | Stringent |
|--------|---------|-----------|
| ABC.Score | >= 0.02 | >= 0.04 |
| ENCODE-rE2G | >= 0.5 | >= 0.8 |
| FitHiChIP loops | FDR < 0.05 + count >= 5 | FDR < 0.01 + count >= 10 |
| Cicero coaccess | > 0.25 | > 0.5 |

## Tips

- ABC requires matched Hi-C from the same cell type; cross-celltype Hi-C performs degraded but is documented as an acceptable fallback (Fulco 2019).
- ENCODE-rE2G's logistic regression weights are cell-type-specific; using a wrong cell type's weights produces calibration errors. Match by tissue similarity.
- ABC's Activity term combines ATAC and H3K27ac as the geometric mean of normalized signals; both must be on the same scale (use RPGC for ATAC, RPM-IP for H3K27ac).
- Cicero co-accessibility is NOT enhancer-gene linking. Cicero predicts co-varying peak pairs; ABC predicts target genes. Concordance is ~30-50%.
- HiChIP H3K27ac directly measures enhancer-promoter loops; orthogonal to accessibility-based methods.
- For high-confidence reporting, agree across two methods. Single-method calls are exploratory.
- CRISPRi-FlowFISH ground truth catalogs (Fulco 2019, Gasperini 2019, Schraivogel 2020) cover ~10 cell types with ~thousands of enhancer-gene pairs; benchmark predictions against these.
- For GWAS variant interpretation, combine with deep-learning-atac (chromBPNet) variant effect predictions; SNP at predicted enhancer + ABC target gene + chromBPNet effect = strong functional hypothesis.
- ENCODE 4 published pre-computed ENCODE-rE2G predictions for many cell types -- check availability before retraining.
- EpiMap and GeneHancer are useful baselines but cell-type-agnostic; do not use as primary calls for cell-type-specific biology.

## Related Skills

- atac-seq/co-accessibility - Cicero (ATAC-only fallback)
- atac-seq/atac-peak-calling - Enhancer candidate generation
- atac-seq/consensus-peakset - Fixed-width enhancer regions
- atac-seq/deep-learning-atac - chromBPNet variant effect prediction
- atac-seq/single-cell-atac - Per-cell-type scATAC inputs
- hi-c-analysis/loop-calling - Hi-C / Micro-C contact prediction
- hi-c-analysis/contact-pairs - Hi-C / Micro-C input
- chip-seq/peak-calling - H3K27ac peaks
- gene-regulatory-networks/scenic-regulons - TF -> target downstream
