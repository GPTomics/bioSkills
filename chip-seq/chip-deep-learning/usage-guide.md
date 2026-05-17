# Deep Learning for ChIP-seq - Usage Guide

## Overview

Apply base-resolution deep learning models to ChIP-seq, ChIP-nexus, and CUT&RUN data for variant-effect prediction, motif syntax discovery, and sequence-only signal prediction. Covers BPNet (Avsec 2021), chromBPNet (Pampari Avsec 2024/2025), EnFormer (Avsec 2021), DeepSEA, and the JASPAR 2026 Deep Learning collection (1259 precomputed BPNet ChIP models). Embeds in silico mutagenesis (counterfactual variant prediction), TF-MoDISco motif discovery from attribution scores, and model-ensemble uncertainty estimation.

## Prerequisites

```bash
# chromBPNet (modern; handles ChIP and ATAC)
pip install chrombpnet
# or conda for full env: conda install -c bioconda chrombpnet

# BPNet (canonical for ChIP-nexus)
pip install bpnet

# EnFormer
pip install enformer-pytorch

# TF-MoDISco
pip install tfmodisco-lite

# DeepLIFT / SHAP for attribution
pip install shap
```

```bash
# CUDA-capable GPU recommended (16+ GB VRAM for training)
# CPU inference possible for precomputed models
```

## Quick Start

Tell the agent what to do:
- "Train chromBPNet on my CTCF ChIP-seq peaks for variant effect prediction"
- "Predict the effect of a SNP on TF binding using in silico mutagenesis with chromBPNet"
- "Run TF-MoDISco on attribution scores from a trained BPNet model to discover motif syntax"
- "Use a precomputed JASPAR 2026 BPNet model to score a list of GWAS variants for TF binding effects"
- "Apply EnFormer for variant effect prediction with 200 kb receptive field (captures distal regulatory effects)"
- "Compare chromBPNet and EnFormer variant predictions for a list of fine-mapped SNPs"

## Example Prompts

### Variant effect prediction
> "For a list of fine-mapped SNPs from GWAS, predict the effect on TF binding using chromBPNet. Report |log2_fc| > 1 as strong-effect, 0.3-1 as moderate, <0.3 as weak. Concordance with EnFormer increases confidence."

### Train chromBPNet from scratch
> "Train chromBPNet on FOXA1 ChIP-seq in HepG2 cells. First train a bias model on non-peak regions, then train the main model. Use ENCODE-matched non-peak background."

### TF-MoDISco motif discovery
> "Compute DeepLIFT attribution scores from a trained BPNet model, then run TF-MoDISco to discover motif syntax including cooperative TF interactions. Compare to PWM motifs from MEME-ChIP."

### Use JASPAR precomputed
> "Download the BPNet model for FOXA1 from JASPAR 2026 Deep Learning collection. Use it to predict variant effects on FOXA1 binding for my list of SNPs without training."

### Cross-model validation
> "Predict the same variant effect with chromBPNet and EnFormer. If they agree (both |log2_fc| > 1 in same direction), high-confidence functional variant. If they disagree, low-confidence."

### Long-range prediction
> "Use EnFormer's 200 kb receptive field to predict whether a distal enhancer SNP affects target gene expression. Native EnFormer outputs include CAGE-seq signal at TSS."

### Ensemble uncertainty
> "Train 5 chromBPNet replicate models with different seeds. Report ensemble mean + std for each variant. High std indicates model extrapolation; low std with consistent prediction is high-confidence."

## What the Agent Will Do

1. **Choose model** based on goal:
   - Single TF, ChIP-nexus data: BPNet
   - TF or histone mark, regular ChIP, with bias factorization: chromBPNet
   - Long-range variant effects: EnFormer
   - Pre-trained models without training cost: JASPAR 2026 deep-learning collection
2. **For training:**
   - Verify ≥5000 high-confidence peaks per TF
   - Train bias model on non-peak background
   - Train main model with cross-validation folds (chromBPNet uses 5 folds by default)
   - Monitor training: Spearman correlation between predicted and observed profiles should be ≥0.6
3. **For variant prediction:**
   - Encode ref and alt sequence (2114 bp window for chromBPNet)
   - Predict counts and profile for each
   - Compute log2_fc = log2(alt_counts / ref_counts)
   - Apply ensemble for uncertainty
4. **For motif discovery:**
   - Compute attribution scores (DeepLIFT or SHAP)
   - Run TF-MoDISco with appropriate parameters
   - TOMTOM the discovered motifs against JASPAR / HOCOMOCO for known matches
   - Identify motifs not in known DBs as candidate novel syntax
5. **Output:**
   - Variant-effect TSV with predicted log2_fc per variant
   - Attribution score plots (sequence logo of importance)
   - TF-MoDISco motif logos
6. **Document**: model version, training data, hyperparameters, ensemble size, validation metrics

## Tips

- **|log2_fc| > 1 is the strong-effect threshold per chromBPNet 2024 paper.** Use as initial filter.
- **Concordance between chromBPNet and EnFormer increases confidence.** Independent architectures with different receptive fields agreeing is strong evidence.
- **Ensemble of 5-10 models for uncertainty.** Single-model predictions can be unstable; ensemble disagreement flags low-confidence variants.
- **Bias model is critical for chromBPNet.** Must be trained on matched assay (ATAC bias for ATAC; ChIP bias for ChIP). Wrong bias = noisy predictions.
- **TF-MoDISco background matters.** Random genomic sequences contain other TF sites; shuffled-input or IgG-control sequences are better.
- **Training cost: 1-3 GPU days for chromBPNet, longer for EnFormer.** Use precomputed JASPAR models when possible.
- **JASPAR 2026 has 1259 BPNet ChIP models.** Lowest-effort variant prediction; covers 240 TFs from ENCODE.
- **EnFormer is tissue-aggregated.** Cell-line-specific predictions require fine-tuning OR chromBPNet trained on cell-line data.
- **Variant must be in ACGT only.** Models cannot handle ambiguous nucleotides; filter N positions.
- **Need ≥5000 peaks for stable BPNet training.** Fewer peaks → fall back to PWM-based motif analysis.

## Troubleshooting

### Training accuracy low (Spearman < 0.5)

1. Insufficient peaks → need ≥5000 high-confidence
2. Wrong bias model → train on matched assay
3. Library quality poor → check chipseq-qc upstream
4. Hyperparameters → check chromBPNet recommended config

### Variant predictions noisy / unstable

1. Single model; train ensemble of 5-10 with different seeds
2. Variant outside training distribution; extrapolation; experimental validation needed
3. Variant in low-quality region (low mappability, repeat); predictions unreliable

### TF-MoDISco produces empty motif list

1. FDR too strict; lower `target_seqlet_fdr` to 0.10
2. Attribution scores too noisy; train model longer; use ensemble attribution
3. Window size mismatch; adjust `sliding_window_size`

### GPU out of memory

1. Reduce batch_size
2. Use gradient checkpointing
3. For EnFormer, use the 5x smaller architecture variant

### `import chrombpnet` fails

TensorFlow version mismatch. Use the conda env from chrombpnet repo with pinned `tf==2.13`.

### EnFormer predictions for short variant don't capture biology

Variant >100 kb from any annotated regulatory element; effect may be too distal or not captured by training data. Combine with chromBPNet for local effects.

### Chromosome naming mismatch with model

JASPAR / EnFormer models expect specific genome assembly and chromosome conventions. Verify hg38 chromosome naming (chr1 vs 1).

## Related Skills

- chip-seq/peak-calling - Source peaks for training
- chip-seq/motif-analysis - PWM-based motif analysis (complementary)
- chip-seq/allele-specific-binding - Validate DL predictions against ASB
- atac-seq/deep-learning-atac - ATAC-specific chromBPNet workflow
- atac-seq/footprinting - TOBIAS footprints as comparison
- causal-genomics/fine-mapping - Variant prioritization with functional scoring
- machine-learning/biomarker-discovery - DL variant scores as biomarker features
- machine-learning/model-validation - Cross-validation and ensemble methods
