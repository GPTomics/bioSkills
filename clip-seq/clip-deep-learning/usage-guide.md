# CLIP Deep Learning - Usage Guide

## Overview

Predict RBP binding from RNA sequence with deep neural networks. Modern models: RBPNet (2024, sequence-to-CL distribution at single-nt), RNAProt (RNN classifier; AUC 87-89%), GraphProt2 (GCN with structure), DeepCLIP, DeepRiPe. Use for variant-effect prediction, in silico binding-site discovery, and systematic RBP comparison. Training requires GPU; chromosome-split prevents data leakage; background must match foreground transcript region. Pretrained models exist for ~150 ENCODE RBPs.

## Prerequisites

```bash
conda install -c pytorch pytorch torchvision
pip install rbpnet rnaprot graphprot2 deepclip biopython
```

## Quick Start

Tell your AI agent:
- "Predict per-base RBP binding from sequence with RBPNet"
- "Train custom RNAProt on my eCLIP peaks with GC-matched 3' UTR background"
- "Variant effect at heterozygous SNP using RBPNet log2 FC"
- "Apply pretrained DeepRiPe to a list of GWAS variants"
- "GraphProt2 with structure for RBPs that prefer specific structural contexts"
- "Why is my model AUC 0.95 but predictions fail on novel sequences? Data leakage"
- "Saliency map to recover motif - use TF-MoDISco for clean output"

## Example Prompts

### Variant-Effect Prediction

> "RBPNet for TARDBP on heterozygous SNP rs12345 - log2 FC > 1 = strong effect"

> "Genome-wide variant scoring for GWAS variants in eCLIP regions"

### Custom Training

> "RNAProt training on my new CLIP data - 50 epochs, 1:1 balanced sampling, chromosome-split"

> "Use train chr1-20, test chr21-22 to prevent gene-neighbor leakage"

### Production Inference

> "DeepCLIP for fast inference on 1M sequences"

> "RBPNet per-base score for a list of 3' UTR sequences"

### Interpretation

> "Saliency on trained DeepRiPe - extract motif with TF-MoDISco"

> "GraphProt2 with RNAfold structure ensemble for structure-dependent RBP"

### Transfer Learning

> "RNA-FM foundation model fine-tuned on my eCLIP - emerging approach"

> "Pretrained RBP not available - try transfer from closest paralog"

### Diagnostics

> "Variant effect different across windows - use model native 256 nt"

> "GraphProt2 underperforms - structure ensemble not provided"

> "Pretrained for related RBP works - paralog transfer valid"

## What the Agent Will Do

1. Choose model: RBPNet (per-base), RNAProt (binary), GraphProt2 (structure), DeepCLIP (fast), DeepRiPe (multi-modal)
2. Verify pretrained model available; if not, train custom on user's CLIP peaks
3. For custom training: prepare GC-matched background; chromosome-split train/test
4. Apply to query sequences: per-base for RBPNet; per-sequence for binary models
5. For variant effect: predict reference + alt; compute log2 FC summed over 50 nt window
6. Interpret with saliency + TF-MoDISco for motif extraction
7. Cross-validate with eCLIP overlap, RBNS Kd, related RBP predictions
8. Flag pitfalls: data leakage, GPU memory, structure input, class imbalance

## Tips

- **RBPNet is the modern per-base model.** Single-nt resolution.
- **Chromosome-split prevents leakage.** Random shuffle leaks gene-neighbor context.
- **Background must match foreground region.** 3' UTR peaks - 3' UTR background.
- **Balanced sampling matters.** RNAProt does this automatically.
- **GPU required.** CPU is 100x slower.
- **Pretrained models are ~150 ENCODE RBPs.** Not all proteins.
- **Variant-effect window 256 nt for RBPNet.** Native model size.
- **GraphProt2 needs RNAfold pre-compute.** Not automatic.
- **Vanilla saliency is noisy.** Use TF-MoDISco or DeepLIFT.
- **Variant log2 FC > 1.** Strong effect threshold.
- **RNA foundation models emerging.** RNA-FM, RNAErnie for transfer learning.

## Related Skills

- clip-seq/clip-motif-analysis - Classical motif alternative
- clip-seq/crosslink-site-detection - CL sites for sequence-to-signal
- clip-seq/clip-peak-calling - Peak BEDs for training
- causal-genomics/mendelian-randomization - Variant-effect feeds MR
- causal-genomics/fine-mapping - Variant prioritization
- machine-learning/model-validation - Train/test methodology
- machine-learning/prediction-explanation - Saliency methods
- machine-learning/biomarker-discovery - ML framework
