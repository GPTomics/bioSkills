# Deep Learning for ATAC-seq - Usage Guide

## Overview

Sequence-based neural networks for ATAC-seq: chromBPNet (bias-corrected per-base accessibility), BPNet (foundational profile model), scBasset (per-cell sequence model), Enformer (long-context Transformer), and Borzoi (multi-tissue chromatin + RNA). Covers in silico variant effect prediction at GWAS / rare-variant SNPs, motif discovery via DeepLIFT + TF-MoDISco, deep bias correction beyond TOBIAS k-mer model, and per-cell-type accessibility prediction for unobserved states.

## Prerequisites

```bash
# Core deep learning frameworks
pip install tensorflow==2.13 torch torchvision torchaudio

# Tools
pip install chrombpnet bpnet-lite tangermeme modisco-lite captum kipoi
```

GPU with >= 16 GB VRAM recommended for training (A100/V100); CPU acceptable for inference with pre-trained models. Training data: deduplicated BAM, peaks/non-peaks BED, reference FASTA, chrom sizes.

## Quick Start

Tell your AI agent what you want to do:
- "Score 100 GWAS SNPs for chromatin effect using a pre-trained chromBPNet K562 model"
- "Run chromBPNet pipeline (bias + accessibility) on a new ATAC dataset"
- "Compute DeepLIFT contributions and run TF-MoDISco-lite to discover de novo motifs"
- "Use chromBPNet bias-corrected bigWig as input to TOBIAS instead of ATACorrect"
- "Run Enformer pre-trained on a 196 kb window to predict variant effects in distal regulation"
- "Train scBasset on my scATAC for per-cell accessibility prediction"

## Example Prompts

### Variant Effect Prediction
> "I have 50 GWAS lead SNPs in autoimmune loci. Use chromBPNet pre-trained models for the closest cell types (CD4 T, B cells), score ref vs alt with tangermeme marginal prediction, and report log2FC magnitudes. SNPs with abs(log2FC) > 1 are strong-effect candidates."

### chromBPNet Training
> "Train a per-cell-type chromBPNet on my 80M-read primary T cell ATAC. Generate splits, train bias model from non-peaks, train accessibility model with bias correction, then use the corrected bigWig for downstream footprinting."

### TF-MoDISco Motif Discovery
> "Compute DeepLIFT contributions for all peaks in my chromBPNet model output, run TF-MoDISco-lite with sliding_window_size=20 and target_seqlet_fdr=0.05, then check motif clusters against JASPAR for known TFs."

### chromBPNet Bias Correction Replacement
> "Use chromBPNet's bias-corrected bigWig output as input to `TOBIAS ScoreBigwig` (formerly `FootprintScores`), skipping `TOBIAS ATACorrect`. Compare aggregate CTCF footprint between chromBPNet-corrected and TOBIAS-corrected to verify both produce clean V-shapes."

### Enformer Variant Effect
> "For a SNP in a non-coding distal element, use Enformer pre-trained model: predict for ref and alt 196 kb windows centered on the SNP, report differential predictions across HepG2/K562/GM12878 tracks for cell-type-specific effects."

### scBasset for scATAC
> "Train scBasset on my 50k-cell scATAC-seq dataset; use the per-cell projection layer to identify cells with anomalous accessibility patterns and compare cluster discrimination to chromVAR."

## What the Agent Will Do

1. Determine if classical pipelines suffice, or whether DL is genuinely needed (variant effect, cross-cell-type, fine-grained motif)
2. Check for pre-trained model availability at encodeproject.org/atac-seq before training new
3. Verify input data: dedup, peaks, chrom sizes, reference FASTA
4. Prepare splits (train / val / test chromosomes)
5. Train bias model first (from naked-DNA control if available; otherwise non-peak background)
6. Train accessibility model with bias correction
7. For variant effect: tangermeme marginal prediction at SNPs of interest
8. For motif discovery: DeepLIFT contributions -> TF-MoDISco-lite clustering
9. For sc state interpolation: scBasset training on metacells
10. Document model versions and pre-trained model sources for reproducibility

## Tool Decision Quick Reference

| Goal | Tool |
|------|------|
| Score GWAS SNPs in known cell types | Pre-trained chromBPNet + tangermeme |
| Score GWAS SNPs in novel cell types | Train chromBPNet first |
| Long-range distal regulation variants | Enformer or Borzoi (196 kb context) |
| Motif discovery from new ATAC data | chromBPNet + DeepLIFT + TF-MoDISco-lite |
| Bias correction for footprinting | chromBPNet bias-corrected bigWig as TOBIAS input |
| Per-cell accessibility prediction | scBasset |
| Cross-cell-type accessibility | Enformer pre-trained |

## Compute Quick Reference

| Task | Hardware | Wall time |
|------|---------|-----------|
| chromBPNet training (1 cell type) | A100, 80 GB | ~24 h |
| chromBPNet inference 1M variants | A100 | ~4 h |
| Enformer pre-trained inference | V100+ | ~30 min / 100k variants |
| Borzoi training | A100, 250 GB | ~7 days |
| scBasset training (10k cells) | V100, 32 GB | ~12 h |
| TF-MoDISco on 1M peaks | CPU 32 cores | ~6 h |

## Tips

- For most ATAC analysis, classical pipelines (MACS3 + DiffBind + TOBIAS) remain primary. Deep learning is justified when variant interpretation, cross-cell-type prediction, or fine motif discovery are the goal.
- Pre-trained models exist for ENCODE cell types (K562, GM12878, HepG2, HUVEC, IMR90 ATAC). Check the chromBPNet model zoo before training.
- Bias model and accessibility model are trained sequentially; the bias must be trained first because the accessibility model uses bias-corrected residuals.
- chromBPNet's bias model is mostly cell-type-invariant. A K562 bias model used on T cells works with degraded performance; better than nothing.
- DeepLIFT vs Integrated Gradients give different attribution patterns. Use DeepLIFT (chromBPNet default) unless documented reasons exist for IG.
- Variant effect prediction window matters: short windows (~1-2 kb) for proximal regulation; Enformer/Borzoi (196 kb) for distal.
- TF-MoDISco-lite (its interface is TF-MoDISco v2) is the maintained version; the original kundajelab/tfmodisco (v0.x) is for replication only.
- scBasset's per-cell projection is unstable below 100 cells per cluster; aggregate cells before training.
- Marginal effects (ref vs alt at the SNP only) are different from in silico mutagenesis (saturation across all positions). Choose by question.
- Borzoi is newer than Enformer with better chromatin + RNA joint prediction; use for variant effects on RNA via chromatin linkage.
- For high-confidence variant prediction, agree across two approaches (chromBPNet + Enformer/Borzoi). Single-tool calls are exploratory.

## Related Skills

- atac-seq/atac-peak-calling - Classical input pipeline
- atac-seq/footprinting - Use chromBPNet bias correction as TOBIAS alternative
- atac-seq/motif-deviation - chromVAR vs scBasset for sc TF activity
- atac-seq/single-cell-atac - scBasset integration
- atac-seq/enhancer-gene-linking - Variant effect feeds enhancer scoring
- atac-seq/allele-specific-accessibility - DL predictions vs observed allelic imbalance
- causal-genomics/fine-mapping - Downstream use of variant effect scores
- machine-learning/biomarker-discovery - General ML patterns
- gene-regulatory-networks/scenic-regulons - Motif discovery + TF networks
