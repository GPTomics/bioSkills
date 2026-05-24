---
name: bio-chipseq-chip-deep-learning
description: Trains and applies base-resolution deep learning models on ChIP-seq / ChIP-nexus / CUT&RUN data. Uses BPNet (Avsec 2021 Nat Genet 53:354; soft motif syntax from ChIP-nexus), chromBPNet (Pampari A et al 2025 Nat Genet; bias-factorized base-resolution profiles), EnFormer (Avsec 2021 Nat Methods 18:1196; 196 kb input, ~100 kb effective receptive field), DeepSEA (Zhou 2015; multi-task CNN), and JASPAR 2026 deep-learning collection (1259 BPNet ChIP models). Performs in silico mutagenesis for variant-effect prediction, DeepLIFT/Grad attribution, and TF-MoDISco motif discovery from attribution scores. Use when predicting variant effects on TF binding, discovering soft motif syntax / cooperativity, integrating ChIP-seq with sequence-only predictions, or applying precomputed JASPAR Deep Learning models to new variants.
tool_type: python
primary_tool: chrombpnet
---

## Version Compatibility

Reference examples tested with: chrombpnet 0.1.7+, BPNet 0.0.23+, TF-MoDISco-lite 2.0+, EnFormer (Avsec lab Colab + DeepMind release), tensorflow 2.13+, pytorch 2.0+, JASPAR 2026 deep-learning collection (released 2025).

# Deep Learning for ChIP-seq

**"Predict TF binding from sequence and quantify variant effects on binding"** -> Train base-resolution convolutional / transformer models on ChIP-seq / ChIP-nexus / CUT&RUN profiles; predict reference and alternate-allele binding profiles for variants; extract motif syntax via TF-MoDISco from sequence-attribution scores.

- Python (modern): chrombpnet (bias-factorized; ATAC/DNase/ChIP)
- Python (canonical TF ChIP): BPNet (originally for ChIP-nexus; soft motif syntax)
- Python (long-range): EnFormer (Avsec 2021 Nat Methods 18:1196; 196 kb input window, ~100 kb effective receptive field; tissue-aggregated training)
- Python (multi-task): DeepSEA (Zhou 2015; older but still used)
- Precomputed: JASPAR 2026 Deep Learning collection (1259 BPNet ChIP models from ENCODE; 240 TFs)

Deep-learning ChIP-seq models predict signal from sequence; their power is in counterfactual variant prediction (effect on binding from a SNP) and discovery of soft motif syntax that PWMs miss (cooperativity, spacing).

## Model Taxonomy

| Model | Year | Architecture | Receptive field | Best for |
|-------|------|--------------|------------------|----------|
| **BPNet** (Avsec 2021 Nat Genet 53:354) | 2021 | CNN with dilated convolutions | ~1 kb | TF ChIP-nexus / ChIP-exo; base-resolution profile prediction; soft motif syntax |
| **chromBPNet** (Pampari A et al 2025 Nat Genet) | 2025 | Bias-factorized CNN | ~1-2 kb | ATAC/DNase + ChIP base-resolution; bias-corrected variant effects |
| **EnFormer** (Avsec 2021 Nat Methods 18:1196) | 2021 | Transformer | ~100 kb effective receptive field (input window 196 kb) | Long-range regulatory predictions; cross-tissue; variant effects spanning enhancer-gene |
| **DeepSEA** (Zhou 2015) | 2015 | CNN multi-task | 1 kb | Predicts presence/absence across many chromatin features simultaneously |
| **DeepBind** (Alipanahi 2015) | 2015 | CNN binary classifier | ~50-200 bp | TF binding presence (older, less precise than BPNet) |
| **Basset** (Kelley 2016) | 2016 | CNN | ~600 bp | DNase / ATAC accessibility prediction |
| **JASPAR 2026 Deep Learning collection** | 2025 | Precomputed BPNet | ~1 kb | 1259 ENCODE TF ChIP-seq models; 240 TFs; ready-to-use |

## Decision Tree: Which Model

| Goal | Model | Why |
|------|-------|-----|
| Predict variant effect on TF binding (cis-pQTL fine-mapping) | chromBPNet or EnFormer | Both predict ref/alt counterfactuals; chromBPNet base-resolution, EnFormer long-range |
| Discover motif syntax / TF cooperativity from existing ChIP | BPNet (ChIP-nexus data) or chromBPNet (regular ChIP) + TF-MoDISco | Attribution-based motif discovery captures soft syntax PWMs miss |
| Use precomputed model on new variant | JASPAR 2026 deep-learning collection | 1259 BPNet ChIP models ready; no training needed |
| Predict ChIP signal from sequence in a new cell type | EnFormer (cross-tissue training) | Long-range receptive field; multi-tissue training |
| Integrate ATAC + ChIP into single model | chromBPNet | Bias-factorized handles both assays |
| TF binding presence/absence multi-task | DeepSEA | Older but simple multi-output |

## In Silico Mutagenesis Workflow

**Goal:** Predict whether a single-nucleotide variant changes transcription factor or histone modification binding.

**Approach:** Encode reference and alternate-allele sequences in the model's expected window (2114 bp for chromBPNet, centered on variant), predict per-base profile + total counts for each, compute log2 fold change in counts as the variant-effect score. Apply ensemble of 5-10 models for uncertainty.

The most clinically/translationally useful application: predict whether a variant changes TF binding.

```python
import chrombpnet
import numpy as np
import tensorflow as tf

# Load trained chromBPNet model
model = tf.keras.models.load_model('chrombpnet_model.h5', compile=False)

# Reference and alternate-allele sequence around variant (2114 bp window typical)
ref_seq = encode_dna_one_hot('NNN...CCATGNNN...')   # 2114 bp; variant position central
alt_seq = encode_dna_one_hot('NNN...CCAAGNNN...')   # G->A at central position

# Predict base-resolution profiles
ref_profile, ref_counts = model.predict(ref_seq[None, ...])
alt_profile, alt_counts = model.predict(alt_seq[None, ...])

# Variant effect: log2 fold change in predicted total counts
log2_fc = np.log2(alt_counts / ref_counts)
print(f'Variant effect: log2_fc = {log2_fc}')

# |log2_fc| > 1 indicates strong-effect SNP per chromBPNet 2024 paper
# Concordance with EnFormer increases confidence for clinical interpretation
```

**Variant effect interpretation:**
- |log2_fc| > 1: strong effect; binding likely affected
- 0.3 < |log2_fc| < 1: moderate effect; investigate further
- |log2_fc| < 0.3: weak / no effect predicted
- Concordance between chromBPNet and EnFormer increases confidence

## TF-MoDISco for Soft Motif Syntax

Standard PWM-based motif discovery misses:
- Cooperative motif interactions (TF dimers, ETS-RUNX, GATA-TAL)
- Soft motif syntax (variable spacing, weak co-binding)
- Long-range dependencies

TF-MoDISco extracts motifs from deep-learning attribution scores:

```python
import tfmodisco
import shap

# Compute DeepLIFT / DeepSHAP attribution scores
explainer = shap.DeepExplainer(model, background_seqs)
attribution_scores = explainer.shap_values(test_seqs)

# Run TF-MoDISco
modisco_results = tfmodisco.workflow.TfModiscoWorkflow(
    sliding_window_size=21,
    flank_size=10,
    target_seqlet_fdr=0.05,
    seqlets_to_patterns_factory=...
)(task_names=['task1'], contrib_scores={'task1': attribution_scores}, ...)

# Output: motif patterns discovered from attribution (not from PWM matching)
# Often more interpretable than PWM motifs for cooperative TF binding
```

## Training chromBPNet from Scratch

```bash
# Install
pip install chrombpnet

# Train bias model (control regions without TF binding)
chrombpnet bias pipeline \
    -ibam input.bam -d ATAC \
    -g hg38.fa -c chrom.sizes -p peaks.bed \
    -n nonpeaks.bed -fl fold_0.json \
    -b bias_model_h5

# Train main chromBPNet model
chrombpnet pipeline \
    -ibam chip.bam -d ChIP \
    -g hg38.fa -c chrom.sizes -p peaks.bed -n nonpeaks.bed \
    -fl fold_0.json \
    -b bias_model_h5 \
    -o output_dir/

# Output: trained model + per-locus base-resolution predictions
```

Training cost: 1-3 GPU days for a single chromBPNet model; multiple GPUs for EnFormer.

## EnFormer Application

```python
from enformer_pytorch import Enformer

model = Enformer.from_hparams(dim_factor=32).from_pretrained('EleutherAI/enformer-official-rough')

# 196 kb input window
seq = torch.tensor(one_hot_encode(reference_seq))[None, ...]
predictions = model(seq)
# Predictions: per-bin signal across 5,313 ENCODE tracks
# Each variant effect = difference in target track prediction

# Variant effect at SNP
ref_pred = model(encode(ref_seq))[..., :, target_track_idx]
alt_pred = model(encode(alt_seq))[..., :, target_track_idx]
variant_effect = (alt_pred - ref_pred).mean()
```

EnFormer's 196 kb input (~100 kb effective receptive field) captures distal regulatory effects; useful when variant is far from TSS.

## Using JASPAR 2026 Deep Learning Models (Precomputed)

JASPAR 2026 (released 2025) added 1259 BPNet models trained on ENCODE TF ChIP-seq:

```python
from pyjaspar import jaspardb

jdb = jaspardb(release='JASPAR2026')

# Get the BPNet model for a specific TF
bpnet_model_info = jdb.fetch_matrix_by_collection('BPNET')
# Each entry has a downloadable model URL and training metadata

# Use a model for in silico mutagenesis on new variants
# (Models are typically Keras H5 or PyTorch state dicts)
```

This is the lowest-effort path for variant-effect prediction on canonical TFs (no training required).

## Per-Tool Failure Modes

### chromBPNet -- Bias model trained on wrong assay

**Trigger:** Using ATAC bias model on ChIP data.

**Mechanism:** ATAC bias model captures Tn5 sequence preferences; ChIP has different bias structure (sonication, fragmentation, antibody-driven).

**Symptom:** Variant effect predictions noisy; attribution scores dominated by Tn5 sequence preferences.

**Fix:** Train bias model on ChIP non-peak regions, not ATAC; or use chromBPNet's ChIP-specific bias correction.

### BPNet -- Trained on insufficient peaks

**Trigger:** Training BPNet on a TF with <5000 high-confidence peaks.

**Mechanism:** Base-resolution profile prediction needs many examples per motif context.

**Symptom:** Model accuracy <0.6 Spearman correlation between predicted and observed profiles.

**Fix:** Combine replicates; use more permissive peak threshold; or fall back to PWM-based motif analysis.

### TF-MoDISco -- Background sequences not representative

**Trigger:** Using random genomic sequences as background for attribution.

**Mechanism:** Genomic sequences include other TF binding sites; attribution conflates target TF with background TFs.

**Fix:** Use shuffled-input sequences as background (preserves dinucleotide); OR use peaks from a control ChIP (e.g., IgG) as background.

### In silico mutagenesis -- Variant outside training distribution

**Trigger:** Predicting effect of a variant in a sequence context the model never saw (e.g., new TF site arrangement).

**Mechanism:** Deep-learning models extrapolate poorly; counterfactual prediction for novel contexts is unreliable.

**Symptom:** Variant effect estimate has huge variance across model replicates.

**Fix:** Train ensemble of 5-10 models; use disagreement as uncertainty estimate; for high-stakes claims, validate experimentally.

### EnFormer -- Tissue-aggregated predictions

**Trigger:** Predicting variant effect for a specific cell line that has its own ChIP-seq.

**Mechanism:** EnFormer was trained on tissue-aggregated tracks; cell-line-specific resolution is limited.

**Fix:** For cell-line-specific predictions, fine-tune EnFormer on cell-line ChIP-seq OR use chromBPNet trained on cell-line data.

### Memory / GPU requirements

**Trigger:** Training chromBPNet or EnFormer on a single GPU with insufficient memory.

**Mechanism:** chromBPNet (~10 GB GPU memory), EnFormer (~16-24 GB), batch sizes / sequence lengths affect memory.

**Fix:** Reduce batch size; use gradient checkpointing; for EnFormer, use the smaller 5x architecture variant.

## Reconciliation with PWM-Based Analysis

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| TF-MoDISco motif matches JASPAR PWM | DL model recovered canonical motif | Confidence in DL model |
| TF-MoDISco motif doesn't match any PWM | Novel motif syntax OR DL model overfit | Validate with TOMTOM against larger DBs; check model ensemble agreement |
| chromBPNet predicts strong variant effect, JASPAR PWM scan does not | Soft syntax / cooperativity; DL captures more than PWM | DL prediction often more accurate; experimental validation ideal |
| EnFormer and chromBPNet disagree on variant effect | Different receptive fields capture different biology | Trust EnFormer for distal effects, chromBPNet for local |
| Variant effect ensemble disagrees within model | Training instability OR variant in extrapolation regime | Treat as low-confidence; do not publish without validation |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `chrombpnet` import fails | TensorFlow / Keras version mismatch | Use chrombpnet conda env with pinned tf 2.13 |
| `cuda out of memory` | Batch size too large | Reduce batch_size in training config |
| Predictions all zero | Bias model corrupted or wrong assay | Re-train bias on matched assay |
| TF-MoDISco produces empty motif list | FDR too strict or attribution noisy | Lower `target_seqlet_fdr` to 0.10; increase model training depth |
| EnFormer "shape mismatch" | Sequence not exactly 196608 bp (default) | Pad or truncate to expected length |
| Variant effect for nucleotide outside ACGT | Model only handles ACGT | Skip variants with N or ambiguous nucleotides |

## References

- Avsec Ž et al 2021 Nat Genet 53:354 (BPNet; ChIP-nexus base-resolution model)
- Pampari A et al 2025 Nat Genet (chromBPNet; bias-factorized base-resolution for ATAC and ChIP; consult current publication for exact volume/pages)
- Avsec Ž et al 2021 Nat Methods 18:1196-1203 (EnFormer; ~100 kb effective receptive field, 196 kb input window)
- Zhou J & Troyanskaya OG 2015 Nat Methods 12:931 (DeepSEA)
- Alipanahi B et al 2015 Nat Biotechnol 33:831 (DeepBind)
- Kelley DR et al 2016 Genome Res 26:990 (Basset)
- Shrikumar A et al 2020 bioRxiv (TF-MoDISco)
- Shrikumar A et al 2017 ICML (DeepLIFT)
- Lundberg SM & Lee SI 2017 NeurIPS (SHAP)
- JASPAR Project 2026 (deep-learning collection)
- Karbalayghareh A et al 2022 bioRxiv (deep-learning generalization across cell types)

## Related Skills

- chip-seq/peak-calling - Source peaks for training DL models
- chip-seq/motif-analysis - PWM-based analysis (complementary to DL)
- chip-seq/allele-specific-binding - Validate DL variant predictions against ASB data
- atac-seq/deep-learning-atac - ATAC-specific chromBPNet workflow
- atac-seq/footprinting - TOBIAS footprints as comparison
- causal-genomics/fine-mapping - Variant-level functional annotation
- machine-learning/biomarker-discovery - DL variant scores as features
- machine-learning/model-validation - Ensemble agreement, cross-validation
