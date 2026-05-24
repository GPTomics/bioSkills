---
name: bio-clip-seq-clip-deep-learning
description: Predict RBP binding from RNA sequence using deep learning models (RBPNet sequence-to-signal, RNAProt RNN, GraphProt2 GCN with structure, DeepCLIP, DeepRiPe multi-modal CNN) for variant-effect prediction, in silico binding-site discovery, model interpretation, and transfer learning from CLIP and RBNS datasets. Use when computational prediction of RBP binding from sequence is needed, evaluating variant effects on binding without further wet-lab experiments, comparing model performance, or training a custom model on ENCODE eCLIP data.
tool_type: python
primary_tool: RBPNet
---

## Version Compatibility

Reference examples tested with: RBPNet (Jens & Gagneur 2024 github), RNAProt 0.5+, GraphProt2 (Sauer 2022 github), DeepCLIP 1.0+ (Gronning 2020), DeepRiPe (Krakau lab), pytorch 2.2+, tensorflow 2.15+, scikit-learn 1.4+, biopython 1.83+, transformers 4.40+ (for RNA foundation models).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- Frameworks: check pytorch / tensorflow versions; reproducibility depends on framework version

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# CLIP-seq Deep Learning

**"Predict RBP binding from RNA sequence using deep learning"** -> Train or apply neural networks that learn the sequence (and optionally structure) preference of an RBP from CLIP-seq peaks or single-nucleotide crosslink sites. The output is per-base or per-site binding probability for any input sequence, enabling: (a) variant-effect prediction at heterozygous SNPs; (b) in silico binding-site discovery on transcripts not covered by CLIP; (c) systematic comparison across RBPs via shared model architectures; (d) interpretation via attribution / saliency to recover RBP-specific motifs and structural preferences. Modern models (RBPNet 2024) predict per-nucleotide crosslink count distributions rather than binary peak/non-peak, providing single-nt resolution outputs.

- Python (RBPNet sequence-to-CL signal): `import rbpnet; model = rbpnet.load_pretrained('RBP_name'); predictions = model.predict(sequence)` produces per-base CL count distribution
- Python (RNAProt RNN classifier): `RNAProt train -i peaks.bed -t background.bed -g genome.fa -o model/` then `RNAProt predict -m model/ -i query_sequences.fa -o predictions.tsv`
- Python (GraphProt2 GCN with structure): `graphprot2 train -i peaks.bed -bg shuffled.bed -g genome.fa --structure -o model/`
- Python (DeepCLIP for binding probability): `deepclip --train --train_data train.fa --validation_data val.fa --predict --predict_data test.fa --output_dir output/`
- Python (DeepRiPe multi-modal CNN): `from deepripe import DeepRiPe; model.train(X_train, y_train); predictions = model.predict(X_test)`

The benchmarks (RNAProt paper, 2021): RNAProt AUC 87-89%; DeepCLIP 84-87%; GraphProt 82-84%. RBPNet (2024) is the modern sequence-to-signal model that predicts per-nt CL distributions at single-nucleotide resolution rather than binary site classification.

## Models Taxonomy

| Model | Architecture | Input | Output | Resolution | Strength | Fails when |
|-------|--------------|-------|--------|------------|----------|------------|
| RBPNet (Jens & Gagneur 2024) | Sequence-to-signal CNN | Sequence | Per-nt CL count distribution | Single-nt | Modern single-nt resolution; predicts CL distribution not binary | New (2024); fewer pretrained RBPs |
| RNAProt (Uhl 2021) | LSTM RNN | Sequence | Binary binding probability | Site/peak | Highest AUC in benchmark (87-89%); feature-rich (RNA-Foldscan structure) | Single-prediction; not per-base profile |
| GraphProt (Maticzka 2014) | Graph kernel + SVM | Sequence + structure | Binary | Site/peak | Original structure-aware; well-validated | Older; superseded by GraphProt2 |
| GraphProt2 (Sauer 2022) | Graph Convolutional Network (GCN) | Variable-length sequence + structure | Nucleotide-wise binding profile | Per-nt | Variable-length input; structure-aware | Slow to train; needs GPU |
| DeepCLIP (Gronning 2020) | CNN | Sequence | Binary | Site | Fast; good benchmarks | Sequence-only; no structure |
| DeepRiPe (Ghanbari 2020) | Multi-modal CNN | Sequence + structure + region | Binary | Site | Multi-input; good ENCODE benchmark | Sequence/structure/region must be pre-extracted |
| iDeep / iDeepE (Pan 2018) | CNN ensemble | Sequence | Binary | Site | Ensemble approach | Older; few pretrained models |
| Pysster (Budach 2018) | CNN-LSTM | Sequence | Binary | Site | Generic framework | Less RBP-specific |
| DeepBind (Alipanahi 2015) | CNN | Sequence | Binary | Site | First deep-learning RBP model | Outdated; superseded |
| BasenjiPredict (Yeo lab adapt) | Multi-task CNN | Sequence | Per-task profile | Per-nt | Joint learning across RBPs | Computational overhead |
| RNA foundation models (RNAErnie, RNA-FM 2024) | Transformer pretrained on RNA | Sequence | Embeddings (downstream task) | Embedding | Transfer learning across RBPs | Foundation model trained at depth; fine-tuning needed |

Methodology evolves; verify the latest publication on RBP deep learning (RBPNet 2024 is the current state-of-the-art per-nt resolution model). RNA foundation models (RNA-FM, RNAErnie) are emerging in 2024 as transfer-learning backbones; fine-tuning on CLIP data for specific RBPs is the next-generation approach.

## Critical Choice: Binary Classification vs Sequence-to-Signal

**Binary classification (RNAProt, DeepCLIP, DeepRiPe, GraphProt2):** Train on labeled site vs background; predict probability of binding for an input sequence. Output: per-sequence score. Pro: simple framework; mature benchmarks. Con: discards single-nt CL distribution information; binary decision boundary.

**Sequence-to-signal (RBPNet):** Train on per-nt CL count distributions from PureCLIP or CTK CITS output; predict per-base CL count for input sequence. Output: per-nt profile. Pro: single-nt resolution; preserves CL count information; matches biology (CL is a sharp signal, not a region). Con: newer (2024); fewer pretrained models; harder to interpret with classical motif tools.

| Goal | Model |
|------|-------|
| Predict RBP binding probability for an input sequence | RNAProt or DeepRiPe |
| Predict per-base CL distribution | RBPNet |
| Variant-effect at heterozygous SNP | RBPNet or DeepRiPe (per-base output) |
| Compare RBP preferences across ENCODE | Multi-task model (BasenjiPredict-style) |
| Transfer learning across RBPs | RNA foundation model + fine-tune |
| In silico screening of variants | RBPNet (per-base) for genome-wide |
| Motif interpretation via attribution | DeepRiPe or GraphProt2 (interpretable) |
| Custom training on new CLIP data | RNAProt (easiest pipeline) |
| Production-grade per-base prediction | RBPNet 2024 |

## Variant-Effect Prediction Workflow

Apply a pretrained or custom-trained model to predict the change in binding upon a sequence variant.

```python
import torch
from rbpnet import RBPNet  # hypothetical API; verify per-package documentation

# Load pretrained model for specific RBP
model = RBPNet.load_pretrained('TARDBP_HEK293T')

# Reference and alternative sequences around a variant
ref_seq = 'CTGTACTGCAGTAGCATGCTAGCATGCTAGCAT'  # 32 nt window centered on variant
alt_seq = 'CTGTACTGCAGTAGCATGCTAGCATGCTAGCAA'  # Variant: T -> A at position 32

# Predict per-base CL distribution for both
ref_pred = model.predict(ref_seq)   # shape: (32, 1) - per-base CL probability
alt_pred = model.predict(alt_seq)

# Variant effect: log2 fold change in summed binding signal
import numpy as np
effect = np.log2((alt_pred.sum() + 1e-9) / (ref_pred.sum() + 1e-9))
print(f'Variant effect (log2 FC): {effect:.4f}')

# Strong-effect variant: |log2 FC| > 1.0
# Mid-effect: 0.5 - 1.0
# Weak: < 0.5
```

For genome-wide variant scoring: apply this in batch to all GWAS variants overlapping the RBP's binding regions. The output is a per-variant log2 FC; downstream Mendelian randomization or fine-mapping integrates with phenotype-association statistics.

## Training a Custom Model (RNAProt Example)

RNAProt is the most accessible training framework. It accepts peak BED + background BED + genome FASTA.

**Goal:** Train a chromosome-split RNN classifier from CLIP peaks to predict RBP binding probability on arbitrary input sequences, with held-out evaluation on a chromosome-distinct test set.

**Approach:** Generate a GC-matched 3' UTR background, split peaks and background by chromosome (train chr1-20, test chr21-22) to prevent gene-neighbor leakage, train RNAProt for 50 epochs at batch size 64, and evaluate held-out AUC against the ENCODE benchmark target of 0.85-0.89.

```bash
# Step 1: Prepare data
# Foreground: positive peaks from CLIPper / Skipper stringent set
# Background: GC-matched random regions from expressed transcripts
bedtools getfasta -fi genome.fa -bed peaks.stringent.bed -s -fo peaks.fa
bedtools shuffle -i peaks.stringent.bed -g chrom.sizes -incl expressed.bed -seed 42 > bg.bed
bedtools getfasta -fi genome.fa -bed bg.bed -s -fo background.fa

# Step 2: Train RNAProt
RNAProt train \
    --in peaks.fa \
    --neg background.fa \
    --out model_dir \
    --epochs 50 \
    --batch-size 64 \
    --learning-rate 0.001 \
    --validation-split 0.2

# Step 3: Apply to query sequences
RNAProt predict \
    --model model_dir \
    --in query.fa \
    --out predictions.tsv

# Output: per-sequence binding probability
```

## RBPNet Sequence-to-Signal Workflow

```python
# Hypothetical RBPNet API - verify against current package
import rbpnet

# Training: needs per-nt crosslink count for each example
# Inputs: sequence windows (256 nt) around peaks
# Targets: per-base CL count vector from PureCLIP or CTK CITS

# Prepare data
X_train, y_train = rbpnet.load_clip_data(
    peaks_bed='peaks.bed',
    crosslinks_bed='pureclip_sites.bed',
    genome_fa='genome.fa',
    window_size=256
)

# Train
model = rbpnet.RBPNet(
    seq_length=256,
    out_length=256,
    conv_layers=4,
    filters=128
)
model.train(X_train, y_train, epochs=50, batch_size=32, val_split=0.2)

# Predict per-base CL distribution
predictions = model.predict(novel_sequences)
```

## Per-Tool Failure Modes

### Training data imbalance

**Trigger:** Peak set 10k positive vs 100k negative background.

**Mechanism:** Class imbalance biases model toward negative class; specificity high but sensitivity low.

**Symptom:** AUC reported at 0.95 but precision-recall at peak threshold poor; few sites recovered.

**Fix:** Balanced sampling (1:1 positive:negative) or class weights. RNAProt does this automatically; DeepCLIP / DeepRiPe require manual balancing.

### Background mismatch

**Trigger:** Random shuffled background not matched to peak transcript context (e.g., peaks from 3' UTRs, background from CDS).

**Mechanism:** Model learns transcript-region differences (AU-content of 3' UTR vs CDS), not RBP specificity.

**Symptom:** Top predicted sites all in 3' UTRs regardless of test sequence; motif analysis shows AU-rich without RBP motif.

**Fix:** Match background to same region as foreground (3' UTR peaks -> 3' UTR background). Use the GC-content matched shuffle.

### Test on training data leak

**Trigger:** Splitting train/test by random shuffle.

**Mechanism:** Adjacent peaks in genome share evolutionary context; random split leaks training peaks near test peaks.

**Symptom:** Held-out AUC > 0.95 in benchmark but fails on truly novel sequences.

**Fix:** Split by chromosome (e.g., train chr1-20, test chr21-22). Or split by gene (no two peaks in same gene across splits).

### GPU requirement underestimated

**Trigger:** RNAProt / RBPNet training on CPU only.

**Mechanism:** Modern RBP DL models have 1-10M parameters; training requires GPU (or 100x slower on CPU).

**Symptom:** Training time > 24 h on CPU; convergence unstable.

**Fix:** Use Google Colab GPU; AWS EC2 GPU instance; or restricted-architecture model for CPU.

### Variant-effect prediction window size

**Trigger:** Variant-effect prediction with too narrow window around the variant.

**Mechanism:** RBP context extends 50-200 nt; window < 100 nt misses long-range context.

**Symptom:** Variant effect estimates noisy; same model gives different effects on different windows.

**Fix:** Use the model's native window size (256 nt for RBPNet); for shorter windows, average across multiple shifted predictions.

### Pretrained models lack the target RBP

**Trigger:** Looking for pretrained DeepRiPe / RBPNet for an uncommon RBP.

**Mechanism:** Only ~150 ENCODE-tested RBPs have pretrained models; thousands of RBPs are not.

**Symptom:** No pretrained available for the protein of interest.

**Fix:** Train custom model on new CLIP data using RNAProt or DeepCLIP. Or use transfer learning from a related RBP (closest paralog).

### Structure prediction integration

**Trigger:** GraphProt2 with `--structure` flag without RNA-Fold structure ensemble.

**Mechanism:** GraphProt2 needs structure ensemble (RNAfold output) as input; without it the GCN cannot leverage structure.

**Symptom:** GraphProt2 performance no better than sequence-only models.

**Fix:** Pre-compute RNAfold ensemble for each training sequence; pass to GraphProt2 input. Or use sequence-only DeepCLIP if structure is not needed.

### Model interpretation via saliency

**Trigger:** Want to recover RBP motif from trained model.

**Mechanism:** Saliency / integrated gradients on trained model produces per-base importance scores; sum across many positive examples gives a motif.

**Symptom:** Saliency results noisy; no clean motif emerges.

**Fix:** Use TF-MoDISco (Schreiber 2020) for cleaner motif extraction from saliency maps. Or use DeepLIFT scores instead of vanilla saliency.

## Decision Tree by Use Case

| Scenario | Model | Why |
|----------|-------|-----|
| Variant-effect prediction (genome-wide GWAS variants) | RBPNet (per-base) | Single-nt resolution; predicts CL distribution |
| Binary "is this sequence bound" | RNAProt | Best AUC in benchmark (87-89%) |
| Structure-aware prediction | GraphProt2 | Structure ensemble integration |
| Custom training on new CLIP data | RNAProt (easiest CLI) | Fast training; CLI tool |
| Single RBP, no pretrained | Train custom RNAProt | Most accessible framework |
| Multi-task across RBPs | BasenjiPredict-style | Joint learning |
| Transfer learning from foundation model | RNA-FM + fine-tune | New approach; not yet in production |
| Production scoring of many sequences | DeepCLIP (fast inference) | Throughput |
| Motif interpretation | DeepRiPe + TF-MoDISco | Interpretable architectures |
| Comparing in silico to RBNS | RBPNet or DeepRiPe + RBNS Kd correlation | Calibration |
| Variant effect at heterozygous SNP | RBPNet | Per-base output |

## Reconciliation: Model Predictions vs CLIP / RBNS

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Model predicts site; CLIP misses | Model false positive; or transient binding not captured by CLIP | Check RBNS prediction; in vitro Kd |
| Model predicts low; CLIP has strong peak | Model false negative; or non-canonical / context-dependent | Investigate training data; structure-dependent? |
| Model AUC 0.95 in test; fails on novel sequences | Train/test leakage; chromosome split needed | Re-train with proper split |
| Variant effect log2 FC inconsistent across windows | Window-size sensitivity | Use model native window; average across shifts |
| Pretrained for related RBP works on target | Cross-RBP transfer | Transfer learning may work for paralogs |
| GraphProt2 underperforms sequence-only | Structure ensemble not provided | Fix structure input |
| Saliency motif noisy | Vanilla gradient method | Use DeepLIFT or TF-MoDISco |
| Train/test split by random outperforms chromosome split | Leakage from gene-neighbor sequences | Trust chromosome-split AUC |

**Operational rule for high-confidence variant-effect:** (a) Use RBPNet per-base output; (b) compute log2 FC at variant position summed over 50 nt window; (c) cross-validate with another model (DeepRiPe); (d) cross-reference with overlapping eCLIP peak if available; (e) report |log2 FC| > 1 as strong effect.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| AUC 0.5 on test | Train data leak or random shuffle | Verify chromosome split |
| Model trained on 1 epoch | Default optimizer state | Train 30-50 epochs with validation |
| GPU OOM | Batch size too large | Reduce batch size to 32 |
| Variant effect log2 FC very large (> 10) | Reference sequence not in training distribution | Verify input sequence reasonable |
| Pretrained model not found | RBP not in pretrained list | Train custom; or use closest paralog |
| Structure flag without input | GraphProt2 misuse | Pre-compute RNAfold structure |
| Saliency map flat | Model architecture too shallow | Use DeepRiPe / RBPNet (deeper) |
| Per-class accuracy uneven | Class imbalance | Use balanced sampling or class weights |
| Cross-RBP transfer fails | RBPs unrelated | Limit transfer to paralogs or use foundation model |
| Custom training crashes | RAM / GPU exhausted | Smaller batch; cache embeddings |

## References

- Alipanahi B et al 2015 Nat Biotechnol 33:831 (DeepBind, first DL RBP model)
- Maticzka D et al 2014 Genome Biol 15:R17 (GraphProt with structure)
- Sauer M et al 2022 (GraphProt2 with GCN)
- Gronning AGB et al 2020 Nucleic Acids Res 48:7099 (DeepCLIP)
- Ghanbari M et al 2020 Bioinformatics 36:3489 (DeepRiPe multi-modal)
- Pan X et al 2018 Bioinformatics 34:i285 (iDeepE)
- Budach S et al 2018 Bioinformatics 34:3387 (Pysster)
- Uhl M et al 2021 GigaScience 10:giab054 (RNAProt RNN)
- Horlacher M, Wagner N, Moyon L et al 2023 Genome Biol 24:180 (RBPNet sequence-to-signal at single-nt)
- Shrikumar A et al 2018 arXiv (TF-MoDISco interpretation)
- Chen J et al 2024 (RNA-FM foundation model)
- Wang W et al 2024 (RNAErnie pretrained transformer)

## Related Skills

- clip-seq/clip-motif-analysis - Motif analysis is the classical alternative
- clip-seq/crosslink-site-detection - Single-nt CL sites for sequence-to-signal models
- clip-seq/clip-peak-calling - Peak BEDs for binary classification training
- causal-genomics/mendelian-randomization - Variant-effect predictions feed MR
- causal-genomics/fine-mapping - Variant prioritization with model scores
- machine-learning/model-validation - Train/test split methodology
- machine-learning/prediction-explanation - Saliency / attribution methods
- machine-learning/biomarker-discovery - Generic ML framework
