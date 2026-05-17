---
name: bio-atac-seq-deep-learning-atac
description: Sequence-based deep learning for ATAC-seq using chromBPNet, BPNet, scBasset, or EnFormer. Use when correcting Tn5 bias with neural networks beyond k-mer models, predicting per-base accessibility profiles, scoring in silico variant effects at GWAS or rare-variant SNPs, discovering motifs via DeepLIFT/TF-MoDISco from a trained model, or generating cell-type-specific accessibility predictions for unobserved cell states.
tool_type: python
primary_tool: chrombpnet
---

## Version Compatibility

Reference examples tested with: chrombpnet 0.1.7+, bpnet-lite 0.6+ (DOI 10.5281/zenodo.7011327), scBasset 0.1.0+ (basenji2 fork), tangermeme 0.1+, tfmodisco-lite 2.2+, DeepLIFT 0.6+, captum 0.7+, tensorflow 2.13+, pytorch 2.1+, kipoi 0.8+.

Verify before use:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt rather than retrying. Deep-learning tooling evolves rapidly; method papers post 2023 may have superseded reference implementations.

# Sequence-Based Deep Learning for ATAC-seq

**"Score the effect of a GWAS SNP on chromatin accessibility"** -> Train (or use pre-trained) sequence-to-accessibility CNNs that take 1-5 kb DNA windows and predict per-base Tn5 cleavage profiles. Outputs include: bias-corrected accessibility, single-base mutation effect predictions, and DeepLIFT contribution scores convertible to motifs via TF-MoDISco.

- CLI: `chrombpnet pipeline --bigwig signal.bw --bigwig-bias bias.bw ...`
- Python: `bpnet-lite` for custom architectures; `tangermeme` for fast scoring
- Python (single-cell): `scBasset` for per-cell sequence-based predictions
- Python (long-context): EnFormer pre-trained models via Kipoi

Sequence models are NOT a replacement for MACS+TOBIAS at every step. They excel at three specific tasks where classical pipelines struggle: (1) Tn5 bias correction in low-complexity sequence contexts, (2) variant effect prediction in non-genic regions, (3) cell-type-specific motif discovery beyond what JASPAR provides.

## Algorithmic Taxonomy

| Tool | Architecture | Training | Output | Strength | Fails when |
|------|-------------|----------|--------|----------|------------|
| chromBPNet (Pampari 2024) | Two-track CNN: bias model + accessibility model; bias trained on naked-DNA control or k-mer baseline, accessibility trained on chromatin signal | Per-cell-type, paired bias track | Bias-corrected per-base profile + total counts | Best-in-class bias correction; established in Kundaje lab pipelines | Requires GPU, ~24h training per cell type; needs >= 50M reads |
| BPNet (Avsec 2021) | Original counts + profile dual-head CNN | TF ChIP-seq or ATAC | Per-base profile prediction | Foundational; widely cited; bpnet-lite reimpl maintained | Less polished than chromBPNet for ATAC; bias correction needs separate model |
| scBasset (Yuan & Kelley 2022) | Basenji2-derived CNN, per-cell projection layer | Single-cell ATAC | Per-cell sequence-derived peak score | First sequence model that predicts per-cell accessibility; outperforms chromVAR for cluster discrimination | Fixed architecture, hard to extend; benchmarks evolving |
| EnFormer (Avsec 2021 Nat Methods) | Long-context Transformer (196 kb input) | Reference epigenome (DNase + histones + CAGE) | Per-bin epigenome prediction | Best for distal regulation modeling; pre-trained available | Pre-trained models cell-line specific; finetuning on custom data is expensive |
| Borzoi (Linder 2025) | EnFormer extension trained on RNA + ATAC | Multi-tissue paired data | Sequence -> RNA + chromatin | Recent SOTA for variant effect on RNA via ATAC linkage | Newer; benchmarks still emerging |
| DeepATAC / Basset (legacy) | Earlier CNN architectures | -- | Binary peak prediction | Historical context; cited in older literature | Superseded by chromBPNet + EnFormer; do not use for new work |
| tangermeme | Inference-only fast wrapper | Use any saved model | Marginal scoring of variants | Speeds up variant effect prediction 100x; works with chromBPNet/BPNet outputs | Inference only; cannot train |

Methodology evolves; verify against current Kundaje lab pipelines (chrombpnet GitHub), Greenleaf lab (scBasset), and Avsec / Linder publications before locking pipelines.

## When Deep Learning Helps vs When Classical Pipelines Suffice

| Task | Classical | Deep Learning |
|------|-----------|---------------|
| Peak calling | MACS3 / Genrich (sufficient) | chromBPNet (overkill unless variant downstream) |
| Tn5 bias correction at TF motifs | TOBIAS ATACorrect (good) | chromBPNet (better at hard cases: low-complexity flanks, deep TF footprints) |
| Differential accessibility | DiffBind / DESeq2 (sufficient) | -- (no clear DL advantage) |
| GWAS variant effect prediction at causal SNPs | Limited (overlap heuristics) | chromBPNet / EnFormer (essential) |
| Motif discovery from de novo data | MEME / HOMER (good) | chromBPNet + TF-MoDISco (better; finds composite + cooperative motifs) |
| Per-cell TF activity | chromVAR (sufficient at the cluster level) | scBasset (better at fine-grained cell states) |
| Cross-cell-type accessibility prediction | -- | EnFormer / Borzoi (only option) |
| Predicting cell-type-specific enhancer activity from sequence | -- | chromBPNet / EnFormer (essential) |

For most standard ATAC analysis, classical pipelines remain primary. Deep learning enters when (a) variant interpretation is the goal, (b) cell-type prediction is needed beyond observed data, or (c) bias correction quality is paramount (low-input, FFPE, transcription factors with weak motifs).

## Per-Tool Failure Modes

### chromBPNet -- Bias model mismatch

**Trigger:** Training the bias model on a dataset different from the accessibility dataset (e.g. K562 bias model used on primary T cells).

**Mechanism:** chromBPNet's bias model captures sequence-specific Tn5 preference, which is mostly cell-type-invariant BUT contributions of chromatin context at cuts can vary. Cross-celltype bias models work but with degraded performance.

**Symptom:** Predicted footprints look correct at known TFs (CTCF) but fail on cell-type-specific regulators.

**Fix:** Train a per-cell-type bias model from naked-DNA control if available, OR use the chromBPNet authors' pre-trained k562 / GM12878 / HepG2 bias as a fallback (acknowledged degradation).

### chromBPNet -- Insufficient training data

**Trigger:** Training on < 50M deduplicated nuclear reads, or < 30k peaks.

**Mechanism:** CNN training needs enough peaks for stable gradient updates and enough background regions for the dual-task loss.

**Fix:** Pool replicates before training; reduce model capacity (`--num-filters`); use pre-trained model on closest cell type and skip retraining.

### BPNet / chromBPNet -- DeepLIFT vs Integrated Gradients confusion

**Trigger:** Computing per-base contributions for motif discovery.

**Mechanism:** DeepLIFT (RevealCancel rule) and Integrated Gradients (50 baseline samples) give different attribution patterns. DeepLIFT preserves additivity; IG is stochastic.

**Fix:** Use DeepLIFT rescale-rule (chromBPNet default) for TF-MoDISco. IG only when DeepLIFT fails on saturating activations. Document the choice.

### scBasset -- Cell projection layer instability

**Trigger:** Few cells per cluster; sparse training data.

**Mechanism:** scBasset learns a per-cell projection vector; with < 100 cells per cluster the projection is noisy.

**Fix:** Aggregate cells to clusters before training, OR use chromBPNet trained on pseudobulks per cluster instead.

### EnFormer -- Pre-trained models lack target cell type

**Trigger:** Using EnFormer for variant effects in a cell type not in its training set (e.g. GTEx tissues are covered; novel primary cell types are not).

**Mechanism:** EnFormer's outputs are per-track predictions; if the target cell type wasn't trained, the agent can use a similar track as proxy but accuracy degrades.

**Fix:** Use a similar tissue track as proxy (HepG2 for liver biology; GM12878 for B-cell-like) OR fine-tune EnFormer on custom data (expensive). Document the proxy.

### tangermeme -- Marginal vs in silico mutagenesis confusion

**Trigger:** Asking for a "variant effect score" without specifying the formula.

**Mechanism:** Marginal effects = ref vs alt at the SNP only. ISM = saturation across all positions in the window (every base mutated). Different magnitudes; different questions.

**Fix:** Define which calculation. For GWAS variant prediction, use marginal at the SNP (matches phenotype-genotype coupling). For motif discovery, use ISM.

## Decision Tree by Goal

| Goal | Recommended approach |
|------|---------------------|
| Score 100 GWAS SNPs for chromatin effects | Pre-trained chromBPNet model on closest cell type; tangermeme for fast scoring |
| Score 1 lead SNP at high resolution | chromBPNet + tangermeme + ISM saturation map |
| Identify TF binding motifs from a new cell type's ATAC | chromBPNet train + DeepLIFT contributions + TF-MoDISco-lite |
| Predict accessibility in a cell type not in training | EnFormer pre-trained (best for ENCODE cell types) or scBasset for sc state interpolation |
| Bias-correct a low-input ATAC library before footprinting | chromBPNet bias model output as `--bias` to TOBIAS or directly use chromBPNet corrected track |
| Cell-type-specific enhancer prediction | chromBPNet trained on each cell type; per-cell-type ISM at candidate loci |
| Replace TOBIAS bias correction | chromBPNet corrected bigWig as input to TOBIAS ScoreBigwig; skip ATACorrect |

## chromBPNet Standard Pipeline

**Goal:** Train a bias-corrected CNN that predicts per-base accessibility from sequence, then score variants with it.

**Approach:** Generate train/valid/test chromosome splits, train the Tn5 bias model on background regions, train the accessibility model with bias correction, then run the standalone variant-scorer repo to predict ref-vs-alt effects at SNPs.

```bash
# 1. Generate train / valid / test chromosome splits (output is a JSON file with chrom assignments)
chrombpnet prep splits \
    -c hg38.chrom.sizes \
    -tcr chr1 chr3 chr6 \
    -vcr chr8 chr20 \
    -op splits/fold_0
# Train chromosomes are auto-inferred (whatever is not in -tcr/-vcr). The `-tecr` flag does NOT exist.

# 2. Train bias model from background regions
chrombpnet bias pipeline \
    -ibam atac.bam \
    -d ATAC \
    -g hg38.fa \
    -c hg38.chrom.sizes \
    -p peaks.narrowPeak \
    -n nonpeaks.bed \
    -fl splits/fold_0.json \
    -b 0.5 \
    -o bias_model/

# 3. Train accessibility model with bias correction
chrombpnet pipeline \
    -ibam atac.bam \
    -d ATAC \
    -g hg38.fa \
    -c hg38.chrom.sizes \
    -p peaks.narrowPeak \
    -n nonpeaks.bed \
    -fl splits/fold_0.json \
    -b bias_model/bias.h5 \
    -o output/

# 4. Variant effect prediction at GWAS SNPs uses the SEPARATE kundajelab/variant-scorer repo
# (the `chrombpnet snp_score` subcommand is commented out in current chrombpnet/parsers.py)
git clone https://github.com/kundajelab/variant-scorer
python variant-scorer/src/variant_scoring.py \
    --model output/chrombpnet_no_bias.h5 \
    --list variants.tsv \
    --genome hg38.fa \
    --chrom_sizes hg38.chrom.sizes \
    --out_prefix variants_predicted
# Output: variants_predicted.variant_scores.tsv with per-SNP log2FC magnitudes
```

`-b 0.5` is the bias scaling factor; chromBPNet docs recommend 0.5-1.0 depending on enrichment. For variant scoring, use the standalone `kundajelab/variant-scorer` companion repo, NOT a chrombpnet subcommand. Verify exact flag names with `python variant_scoring.py --help` because the API evolves.

## DeepLIFT + TF-MoDISco for Motif Discovery

The maintained version is `tfmodisco-lite` (jmschrei/tfmodisco-lite, `pip install modisco-lite`), which exposes a CLI rather than the deprecated v1 `TfModiscoWorkflow` Python API. The original `kundajelab/tfmodisco` package (with `tfmodisco.tfmodisco_workflow.workflow.TfModiscoWorkflow`) is unmaintained and incompatible with `modisco-lite`.

**Goal:** Discover de novo motifs from a trained chromBPNet model using per-base contribution scores.

**Approach:** Extract one-hot sequences and DeepLIFT/SHAP contributions from chromBPNet, run modisco-lite to cluster seqlets into motif patterns, then generate an annotated HTML report matched against a known motif database.

```bash
# Generate one-hot sequence and SHAP / DeepLIFT contribution score arrays from chromBPNet
# (chromBPNet `chrombpnet contribs_bw` writes hypothetical contributions; convert to numpy via shap_to_modisco)

# Run TF-MoDISco-lite via its CLI
modisco motifs \
    -s ohe.npz \
    -a shap.npz \
    -n 2000 \
    -w 500 \
    -o modisco_results.h5

# Generate HTML report with discovered motifs matched to known databases
modisco report \
    -i modisco_results.h5 \
    -o modisco_report/ \
    -m motifs_meme.txt \
    -s modisco_report/
```

`-n 2000` caps seqlets per metacluster; `-w 500` is the sliding-window length. `motifs_meme.txt` (e.g. JASPAR or HOCOMOCO MEME-format) lets `modisco report` annotate clusters against known motifs.

## In Silico Variant Effect Prediction

**Goal:** Score the effect of SNPs on predicted accessibility using a trained chromBPNet model.

**Approach:** Load the bias-free chromBPNet model, build ref and alt one-hot windows centered on each SNP, run tangermeme's substitution_effect to get paired predictions, then compute log2(alt/ref) per variant.

```python
import numpy as np

# tangermeme's variant-effect API: substitution_effect for SNPs, marginalize for motif insertions
from tangermeme.variant_effect import substitution_effect
from tangermeme.predict import predict

# Load pre-trained chromBPNet model (saved as Keras .h5 or PyTorch state_dict).
# chromBPNet wraps Keras; load with tensorflow.keras.models.load_model and wrap for tangermeme.
# `load_chrombpnet_model` below is pseudocode -- substitute the actual loader for the installed version
# (e.g. tf.keras.models.load_model + tangermeme.io.adapter, or torch.load for PyTorch checkpoints).
model = load_chrombpnet_model('output/chrombpnet_no_bias.h5')

# substitution_effect: per-SNP ref vs alt prediction across a sequence window
# X shape (N, 4, L); substitutions is a sparse-COO tensor of shape (-1, 3) where each row is
# [example_idx, position, new_base_idx] (new_base_idx 0-3 for ACGT)
y_ref, y_alt = substitution_effect(model, X, substitutions)
log2fc = np.log2(y_alt.sum(axis=-1) / y_ref.sum(axis=-1))
```

For motif marginalization (testing a candidate motif's effect by inserting it into background sequences), use `tangermeme.marginalize.marginalize(model, X, motif)`. The `motif` argument is a one-hot tensor of shape `(-1, 4, motif_length)`; convert string motifs via `tangermeme.utils.one_hot_encode`. Verify the exact signatures with `help(tangermeme.marginalize.marginalize)` because tangermeme is actively developed; `marginal_predict` is NOT a real function name.

`log2fc` magnitudes are unitless; |log2fc| > 1 typical for strong-effect SNPs in regulatory regions.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| chromBPNet predicts strong effect; MACS does not call peak | Sequence model captures latent regulatory potential | Trust chromBPNet for variant effect; not for peak calling |
| EnFormer prediction differs from chromBPNet at same locus | Different context windows (196 kb vs 1-2 kb); different cell types | Both can be correct at different scales; report both with their context size |
| TF-MoDISco motifs differ from JASPAR | Different methodology (sequence-based vs ChIP-validated) | TF-MoDISco can find composites and cooperative; check JASPAR for confirmation |
| chromBPNet bias correction differs from TOBIAS ATACorrect | Different bias models (CNN vs k-mer) | chromBPNet is more accurate but slower; TOBIAS still publishable for standard use |

**Operational rule:** For high-confidence variant prediction, agree across two approaches: chromBPNet + EnFormer (or Borzoi). Single-tool calls should be reported as exploratory. For motif discovery, validate TF-MoDISco hits against JASPAR/HOCOMOCO before publication.

## GPU and Compute Considerations

| Task | Hardware | Wall time |
|------|---------|-----------|
| chromBPNet training (per cell type) | 1 A100 GPU, 80 GB RAM | ~24 h |
| chromBPNet inference at 1M variants | 1 A100 | ~4 h |
| EnFormer pre-trained inference | 1 V100+ | ~30 min for 100k variants |
| Borzoi training | 1 A100, ~250 GB RAM | ~7 days |
| scBasset training (10k cells) | 1 V100, 32 GB RAM | ~12 h |
| TF-MoDISco on 1M peaks | CPU 32 cores | ~6 h |

For most labs without sustained GPU access: use pre-trained chromBPNet/EnFormer models for inference; only train custom models when the cell type is not in the public model zoo (encodeproject.org/atac-seq pre-trained chromBPNet).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| chromBPNet `bias.h5` missing | Bias model training failed silently | Re-run `chrombpnet bias pipeline` with verbose; check input BAM size |
| Out of memory during training | Default batch size too large for GPU | `--batch-size 64` or smaller; reduce `--num-filters` |
| Predicted profile is constant | Model collapsed (training too short) | Increase epochs; verify input peaks are non-empty |
| TF-MoDISco produces too many small clusters | `target_seqlet_fdr` too loose | Tighten to 0.01; or increase `flank_size` |
| EnFormer prediction has wrong shape | Pre-trained model expects 196 kb input | Pad input to exactly 196,608 bp |
| Variant effect predictions cluster near zero | SNP outside model's effective window | Predict on window-centered sequences (variant at the center) |
| chromBPNet model not converging | Peaks file contains chrM or blacklist | Pre-filter; chromBPNet does not auto-filter |
| scBasset training crashes on Apple Silicon | TensorFlow Metal incompatible with operations | Use CPU mode or run on Linux GPU |

## References

- Pampari A et al 2024 bioRxiv (chromBPNet; Tn5 bias correction with deep learning)
- Avsec Z et al 2021 Nat Methods 18:1196 (BPNet; foundational sequence-to-profile)
- Avsec Z et al 2021 Nat Methods 18:1224 (EnFormer; long-context Transformer)
- Linder J et al 2025 (Borzoi; multi-tissue sequence-to-RNA+chromatin)
- Yuan H & Kelley DR 2022 Nat Methods 19:1088 (scBasset)
- Shrikumar A et al 2017 ICML (DeepLIFT)
- Schreiber J et al 2024 (tangermeme; fast inference utilities)
- Shrikumar A et al 2018 bioRxiv (TF-MoDISco)
- Kelley DR 2020 Genome Res 30:1133 (Basenji2; precursor)

## Related Skills

- atac-seq/atac-peak-calling - Classical peak calling input
- atac-seq/footprinting - Use chromBPNet bias correction as TOBIAS alternative
- atac-seq/motif-deviation - chromVAR vs scBasset for per-cell motif activity
- atac-seq/single-cell-atac - scBasset integration with sc workflow
- atac-seq/enhancer-gene-linking - Variant effect feeds enhancer scoring
- atac-seq/allele-specific-accessibility - DL-predicted variant effects vs observed allelic imbalance
- causal-genomics/fine-mapping - Downstream use of variant effect scores
- machine-learning/biomarker-discovery - General ML patterns
- gene-regulatory-networks/scenic-regulons - Combine motif discovery with TF networks
