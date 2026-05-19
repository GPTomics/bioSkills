# Differential CLIP - Usage Guide

## Overview

Identify regions with differential RBP binding across conditions (KD vs control, treatment vs vehicle). Three scales: window-level (DEWSeq, Flipper), peak-level (DESeq2/edgeR/limma-voom), or crosslink-level (rare). The canonical CLIP statistical design uses the interaction term `~ type + condition + type:condition` testing whether the IP/SMInput ratio shifts across conditions - NOT testing condition on IP counts alone (which confounds binding with expression). DEWSeq + htseq-clip is the EMBL-Hentze framework; Flipper is the Skipper-companion modern alternative.

## Prerequisites

```r
BiocManager::install(c('DEWSeq','DESeq2','edgeR','limma'))
```
```bash
pip install htseq-clip
# Flipper (with Skipper)
git clone https://github.com/algaebrown/skipper
```

## Quick Start

Tell your AI agent:
- "DEWSeq + htseq-clip on eCLIP knockdown vs control with the interaction design"
- "Flipper on my Skipper output for treatment vs vehicle"
- "DESeq2 peak-level on consensus CLIPper peaks with `~ type + condition + type:condition`"
- "Aggregate adjacent significant DEWSeq windows with bedtools merge -d 100"
- "Why are my DE results dominated by housekeeping genes? Forgot the interaction term"
- "KD shows ~0 DE - check siRNA efficiency by WB and avoid TMM normalization"
- "Allele-specific differential with BEAPR"

## Example Prompts

### DEWSeq Window-Level

> "htseq-clip extract with 50nt windows / 20nt step; htseq-clip count on dedup BAM; DEWSeq design ~ type + condition + type:condition"

> "Aggregate adjacent significant windows within 100 nt"

### Flipper (Skipper-Coupled)

> "Run Flipper on the Skipper output directory for treatment vs control contrast"

### Peak-Level

> "Generate consensus peakset across all samples with bedtools merge; featureCounts on consensus; DESeq2 with interaction term"

### Design Discussions

> "Interaction term `type:condition` is the biologically meaningful test - tells me whether IP/SMInput ratio differs across conditions"

> "Don't use TMM normalization on CLIP - it forces global shifts invisible"

> "Use multiple independent siRNAs and require concordance"

### Diagnostics

> "DE results dominated by housekeeping - missing interaction term"

> "0 DE peaks in KD - check KD efficiency"

> "Significant windows scattered - need bedtools merge -d 100"

## What the Agent Will Do

1. Identify upstream peak caller (CLIPper, Skipper) and pick matching differential framework
2. For CLIPper upstream: DEWSeq + htseq-clip windows; or peak-level with consensus peakset
3. For Skipper upstream: Flipper
4. Design `~ type + condition + type:condition` to test SMInput-aware differential binding
5. Apply padj < 0.05 AND |log2FC| > 1 thresholds
6. For window-level: aggregate adjacent significant windows
7. For KD experiments: validate KD efficiency; use multiple siRNAs
8. For rescue: re-introduce siRNA-resistant RBP cDNA

## Tips

- **Interaction term is mandatory.** `~ type + condition + type:condition` tests the actual differential.
- **Consensus peakset across conditions.** Per-condition peaks bias to "treatment-only" peaks.
- **SMInput-aware design.** TMM normalization erases global shifts.
- **n >= 2 reps per condition.** n=2 marginal for DESeq2 dispersion fit.
- **DEWSeq for CLIPper; Flipper for Skipper.** Match tool to upstream caller.
- **Window-level captures narrow shifts.** Peak-level loses sub-peak resolution.
- **Aggregate adjacent windows.** `bedtools merge -d 100` for biological regions.
- **Validate KD efficiency by WB.** > 70% protein loss target.
- **Multiple siRNAs.** Off-target effects vary; require concordance.
- **Spike-in normalization for global shifts.** Drosophila S2 nuclei or similar.

## Related Skills

- clip-seq/clip-peak-calling - Upstream peaks
- clip-seq/clip-qc - QC required
- clip-seq/binding-site-annotation - Annotate DE regions
- clip-seq/clip-motif-analysis - Motifs on DE windows
- clip-seq/ago-clip-mirna-targets - miRNA target DE
- differential-expression/deseq2-basics - NB GLM
- differential-expression/de-results - DE interpretation
- chip-seq/differential-binding - DNA-protein analogue
