# Differential Binding Analysis - Usage Guide

## Overview

Identifies differentially bound ChIP-seq regions between conditions. Distinguishes three distinct normalization problems (composition bias, trended/abundance-dependent bias, global shift) and matches each to the appropriate fix (RLE/TMM on reads-in-peaks, TMM on background bins, or spike-in scaling). Covers DiffBind (BAM + peaks), DESeq2/edgeR/PyDESeq2 (count matrix), csaw (sliding windows for global-shift-robust testing), NormR (control-aware), MAnorm2, and ChIPseqSpikeInFree.

## Prerequisites

```r
BiocManager::install(c('DiffBind', 'DESeq2', 'edgeR', 'csaw', 'apeglm', 'rtracklayer'))
BiocManager::install(c('normr', 'ChIPseqSpikeInFree'))   # optional
```

```bash
pip install pydeseq2 pandas numpy
```

## Quick Start

Tell the agent what to do:
- "Run DiffBind on my sample sheet with default normalization for a typical TF ChIP"
- "I'm comparing H3K27ac in HDAC-inhibitor-treated vs control cells; choose the right normalization"
- "Use csaw with sliding windows to detect global shifts in H3K27me3 after EZH2 inhibition"
- "I have spike-in Drosophila reads; apply ChIP-Rx scaling to my DESeq2 analysis"
- "Run differential binding with both reads-in-peaks and background-bin normalization and reconcile"
- "Suspect global shift but no spike-in; run ChIPseqSpikeInFree for post-hoc detection"
- "Run NormR for control-aware enrichment/depletion calling"

## Example Prompts

### Standard TF ChIP (local changes)
> "Run DiffBind on 6 samples (3 ctrl, 3 treat) for FOXA1 ChIP-seq. Use the default reads-in-peaks RLE normalization since the biology is localized TF rebinding."

### Global shift (HDACi / BETi / EZH2i)
> "Compare H3K27me3 in EZH2-inhibitor-treated vs DMSO. EZH2 inhibition reduces H3K27me3 genome-wide; use spike-in scaling. I have Drosophila S2 read counts per sample."

### Broad mark, no spike-in
> "Run csaw with 1 kb windows and background-bin TMM normalization for H3K27me3 differential. Merge significant windows within 5 kb."

### Suspected composition bias
> "DiffBind with default settings shows almost no differential peaks but the IGV signal is obviously different. Switch to `background=TRUE` and re-analyze."

### Post-hoc global shift detection
> "I didn't run spike-in. Use ChIPseqSpikeInFree to detect whether there's a global shift in my BET-inhibitor experiment."

### Sample sheet preparation
> "Create a DiffBind sample sheet from my 8 BAMs + narrowPeak files. Columns: SampleID, Condition, Replicate, bamReads, bamControl, Peaks, PeakCaller."

### Method reconciliation
> "DiffBind reports 50 differential peaks; csaw windows reports 300. Walk through reconciliation logic and decide which to trust."

## What the Agent Will Do

1. **Diagnose the normalization problem**: ask about expected global vs local changes (HDACi/BETi/EZH2i/dosage → global; standard TF perturbation → local); check MA plot loess curve after a default run
2. **Choose the appropriate tool/normalization** per the three-problem framework
3. **Build sample sheet** with SampleID, Condition, Replicate, bamReads, bamControl, Peaks, PeakCaller
4. **Run differential binding:**
   - DiffBind: `dba.count()` → `dba.normalize()` → `dba.analyze()`
   - DESeq2/edgeR: count matrix → DESeqDataSet → test
   - csaw: `windowCounts()` → `normFactors()` → `glmQLFTest()` → merge
   - NormR: `diffR()` joint binomial mixture on bins
5. **Spike-in integration** if applicable; compute scaling factors from Drosophila/E. coli reads, pass via `sizeFactors()` or `DiffBind` spike-in fields
6. **Verify normalization**: `dba.normalize(obj, bRetrieve = TRUE)` for DiffBind; inspect MA loess curve
7. **Internal-control sanity check**: for spike-in, confirm blacklist regions show no signal change post-scaling
8. **Output**: differential peaks BED + log2FC + padj; volcano + MA + PCA + heatmap plots
9. **Document**: version, normalization choice, design formula, padj threshold, spike-in metadata in methods

## Tips

- **Always check the MA plot loess after default run.** A shifted curve indicates wrong normalization. This is a 60-second sanity check that saves hours of misinterpretation.
- **DiffBind 3.x defaults changed.** `summits=200` recenters narrow peaks (set FALSE for broad histones); blacklist filter on by default; `dba.normalize()` now required.
- **Pre-filtering for ChIP-seq is lighter than RNA-seq.** Peaks are already enriched regions; aggressive filtering removes condition-specific gains.
- **LFC shrinkage is for ranking, not significance.** apeglm-shrunk fold changes are better for volcano plots and downstream ranking but do NOT change padj.
- **edgeR QL-F controls type-I error better than DESeq2 Wald** for small replicate counts (n=2-3).
- **Spike-in is the gold standard for global shifts.** Apply via `sizeFactors()` (read-level), never multiply peak counts.
- **For broad marks, use full peak width (summits=FALSE) or csaw with 1-2 kb windows.** Narrow recentering loses domain-level biology.
- **PCA before differential**: samples should cluster by condition, not batch. If batch dominates, add batch term to design.
- **Document spike-in scaling factors in methods.** Hammond Norris 2024: ~25% of published spike-in papers have flawed normalization detectable from methods.

## Troubleshooting

### Almost no differential peaks despite obvious effect

Usually wrong normalization for the biology:
1. Check MA loess; if shifted off y=0, switch to `background=TRUE` (DiffBind) or spike-in
2. Global-shift drug? Spike-in required; no algorithmic fix works
3. Pre-filtering too aggressive → `rowSums >= 1` only

### Differential peaks have wrong sign

1. Composition bias inverted the signal; switch normalization
2. Contrast direction reversed; verify `contrast = c('condition', 'treat', 'ctrl')` for positive log2FC = up in treat

### High noise / many marginal peaks

1. Loose ENCODE-pattern input peaks; pre-filter to IDR or naive-overlap passing
2. Low FRiP / failing ChIP; see chipseq-qc; failing QC samples produce noisy differential

### DiffBind very slow

`dba.count(obj, bParallel = TRUE)` for parallel counting. For large peak sets (>100k), consider csaw with windows instead.

### Spike-in scaling gives nonsense

1. Applied to peak counts instead of `sizeFactors()` → fix application layer
2. Spike-in reads not deduplicated → re-extract from BAM
3. Spike-in saturated (too many reads) → verify titration linearity
4. Spike-in mismatch with target → confirm spike-in genome (Drosophila for human/mouse; E. coli for CUT&RUN/Tag)
5. Internal control changes → see chipseq-qc and re-validate

### `Error: condition has only one level`

Either factor has only one level or `Condition` column missing from sample sheet. Add at least two distinct conditions with ≥2 replicates each.

## Related Skills

- chip-seq/peak-calling - Generate consensus peaks for DiffBind input
- chip-seq/chipseq-qc - Replicate concordance required before differential testing
- chip-seq/spike-in-normalization - Spike-in workflow for global-shift experiments
- chip-seq/cut-and-run-tag - E. coli spike-in (CUT&RUN) vs Drosophila (ChIP-Rx)
- chip-seq/peak-annotation - Annotate differential peaks to genes/cCREs
- differential-expression/deseq2-basics - DESeq2 patterns
- differential-expression/edger-basics - edgeR QL-F framework
- atac-seq/differential-accessibility - Parallel ATAC workflow
