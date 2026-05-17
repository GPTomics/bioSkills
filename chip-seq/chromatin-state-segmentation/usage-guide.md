# Chromatin State Segmentation - Usage Guide

## Overview

Segment the genome into chromatin states by integrating multiple histone modification and chromatin factor ChIP-seq tracks. Covers ChromHMM (canonical multivariate HMM on binarized signal), Segway (Dynamic Bayesian Network on continuous signal), EpiSegMix (flexible-distribution HMM, 2024), EpiLogos (multi-biosample visualization), IDEAS (cell-type-aware joint), and full-stack ChromHMM (Vu Ernst 2022). Embeds state-count selection logic (15 vs 18 vs 25), Roadmap canonical model interpretation, bin-size trade-offs, and downstream OverlapEnrichment / NeighborhoodEnrichment for functional characterization.

## Prerequisites

```bash
# ChromHMM (canonical)
wget http://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM.zip
# Requires Java 8+

# Segway (continuous signal alternative)
pip install segway

# Other utilities
conda install -c bioconda samtools bedtools
```

## Quick Start

Tell the agent what to do:
- "Run ChromHMM on H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3 BAMs at 15-state resolution for GM12878"
- "Binarize ChIP-seq BAMs with `BinarizeBam`, learn 18-state model, then OverlapEnrichment with CGI / repeats / lncRNA / TSS anchors"
- "Compare ChromHMM 15-state vs 18-state vs 25-state segmentations and pick the smallest N where biology is interpretable"
- "Apply the precomputed Roadmap 25-state model to my new sample for cross-cell-type consistency"
- "Use EpiSegMix instead of ChromHMM for continuous signal handling"
- "Visualize ChromHMM segmentations across 12 cell types with EpiLogos"

## Example Prompts

### Standard 15-state ChromHMM
> "Build cellMarkFileTable for GM12878 with 5 core marks (H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3) and matching input controls. Binarize with `BinarizeBam`, learn a 15-state model, generate emission heatmap and segmented BED."

### Multi-cell-type joint segmentation
> "I have 5 marks in 3 cell types (GM12878, K562, HepG2). Build a joint cellMarkFileTable and learn a 25-state model so states are comparable across cell types. Map state labels to Roadmap conventions where possible."

### Functional characterization
> "After learning the model, run `OverlapEnrichment` against CGI / refseq exons / lncRNA / repeats / ENCODE blacklist anchor files. Run `NeighborhoodEnrichment` against TSS positions."

### Bivalent states
> "Use the 18-state model (extends 15-state with bivalent enhancer + bivalent flanking states) for an ESC-like cell type where bivalent domains are common (Bernstein 2006)."

### Continuous-signal alternative
> "ChromHMM binarized; switch to Segway for continuous-signal segmentation that preserves signal magnitudes."

### Cross-cell-type visualization
> "Visualize learned 25-state segmentations across 12 biosamples in EpiLogos to identify tissue-specific regulatory states."

### Apply precomputed model
> "Apply the Roadmap Epigenomics 25-state model to my new tissue sample without retraining (use `MakeSegmentation` with the precomputed model)."

## What the Agent Will Do

1. **Verify mark panel**: need ≥ 4-5 marks for meaningful states; 5 core (H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3) for canonical Roadmap-compatible models
2. **Build `cellMarkFileTable.txt`**: tab-separated: cell_type, mark, BAM, control_BAM
3. **Binarize BAMs**: `ChromHMM BinarizeBam` produces per-chromosome 200 bp bin matrices
4. **Train model** at multiple N; compare 15, 18, 25 states; choose by emission-matrix interpretability
5. **Inspect emission heatmap**: verify states map to canonical biology (TssA, EnhA, EnhWk, Tx, ReprPC, Het, Quies)
6. **Run downstream enrichment:**
   - `OverlapEnrichment` against CGI, exons, lncRNA, repeats, blacklist
   - `NeighborhoodEnrichment` against TSS, transcription factor binding sites
7. **Cross-cell-type consistency**: use joint model OR apply same precomputed model across samples
8. **Visualize**: EpiLogos for multi-sample; per-state BED for browser
9. **Document**: mark panel, bin size, state count chosen, emission heatmap, cell type(s)

## Tips

- **Need ≥ 5 marks for a canonical 15-state model.** Fewer marks = simpler peak-based annotation (peak-annotation) is more appropriate.
- **5 core marks (Roadmap):** H3K4me3 (active TSS), H3K27ac (active regulatory), H3K4me1 (enhancers), H3K36me3 (transcribed), H3K27me3 (Polycomb). Add H3K9me3 for heterochromatin.
- **Train at multiple N (15, 18, 25).** Compare emission matrices; choose smallest N where states are biologically distinct.
- **Default 200 bp bin works for most uses.** Reduce to 100 or 50 only for sharp boundary studies.
- **Use sonicated input control for binarization.** IgG is not appropriate for histone mark ChromHMM.
- **For cross-cell-type comparison, use joint segmentation or precomputed model.** Independently-trained per-cell-type models produce non-comparable state labels.
- **Segway uses continuous signal (bigWig); ChromHMM uses binarized signal.** Segway preserves magnitudes; ChromHMM simpler and more standardized.
- **EpiLogos is for visualization, not segmentation.** Use after ChromHMM/Segway segmentation across many biosamples.
- **Full-stack ChromHMM (Vu Ernst 2022)** has 100 states across 1032 datasets / 127 reference epigenomes; useful for applying a comprehensive cross-tissue annotation.
- **ChIP-seq normalization affects ChromHMM.** Spike-in-scaled BAMs are appropriate for cross-condition state comparison; standard BAMs for within-cell-type analysis.

## Troubleshooting

### State count too few; biologically distinct regions lumped

Train with more states (try 18 or 25). Inspect emission heatmap to verify distinct states.

### State count too many; emission matrix has redundant states

Some states have nearly identical emission profiles. Use smaller N; or merge states post-hoc.

### Mark panel insufficient

If only 2-3 marks available, segmentation is unreliable. Either:
1. Add more marks
2. Use simpler peak-based annotation (chip-seq/peak-annotation)
3. Use a precomputed model with matching mark panel

### `OutOfMemoryError` in Java

Increase JVM heap: `java -mx32G -jar ChromHMM.jar ...`. For very large samples, use `-mx64G`.

### Binarization shows no signal

1. Wrong control file (IgG instead of input)
2. Mark BAM has no actual signal (failed ChIP)
3. Threshold too strict; verify via `samtools idxstats` that marks have coverage

### State labels don't match Roadmap

Trained model independently. Apply Roadmap precomputed model OR map states by emission similarity to Roadmap states.

### Cross-replicate state assignments differ

Insufficient marks or biological variability. Increase mark panel or use joint training across replicates.

### `MakeSegmentation` error with precomputed model

Mark panel must match the model exactly (same marks in same order as binarized data).

## Related Skills

- chip-seq/peak-calling - Per-mark peak calling
- chip-seq/chipseq-qc - Per-mark QC before integration
- chip-seq/peak-annotation - ENCODE cCRE classification (PLS / pELS / dELS / CTCF-only / DNase-H3K4me3)
- chip-seq/spike-in-normalization - Spike-in for cross-condition state comparison
- atac-seq/single-cell-atac - Multi-modal integration with chromatin states
- machine-learning/model-validation - State-count selection
- data-visualization/genome-tracks - Visualize state segmentations
- gene-regulatory-networks/coexpression-networks - Cross-reference chromatin states with expression
