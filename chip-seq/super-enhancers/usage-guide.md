# Super-Enhancers - Usage Guide

## Overview

Identify super-enhancers (SE) from H3K27ac, MED1, or BRD4 ChIP-seq using ROSE / ROSE2 (canonical), LILY (input-subtraction variant), HOMER `-style super`, or custom hockey-stick analysis. Handles peak stitching, TSS exclusion, hockey-stick inflection, marker choice (H3K27ac vs MED1/BRD4 for BET-inhibitor predictions), spike-in normalization for cross-condition comparison, ENCODE dELS cross-referencing, and core regulatory circuitry (CRC) reconstruction. SE are operationally defined regulatory regions driving cell identity and cancer biology; see Pott & Lieb 2015 for the continuum critique.

## Prerequisites

```bash
# ROSE2 (Python 3 port; ROSE Python 2 unmaintained)
git clone https://github.com/linlabbcm/rose2.git
pip install ./rose2

# HOMER alternative
conda install -c bioconda homer

# General utilities
conda install -c bioconda samtools bedtools

# LILY (input-subtraction variant; optional)
# https://github.com/BoevaLab/LILY
```

## Quick Start

Tell the agent what to do:
- "Run ROSE2 on H3K27ac peaks with default stitching (12.5 kb) and TSS exclusion (2.5 kb)"
- "Compare SE between DMSO and BET-inhibitor treatment using spike-in normalization"
- "Use BRD4 ChIP for SE calling because I want to predict BET-inhibitor response"
- "Build core regulatory circuitry from SE-encoded TFs"
- "Cross-reference my SE against the ENCODE dELS atlas"
- "My H3K27ac is low-quality; try LILY with input subtraction"
- "Hockey-stick inflection looks wrong; verify ROSE2 output and adjust if needed"

## Example Prompts

### Standard SE calling
> "Run ROSE2 on H3K27ac peaks from K562 cells. Default 12.5 kb stitching, 2.5 kb TSS exclusion, with matched input control."

### BET-inhibitor responsiveness
> "I want to predict which SE will lose signal under JQ1 treatment. Use BRD4 ChIP for SE calling, not H3K27ac."

### Cross-condition (treatment vs control)
> "Compare SE between EZH2-inhibitor-treated and DMSO H3K27ac. Apply Drosophila spike-in normalization first because EZH2i changes H3K27ac globally (indirect; via H3K27me3 loss)."

### CRC reconstruction
> "Run the Saint-André 2016 core regulatory circuitry algorithm: identify TFs encoded by SE-associated genes whose motifs appear in their own SE and in other SE-encoded TFs."

### dELS cross-reference
> "Cross-reference my SE constituents against ENCODE dELS BED. What fraction overlap? Identify SE constituents that are NOT in dELS as potentially novel."

### Differential SE
> "Don't call SE per condition then intersect. Instead, build the union SE set, quantify signal per condition at each, and run DiffBind with appropriate normalization."

### Low-quality data
> "H3K27ac FRiP is 4% (borderline). Use LILY with input subtraction instead of ROSE2."

## What the Agent Will Do

1. **Choose marker**: H3K27ac for discovery; MED1 or BRD4 for BET-inhibitor predictions; H3K27ac+MED1 intersection for gold-standard
2. **Pre-filter peaks**: remove blacklist; remove hyper-ChIPable (top-1% input signal); remove promoter peaks if not using ROSE TSS exclusion
3. **Convert to ROSE GFF format** if needed (column order: chr, source, feature, start, end, ., strand, ., ID=N)
4. **Run ROSE2** with `-s 12500 -t 2500` defaults; provide input control for background subtraction
5. **Inspect hockey-stick plot**: verify clean inflection; if degenerate (insufficient peaks), document and consider top-N% by signal as alternative
6. **For cross-condition comparison:** spike-in normalize signal first, build union SE set, quantify at union with DiffBind
7. **For CRC:** Extract SE-encoded TFs, scan their motifs across SE, build network graph, identify connected components
8. **Cross-reference dELS** if ENCODE cCRE atlas available
9. **Assign target genes**: nearest-TSS is naive; use ENCODE-rE2G or ABC for cell-type-specific (see chip-seq/peak-annotation)
10. **Output**: SE BED, hockey-stick plot, ranked SE table, CRC network (if requested), differential SE if cross-condition
11. **Document**: marker choice, stitching/TSS-exclusion parameters, normalization (spike-in if cross-condition), assigned-gene method

## Tips

- **Use H3K27ac for discovery; BRD4 for BET-inhibitor predictions.** H3K27ac SE include inactive-but-acetylated regions; BRD4 SE define drug-responsive subset.
- **Default `-s 12500 -t 2500` is the cross-paper standard.** Deviating requires justification and breaks comparability.
- **Pre-filter hyper-ChIPable regions.** Otherwise SE list can be dominated by housekeeping / ribosomal gene clusters.
- **Cross-condition comparison without spike-in is unreliable.** Global signal shifts (HDACi, BETi, EZH2i) confound threshold-based SE counting.
- **Build a union SE set for differential analysis, then quantify with DiffBind.** Calling SE per condition and intersecting BEDs is the wrong approach.
- **Most SE constituents are not functional.** Hnisz 2017 / Dukler 2017 CRISPR-tiling shows 1-3 essential constituents per SE.
- **CRC algorithm requires SE + TF motif annotations + gene-SE mapping.** Output is the network graph; biology comes from interpreting connected components.
- **ENCODE dELS atlas cross-reference is a sanity check.** Most SE constituents should overlap dELS for high-quality data.
- **Hockey-stick inflection is sensitive to peak count.** ≥ 5000 enhancer peaks recommended.

## Troubleshooting

### Too few SE called

1. Hockey-stick inflection failed -> check `_Plot_points.png`; need ≥ 5000 enhancer peaks
2. Stitching distance too small for biology -> try `-s 15000` or `-s 25000` for tissue samples
3. TSS exclusion too aggressive -> reduce `-t 500` if including promoter-proximal SE
4. Input control too strong -> consider running without input or with LILY's variant

### Too many SE called

1. Hyper-ChIPable artifacts inflated input -> filter blacklist + custom top-1% input signal
2. Insufficient differentiation in signal distribution -> hockey-stick fails; check `_Plot_points.png`
3. Threshold not biologically meaningful -> use top 5-10% by signal as alternative

### SE counts differ between conditions

This is expected when:
- Real biology (gain/loss of regulatory programs)
- Global signal shift (HDACi, BETi, EZH2i) confounding threshold
- Different sequencing depth between conditions

Fix the confounding before interpreting: spike-in normalize and quantify on union SE set.

### CRC network has too many components / too few

1. Too many TFs (default lists may include non-master TFs) -> restrict to TFs with motifs in own SE
2. Too few cross-bindings -> check motif quality (use HOCOMOCO v12 not JASPAR for robust calls)
3. Cell type with no master regulators (rare) -> CRC may not apply

### Differential SE all "novel" gain/loss

Without spike-in normalization, this is often artifact. With spike-in, expect ~30-50% conserved SE between conditions and a smaller "differential" set.

## Related Skills

- chip-seq/peak-calling - Generate H3K27ac / MED1 / BRD4 peaks
- chip-seq/chipseq-qc - Filter hyper-ChIPable peaks before SE calling
- chip-seq/spike-in-normalization - Mandatory for cross-condition SE comparison
- chip-seq/differential-binding - Quantitative differential testing on union SE set
- chip-seq/peak-annotation - Annotate SE-associated genes; cross-reference dELS
- chip-seq/cut-and-run-tag - SE on CUT&RUN/CUT&Tag H3K27ac (E. coli spike-in)
- atac-seq/enhancer-gene-linking - ENCODE-rE2G for SE-target gene assignment
- data-visualization/genome-tracks - SE region visualization
