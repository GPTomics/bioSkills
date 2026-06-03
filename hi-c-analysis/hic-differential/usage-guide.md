# Hi-C Differential Analysis - Usage Guide

## Overview

This skill compares Hi-C contact maps between conditions at the correct scale: differential bin-pair contacts, differential A/B compartments, differential TAD boundaries, and differential loops. The central trap it guards against is that matrix balancing (ICE/KR/SCALE) makes each map internally consistent but does NOT make two maps cross-comparable -- depth, cis/trans ratio, and the shared distance-decay P(s) dominate any naive log2 ratio. The right approach is a between-sample, distance-stratified normalization (the M-D plot, multiHiCcompare's cyclic loess) followed by a replicate-aware test matched to the scale of the question. Because cooltools deliberately provides only feature-extraction parts (no turnkey test), the differential statistics live in R/Bioconductor (multiHiCcompare, diffHic, dcHiC, hicrep), with HiCRep SCC as the mandatory replicate-reproducibility gate and CNV correction required for any cancer comparison.

## Prerequisites

```bash
pip install cooler cooltools bioframe
# Differential statistics live in R/Bioconductor:
# install.packages('BiocManager')
# BiocManager::install(c('multiHiCcompare', 'diffHic', 'edgeR', 'hicrep'))
# dcHiC: clone github.com/ay-lab/dcHiC and follow its conda env setup
```

Coolers must be BALANCED before feature extraction, and a `.mcool` must be addressed by a single-resolution URI (`file.mcool::/resolutions/10000`). A credible differential claim needs n>=2 (ideally 3) biological replicates per condition.

## Quick Start

Tell your AI agent what you want to do:
- "Compare Hi-C between treatment and control with replicates"
- "Gate my replicates with HiCRep SCC before calling differences"
- "Find differential bin-pair contacts with multiHiCcompare"
- "Find differential A/B compartments with dcHiC"
- "Find differential TAD boundaries from delta insulation"
- "CNV-correct my tumor-vs-normal comparison before testing"

## Example Prompts

### Replicate QC gating
> "I have two replicates per condition. Before any differential call, compute HiCRep SCC for the within-condition replicate pairs and the between-condition pairs at 50kb up to 5Mb, and tell me whether within exceeds between."

### Differential bin-pair contacts
> "I have balanced coolers for two conditions, two replicates each. Build a multiHiCcompare experiment from sparse upper-triangular tables, run cyclic loess M-D normalization, run the edgeR exact test, and give me the differential bin-pairs at p.adj < 0.1 with the M-D diagnostic plot."

### Differential compartments
> "Run dcHiC on my four samples (two conditions, two reps) to find differential A/B compartments at 100kb on hg38, including graded shifts that don't flip A/B, and point me at the fdr_result directory."

### Differential boundaries
> "Compute cooltools insulation at a 200kb window for both conditions, difference the scores, and rank the boundaries that change the most."

### Cancer comparison
> "This is tumor vs matched normal. CNV-correct with diffHic normalizeCNV using the marginal log-ratios before the edgeR test so I don't call copy-number blocks as differential contacts."

## What the Agent Will Do

1. Gate reproducibility with HiCRep SCC (within- vs between-condition) before any call.
2. Choose the method by scale: bin-pair (multiHiCcompare/diffHic), compartment (dcHiC), boundary (delta insulation), or loop (diffloop/DiffHiChIP).
3. Apply a distance-stratified between-sample normalization (cyclic loess / diffHic loess offsets), never a naive log2 ratio.
4. Filter low-abundance bin-pairs on a contrast-independent statistic before testing.
5. Run a replicate-aware NB-GLM test and control FDR (distance-aware where possible).
6. For cancer/aneuploid data, CNV-correct (diffHic normalizeCNV / OneD) before testing.
7. Inspect M-D diagnostics and report the comparison scale explicitly.

## Tips

- Balance is within-map only -- normalize BETWEEN samples, distance-stratified, before differencing.
- Report HiCRep SCC, never raw Pearson; the P(s) decay makes Pearson look reproducible for unrelated maps.
- No replicates means no honest FDR; treat Selfish/FIND n=1 output as descriptive only.
- Match the tool to the scale: compartment (dcHiC), TAD boundary (delta insulation), loop (diffloop), bin-pair (multiHiCcompare/diffHic) are different objects.
- A 1Mb fold-change is not the same currency as a 50kb fold-change; pooling distances into one FDR over-calls short-range and under-calls long-range.
- CNV-correct before ANY tumor-vs-normal comparison; balancing makes the CNV artifact worse.
- Verify `get.scc` and cooltools signatures against the installed version before chaining (both have version skew).

## Related Skills

- compartment-analysis - Per-condition A/B eigenvectors that dcHiC differences
- tad-detection - Per-condition insulation scores for delta-insulation boundary tests
- loop-calling - Per-condition loop calls fed to diffloop/DiffHiChIP
- matrix-operations - Balancing and expected/O-E that precede any comparison
- hichip-plac-loops - Peak-anchored loop calls and DiffHiChIP for HiChIP/PLAC-seq
- hic-data-io - Load and convert the cooler files this skill compares
- hic-visualization - Render differential maps and split-view comparisons
- chip-seq/peak-annotation - Annotate differential anchors/boundaries with TF peaks
- genome-intervals/overlap-significance - Permutation test for differential-feature enrichment
- differential-expression/de-results - The edgeR/FDR mental model reused here
