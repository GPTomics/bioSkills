# Quantification - Usage Guide

## Overview
Reconstruct a protein-by-sample abundance matrix from mass spectrometry signals, choosing a summarizer and normalizer that match where the signal physically came from. Label-free (LFQ/MaxLFQ), isobaric (TMT/iTRAQ reporter ions), and metabolic (SILAC) approaches each carry an irreducible error set by their measurement physics, and the peptide-to-protein summarization choice changes the answer more than the downstream statistical test does.

## Prerequisites
```bash
pip install numpy pandas scipy
# R packages:
# BiocManager::install(c("MSstats", "MSstatsTMT", "MSnbase"))
# install.packages("iq")
```

## Quick Start
Tell your AI agent what you want to do:
- "Summarize my MaxQuant peptides to protein level with MSstats"
- "Run the real MaxLFQ algorithm on my peptide intensity matrix"
- "Extract TMT reporter ions and correct for isotope impurity"
- "Bridge my multiple TMT plexes with an IRS reference channel"
- "Compute SILAC heavy/light ratios and check for Arg-to-Pro conversion"

## Example Prompts

### Label-Free Summarization
> "Convert my MaxQuant evidence.txt into normalized protein-level abundances using MSstats Tukey median polish"

> "Run iq::maxLFQ on my peptide quant matrix instead of median centering"

> "Compare TMP and MaxLFQ summarization on the same data and report where the answer moves"

### Normalization
> "Median-center my label-free intensity matrix to correct sample loading"

> "Apply sample-loading normalization then IRS to bridge my three TMT plexes"

### TMT/iTRAQ Processing
> "Extract TMT10 reporter ions from my mzML and apply lot-specific impurity correction"

> "Explain why MS2 reporter quant compresses my fold changes and whether SPS-MS3 helps"

### SILAC
> "Compute SILAC log2 ratios but keep proteins present only in the heavy channel"

> "Check my SILAC labeling efficiency and flag Arg-to-Pro conversion before trusting ratios"

## What the Agent Will Do
1. Identify the labeling strategy (LFQ, TMT/iTRAQ, SILAC) and the signal origin it implies
2. Load peptide/PSM/reporter data and replace MaxQuant zeros with NaN before any transform
3. Apply normalization appropriate to the data (median centering for LFQ; sample-loading then IRS across TMT plexes)
4. Summarize peptides to proteins (MSstats Tukey median polish or iq::maxLFQ), reporting the summarizer chosen
5. For TMT, extract and impurity-correct reporter ions with lot-specific Certificate values
6. For SILAC, verify labeling efficiency and Arg->Pro conversion and preserve on/off biology
7. Run a sensitivity check across at least two summarizers and hand the matrix to differential-abundance

## Tips
- Replace MaxQuant zeros with NaN first; log2(0) is -inf and corrupts every transform
- Median centering is a normalizer, not a summarizer, and it is NOT MaxLFQ; call iq::maxLFQ for the real algorithm
- Report the summarization method as prominently as the statistical test; it moves the answer more
- TMT ratio compression is physical and cannot be normalized away; attack it with SPS-MS3, narrow windows, or ion mobility
- TMT plexes are not comparable without a per-plex reference channel and an IRS bridge
- SILAC: a protein present only in the heavy channel is real biology, not a value to discard as NaN
- Spectral counting / NSAF is largely obsolete; keep it as historical context, not a recommendation

## Related Skills
- data-import - Parse MaxQuant/DIA-NN outputs and pick the right intensity column before quantifying
- protein-inference - Razor/shared-peptide group assignment that determines which protein a peptide counts toward
- differential-abundance - Statistical testing, missing-value modeling, and the downshift false-positive trap
- proteomics-qc - CV, correlation, and PCA checks that confirm normalization worked
- dia-analysis - DIA fragment-level MaxLFQ and DIA-NN execution
- differential-expression/de-results - Shared empirical-Bayes and FDR conventions for expression matrices
- data-visualization/heatmaps-clustering - Visualize the normalized abundance matrix
- workflows/proteomics-pipeline - End-to-end pipeline that calls this skill for the quant step
