# DIA Analysis - Usage Guide

## Overview
Identifies and quantifies proteins from data-independent acquisition (DIA) mass spectrometry by scoring reconstructed fragment-chromatogram peak groups against a decoy null, then filtering the output at the correct q-value level and context. Covers DIA-NN library-free (directDIA), predicted-library, and library-based routes, with notes on Spectronaut, OpenSWATH, and EncyclopeDIA. The crux: every wide-window MS2 is chimeric, so the engine deconvolves rather than matches, and "1% FDR" only means something once the level (precursor/peptide/protein-group) and context (run vs experiment-wide/global) are named.

## Prerequisites
```bash
pip install pandas pyarrow numpy
# CLI: DIA-NN (recommended), MSFragger-DIA/FragPipe, OpenSWATH, EncyclopeDIA
# Commercial: Spectronaut
# Staggered/overlapping data: ProteoWizard msconvert (demultiplexing) before search
```

## Quick Start
Tell your AI agent what you want to do:
- "Run DIA-NN on my mzML files via the predicted-library route at 1% FDR"
- "Search my DIA data against a spectral library I built"
- "Load the DIA-NN report.parquet, filter q-values correctly, and build a protein matrix"
- "My cohort hits do not reproduce -- which q-value should I be filtering on?"
- "Demultiplex my staggered-window data before searching"

## Example Prompts

### Library-Free and Predicted-Library Analysis
> "Run DIA-NN via the predicted-library route against UniProt human, two-pass, with auto mass accuracy"

> "Set up library-free directDIA for my Astral 2-Th narrow-window runs"

> "Generate an in-silico predicted library and search my DIA data against it"

### Library-Based Analysis
> "Search my DIA data against the chromatogram library I built in EncyclopeDIA"

> "Run DIA-NN with my Prosit-predicted library for targeted extraction"

> "Use a previously generated .speclib so I do not redigest the FASTA"

### FDR and Output Filtering
> "Filter my DIA-NN report at precursor and protein-group level, and global PG for the matrix"

> "Why is my pg_matrix protein count lower than the report count?"

> "My per-run 1% FDR cohort has irreproducible hits -- fix the FDR context"

### Results Processing
> "Load report.parquet, drop unquantified zeros, and log2-transform a MaxLFQ protein matrix"

> "Pivot the filtered long report into a wide protein-by-run matrix"

> "Set up auditable run/experiment/global q-values with OpenSWATH and PyProphet"

## What the Agent Will Do
1. Identify the acquisition design (fixed/variable/staggered windows, diaPASEF, narrow-window Astral) and demultiplex staggered data at conversion if needed
2. Choose the route -- predicted-library (default), library-free directDIA, or library-based -- per the decision tree
3. Run DIA-NN (or the chosen engine) with auto mass accuracy and two-pass global FDR
4. Read report.parquet and filter at the correct LEVEL (precursor + protein-group) and CONTEXT (run + global for matrices)
5. Convert unquantified zeros to NA, then hand the log2 matrix to quantification and differential-abundance

## Tips
- The predicted-library route is the modern default; raw directDIA on wide-window data needs careful two-pass FDR.
- "1% FDR" alone is ambiguous; always state the level (precursor/peptide/protein-group) and context (run/experiment/global).
- For cohorts, filter on Global.PG.Q.Value, not the per-run Q.Value, or the experiment-wide error inflates.
- DIA-NN 1.9+ writes report.parquet by default; loaders assuming report.tsv silently break.
- A matrix protein count below the report count is expected (extra 5% run-specific PG filter), not data loss.
- Let DIA-NN auto-optimize tolerances with --mass-acc 0 rather than hard-coding ppm from another instrument.
- Aim for >= 6 MS2 points across each LC peak; if quant is noisy, the window/cycle design may be the cause.
- Convert DIA-NN's 0 (not-quantified) to NA before log2 or normalization.

## Related Skills

- spectral-libraries - Build experimental, chromatogram, or predicted libraries to search against
- quantification - Normalization, MaxLFQ roll-up, and matrix summarization after filtering
- differential-abundance - Moderated statistical testing of the protein matrix
- proteomics-qc - Per-run ID counts, missing-value rates, and acquisition QC
- data-import - Convert and load raw vendor/mzML data before the search
- peptide-identification - DDA spectrum-to-peptide matching (the non-DIA counterpart)
- workflows/proteomics-pipeline - End-to-end DIA-to-differential-abundance orchestration
