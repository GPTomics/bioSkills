# cfDNA Preprocessing - Usage Guide

## Overview
Preprocess plasma cell-free DNA sequencing data so the diagnostic signal survives: choose a consensus strategy (single-strand UMI vs duplex), run the fgbio align->group->consensus->re-align chain with the correct flags, and avoid the cfDNA dedup trap where naive coordinate dedup collapses nucleosome-coincident independent molecules. The fragment-length distribution is treated as structured biological signal and a pre-analytical QC instrument, not noise.

## Prerequisites
```bash
conda install -c bioconda fgbio bwa samtools
pip install pysam numpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Build UMI consensus reads from my targeted cfDNA panel BAM"
- "Set up a duplex consensus pipeline for sub-0.1% VAF detection"
- "Do minimal preprocessing for sWGS tumor-fraction estimation"
- "Check the insert-size distribution of my cfDNA library for gDNA contamination"
- "I have no UMIs - how should I dedup this cfDNA quantitatively?"

## Example Prompts

### Single-Strand UMI Consensus
> "Process my targeted panel cfDNA with single-strand UMIs: extract UMIs, align, group by adjacency, call molecular consensus, and filter."

> "My read structure has a UMI stem - help me write the correct fgbio read structure so the stem does not bleed into the template."

### Duplex Consensus
> "Set up a duplex consensus pipeline so I can detect variants below 0.1% VAF, using paired grouping and a true-duplex filter."

> "Explain why I should call duplex consensus permissively and filter afterward instead of setting min-reads high on the caller."

### Minimal Processing
> "I am running sWGS for tumor fraction - what minimal preprocessing do I need, and should I skip consensus calling?"

### QC and Interpretation
> "Plot my cfDNA insert-size distribution and tell me whether the ~167 bp peak and sawtooth look healthy or contaminated."

> "My fragment mode is about 10 bp short of 167 - is that gDNA contamination or a library-prep artifact?"

## What the Agent Will Do
1. Establish the library prep (dsDNA vs ssDNA/adaptase) and UMI design to set fragment-size and consensus expectations.
2. Extract UMIs into the `RX` tag with the correct read structure (including the skip/stem token).
3. Align with `bwa mem -Y` so tag-bearing short-fragment sequence is preserved.
4. Group reads by UMI with `adjacency` (simplex) or `paired` (duplex, mandatory).
5. Call single-strand or duplex consensus permissively, then RE-align the unmapped consensus reads.
6. Apply the real quality gate with `FilterConsensusReads` using the appropriate `--min-reads` strand specification.
7. Read the insert-size histogram for the mode, sawtooth, and long-fragment fraction as a QC check.

## Tips
- **Prep gates everything** - record dsDNA vs ssDNA/adaptase; a fragmentomics pipeline tuned for one will misread the other.
- **Duplex needs paired grouping** - `--strategy adjacency` cannot reconstruct strand pairing; use `--strategy paired` for any duplex workflow.
- **Call permissively, filter strictly** - `CallDuplexConsensusReads --min-reads` is a pre-filter; do the real filtering in `FilterConsensusReads`.
- **Re-align consensus reads** - consensus output is unmapped by design; skipping the second alignment yields garbage coordinates.
- **Never naive-dedup no-UMI cfDNA** - nucleosome-positioned ends make coordinate dedup delete real independent molecules and deflate VAF.
- **Mind the singleton tax** - at picogram input, requiring two reads per family discards genuine low-VAF evidence; favor recovery for detection.
- **Do not size-select then measure fragmentomics** - selection conditions on length and biases every length-derived feature.

## Related Skills
- analytical-validation - the LoD/molecule-counting framework input quality feeds
- fragment-analysis - fragmentomics consumes the preprocessed fragment ends
- ctdna-mutation-detection - consensus reads feed low-VAF calling
- tumor-fraction-estimation - sWGS minimal-processing path
- alignment-files/duplicate-handling - general dedup vs the cfDNA UMI caveat
- read-qc/quality-reports - upstream read QC
