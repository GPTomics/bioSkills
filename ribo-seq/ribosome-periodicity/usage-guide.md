# Ribosome Periodicity - Usage Guide

## Overview

Validate Ribo-seq library quality by measuring 3-nucleotide periodicity (the signature of genuine elongating ribosomes) and calibrating read-length-specific P-site offsets. Periodicity is both the QC gate that decides whether the data is interpretable at codon level and the prerequisite calibration for ORF detection, translation efficiency, and stalling analysis.

## Prerequisites

```r
# R (primary path)
install.packages("riboWaltz")  # or BiocManager / devtools from GitHub
```

```bash
# Python alternative
pip install plastid numpy scipy pysam
```

## Quick Start

Tell your AI agent what you want to do:
- "Check 3-nucleotide periodicity in my Ribo-seq data"
- "Calibrate P-site offsets per read length with riboWaltz"
- "Decide which read lengths to keep for downstream analysis"
- "Make a metagene/metaheatmap around start codons"

## Example Prompts

### Periodicity Validation

> "Does my Ribo-seq library show 3-nt periodicity, and what is the frame-0 fraction?"

> "Is my library good enough for ORF calling, or only for gene-level counts?"

> "Why does my bacterial library look aperiodic?"

### P-site Calibration

> "Calibrate P-site offsets for each read length with riboWaltz"

> "Get offsets with plastid's metagene and psite scripts instead of R"

> "Which read lengths should I drop because they lack periodicity?"

### Visual QC

> "Make a metaprofile and metaheatmap around start and stop codons"

> "Compare frame-0 fraction across my samples"

> "Score periodicity from the CDS body, not the start peak"

## What the Agent Will Do

1. Load aligned footprints and annotation into riboWaltz (or plastid)
2. Filter read lengths by periodicity strength
3. Calibrate per-length P-site offsets (auto 5'/3' end selection)
4. Report frame-0 fraction per length as the headline metric
5. Produce metaprofile/metaheatmap plots for visual confirmation

## Tips

- **Frame-0 fraction** > 0.6-0.7 is good; ~45-60% marginal; ~33/33/33 uninterpretable at codon level
- **Offsets are per read length** and empirical - the ~12 nt value is a starting expectation, not a constant
- **A-site = P-site + 3** - use the A-site for tRNA/decoding effects
- **Trim start and stop peaks** before scoring body periodicity (they dwarf the body)
- **plastid has no metagene_analysis() function** - use the metagene/psite CLI scripts
- **Bacteria/MNase** - anchor on the 3' end and expect weaker periodicity than RNase I

## Related Skills

- riboseq-preprocessing - Produce the aligned BAM
- orf-detection - Consumes per-length P-site offsets
- translation-efficiency - Needs correct P-site positioning
- ribosome-stalling - Uses the calibrated A-site offset
