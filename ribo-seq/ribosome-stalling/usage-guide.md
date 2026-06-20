# Ribosome Stalling - Usage Guide

## Overview

Detect ribosome pausing and stalling at codon resolution to study elongation dynamics, codon dwell times, pause motifs, and ribosome collisions. The decisive judgment is whether a pause is real biology or a cycloheximide artifact, followed by using local-relative occupancy metrics at the A-site rather than a global z-score.

## Prerequisites

```bash
pip install plastid numpy scipy biopython twobitreader
```

## Quick Start

Tell your AI agent what you want to do:
- "Find ribosome pause sites with a local pause score"
- "Calculate A-site codon occupancy"
- "Check whether my pauses are real or a cycloheximide artifact"
- "Look for polyproline and poly-basic stalling"

## Example Prompts

### Pause Detection

> "Score pauses as occupancy over the gene mean, not a global z-score"

> "Was my library cycloheximide-treated? Can I trust the dwell times?"

> "Which pauses also show up as disome peaks?"

### Codon Analysis

> "Calculate A-site codon occupancy across the transcriptome"

> "How strong is the correlation between codon occupancy and tRNA abundance?"

> "Normalize each gene to its own mean before pooling codons"

### Motif Analysis

> "What amino-acid motifs are enriched at my pause sites?"

> "Find polyproline-associated stalling"

> "Extract the A-site-centered context around pauses"

## What the Agent Will Do

1. Confirm the harvest protocol (drug, freezing) before any dwell claim
2. Map footprints to the calibrated A-site offset
3. Compute per-codon occupancy and local-relative pause scores
4. Require an adequate per-gene coverage floor
5. Extract motif context and cross-check disome/collision evidence

## Tips

- **Cycloheximide flips conclusions** - dwell times are only valid on flash-frozen no-drug data
- **A-site for decoding** - tRNA/codon effects register at the A-site (P-site + 3)
- **Local-relative metrics** - pause score = occupancy / gene mean, not a global z-score
- **Normalize per gene first** - then pool codons (mean-of-ratios)
- **Coverage floor** - need a few hundred in-frame footprints per gene, not 100
- **Disomes confirm pauses** - a collision peak coinciding with a monosome pause is strong evidence

## Related Skills

- ribosome-periodicity - Calibrate the A-site offset
- orf-detection - Locate ORFs containing pause sites
- initiation-site-mapping - Separate initiation drugs from elongation pausing
- translation-efficiency - Gene-level translation context
