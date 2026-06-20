# Translation Initiation Site Mapping - Usage Guide

## Overview

Map translation initiation sites (TIS) at single-nucleotide resolution from initiation-drug ribosome profiling (TI-seq), including non-AUG and alternative starts. This is a distinct analysis from elongation ORF detection: it asks which start codon is used, and it requires a dedicated harringtonine, lactimidomycin (GTI-seq/QTI-seq), or retapamulin (Ribo-RET) library.

## Prerequisites

```bash
pip install ribotish
# PRICE is distributed as GEDI (Java)
```

## Quick Start

Tell your AI agent what you want to do:
- "Map translation initiation sites from my harringtonine Ribo-seq"
- "Find non-AUG and uORF start codons with Ribo-TISH"
- "Detect cryptic initiation with PRICE"
- "Compare initiation between two conditions"

## Example Prompts

### Initiation Detection

> "Run Ribo-TISH predict on my LTM and elongation libraries"

> "Enable near-cognate (CUG/GUG) start codons in the search"

> "Which uORFs initiate at non-AUG starts?"

### Experiment-Specific

> "My data is harringtonine - set the harr flag and width"

> "Detect cryptic translation events with PRICE"

> "I have Ribo-RET bacterial data - how do I map initiation?"

### Differential

> "Compare initiation-site usage between treatment and control with tisdiff"

> "Which start codons change under stress?"

## What the Agent Will Do

1. QC each library with ribotish quality and write per-length offsets
2. Run ribotish predict with the elongation and initiation-drug BAMs
3. Enable alternative start codons for near-cognate detection
4. Report initiation sites with start-codon identity and ORF type
5. Optionally compare initiation across conditions

## Tips

- **Initiation drugs differ** - LTM gives sharper peaks than harringtonine
- **Pair the libraries** - predict needs both elongation (-b) and TIS (-t) BAMs
- **Enable alternative starts** - uORFs often initiate at CUG/GUG; --alt recovers them
- **Report the start codon** - alternative N-terminal and uORF starts are often non-AUG
- **Bacteria use Ribo-RET** - retapamulin, not harringtonine/LTM
- **Initiation is not stable protein** - a mapped uORF start needs downstream validation

## Related Skills

- orf-detection - Call and validate ORF bodies downstream of starts
- ribosome-periodicity - Calibrate P-site offsets
- riboseq-preprocessing - Align both libraries
- ribosome-stalling - Initiation drugs are not for elongation pausing
