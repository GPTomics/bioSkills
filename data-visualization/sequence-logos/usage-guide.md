# Sequence Logos - Usage Guide

## Overview

Sequence logos visualize per-position base/aa composition for an aligned motif set. Position height encodes information content in bits (Schneider-Stephens 1990); letter height is proportional to frequency. The bits encoding is the canonical bioinformatics standard - conserved positions are visually tall, variable positions short. ggseqlogo (R) is the modern default; Logomaker (Python) is the most flexible for non-standard alphabets; WebLogo is the CLI / web original.

## Prerequisites

```r
install.packages('ggseqlogo')
```

```bash
pip install logomaker
# CLI option:
pip install weblogo
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot a sequence logo of these aligned TF binding sites with bits encoding"
- "Apply human-genome background composition to the information calculation"
- "Stack two motifs vertically to compare TF-A vs TF-B"
- "Custom color scheme for a protein motif by amino-acid class"
- "Show enrichment AND depletion using log-odds encoding"

## Example Prompts

### TF binding-site logo

> "Build a sequence logo of CTCF binding sites with bits encoding and human genome background composition (A=0.29, C=0.21, G=0.21, T=0.29). Use ggseqlogo's nucleotide color scheme."

### Splice site composition

> "Logo of 5' splice sites from 200 aligned exon-intron junctions. Bits encoding, uniform background. ggseqlogo."

### Differential motif

> "Use Logomaker weight (log-odds) encoding to show enrichment AND depletion in TF-A motif relative to TF-B background. Flip depleted letters below the axis."

### Protein kinase substrate

> "Logo of phosphorylation-site neighborhoods (15 aa around phospho-S/T). Color by amino-acid class: phospho-acceptor red, basic blue, acidic purple."

### Multi-motif comparison

> "Stack three TF logos (CTCF, REST, GATA1) vertically to compare conservation patterns. Same scale for all."

## What the Agent Will Do

1. Load aligned sequences (FASTA) or PWM (counts/probability/PSSM).
2. Verify alignment: same length for all entries; sequence type (DNA, RNA, protein).
3. Set background composition: uniform default; genome-derived if appropriate.
4. Compute information content per position (bits encoding).
5. Render with ggseqlogo / Logomaker / WebLogo per project preference.
6. Annotate N (number of input sequences) for transparency.
7. Apply CVD-safe color scheme; for proteins use functional-class coloring.

## Tips

- **Default to bits**, not probability. Bits show the conservation gradient; probability shows raw frequency with every position equal-height.

- **Background matters.** Uniform default is wrong for non-uniform genomes. Pass `bg_freq` (ggseqlogo) or `background` (Logomaker) with the actual base composition.

- **N ≥ 20 for credibility.** Below this, even random sequences look conserved (Schneider 1986 small-sample bias). Annotate N in the caption.

- **PWM orientation matters.** Both ggseqlogo and Logomaker expect positions as rows or columns differently - verify with `dim(pwm)` before plotting; transpose if needed.

- **EDLogo / weight encoding** (Dey 2018, Logomaker) supports depleted-residue display below the axis. Useful for differential motif analysis.

- **Stacking logos requires same alphabet and same scale.** DNA (max 2 bits) vs protein (max 4.3 bits) cannot be visually compared without normalization.

- **Custom alphabets** (RNA U, protein extended) need explicit `seq_type` (ggseqlogo) or custom color schemes (Logomaker).

- **For PSSM input** (signed log-odds), use Logomaker `from_type='counts'` then `to_type='weight'`.

- **WebLogo is best for batch CLI use** with `--composition equiprobable` (uniform) or `--composition <species>` for built-in genome compositions.

## Related Skills

- chip-seq/motif-analysis - Motif discovery that produces the PWM
- atac-seq/footprinting - Footprinting motifs to visualize
- clip-seq/clip-motif-analysis - CLIP-derived motifs
- alignment/multiple-alignment - Aligned sequences as input
- data-visualization/color-palettes - Custom alphabet schemes
