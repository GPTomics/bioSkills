---
name: bio-data-visualization-sequence-logos
description: Build sequence logos from aligned DNA, RNA, or protein motifs using ggseqlogo (R), Logomaker (Python), or WebLogo with explicit bits vs probability encoding, background-frequency correction, custom alphabets, and multi-logo stacking. Use when visualizing motif PWMs (TF binding, splice sites, CRISPR spacers), aligned-position composition, or comparing two motif sets.
tool_type: mixed
primary_tool: ggseqlogo
---

## Version Compatibility

Reference examples tested with: ggseqlogo 0.2 (CRAN; per Wagih 2017), Logomaker 0.8+ (Python), WebLogo 3.7+ (CLI), Biopython 1.83+ (motif parsing), MEME suite 5.5+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Sequence Logos

**"Plot a sequence motif"** -> Render a per-position stack of letters whose total height encodes information content (Schneider-Stephens 1990 *Nucleic Acids Res* 18:6097) and individual letter height is proportional to base/aa frequency. The information-content encoding makes conserved positions visually tall and variable positions visually short — the visual is *the conservation profile*.

- R: `ggseqlogo::ggseqlogo` (Wagih 2017 *Bioinformatics* 33:3645)
- Python: `logomaker.Logo`
- CLI: `weblogo` (Crooks 2004 *Genome Res* 14:1188)

## The Single Most Important Modern Insight -- Bits vs Probability Are Different Visualizations

A sequence logo can encode each position as **bits** (information content) or **probability** (raw frequency). They look superficially similar; they communicate different things.

- **Bits (Schneider-Stephens 1990):** position height = `R = log2(K) − H(p)` where K=4 for DNA, H is Shannon entropy. Maximum 2 bits for DNA, 4.3 bits for protein. A fully conserved position is 2 bits; a uniform position is 0. This is the canonical motif encoding.
- **Probability:** position height = 1.0; letter height = frequency. Every position has the same total height. Cannot distinguish "conserved A" from "variable" — both can show 100% A at a position.
- **EDLogo (enrichment-depletion):** Dey et al. 2018 — uses log-odds of observed vs background, supporting depleted-residue display.

**Default to bits unless a specific reason exists otherwise.** Bits is what reviewers expect to see for a TF binding site, splice site, or CRISPR spacer composition.

## Decision Tree by Use Case

| Use case | Encoding | Background | Tool |
|----------|----------|------------|------|
| TF binding motif (JASPAR/CIS-BP PWM) | bits | uniform OR genome composition | ggseqlogo, Logomaker |
| Splice-site motif (5'SS, 3'SS) | bits | uniform | ggseqlogo |
| CRISPR sgRNA position-composition | probability | – | logomaker (custom alphabet) |
| Protein motif (kinase substrate) | bits | proteome composition | Logomaker (matrix_type='counts') |
| Alignment-conservation cartoon | bits OR probability | depends on intent | WebLogo |
| Differential motif (TF-A vs TF-B) | EDLogo log-odds | TF-B | Logomaker (matrix_type='weight') |

## ggseqlogo (R) -- Canonical Bioinformatics Default

**Goal:** Render a sequence motif as a per-position letter stack whose total height encodes information content (Schneider-Stephens 1990) and individual letter heights reflect frequency, optionally corrected for genome background.

**Approach:** Pass a PWM matrix (rows = letters, columns = positions) or vector of aligned same-length sequences to `ggseqlogo()` with `method = 'bits'` and explicit `bg_freq` for the relevant genome composition; stack multiple motifs as a named list.

```r
library(ggseqlogo)

# Input: PWM matrix (rows = positions, columns = nucleotides A/C/G/T)
# or aligned sequence vector

# From a vector of aligned sequences (same length)
seqs <- c('ATGCAA', 'ATGCAC', 'ATGCAG', 'ATGCAT', 'ACGCAA')
ggseqlogo(seqs, method = 'bits')

# From a PWM matrix (probability or counts)
pwm <- matrix(c(0.7, 0.1, 0.1, 0.1,
                0.1, 0.7, 0.1, 0.1,
                0.4, 0.1, 0.4, 0.1), ncol = 3,
              dimnames = list(c('A', 'C', 'G', 'T'), NULL))
ggseqlogo(pwm, method = 'bits')           # 'bits' OR 'probability'

# Multiple logos stacked (e.g., compare TF-A and TF-B)
ggseqlogo(list(TFA = seqs_a, TFB = seqs_b),
          method = 'bits',
          col_scheme = 'nucleotide')
```

```r
# Custom color scheme (protein motif, kinase substrate)
ggseqlogo(protein_pwm,
          method = 'bits',
          seq_type = 'aa',                      # auto-detected usually
          col_scheme = make_col_scheme(
              chars = c('S','T','Y','K','R','H','D','E','A','V','L','I','M'),
              cols  = c('#D55E00','#D55E00','#D55E00',          # phospho-acceptors
                        '#0072B2','#0072B2','#0072B2',          # basic
                        '#CC79A7','#CC79A7',                    # acidic
                        '#009E73','#009E73','#009E73','#009E73','#009E73')))  # hydrophobic
```

## Logomaker (Python) -- Most Flexible

```python
import logomaker
import pandas as pd

# Counts matrix (rows = position, columns = ACGT)
counts_df = pd.DataFrame({'A': [10, 0, 5, 8],
                          'C': [0, 8, 5, 1],
                          'G': [0, 2, 0, 0],
                          'T': [0, 0, 0, 1]})

# Convert counts -> information (bits)
ic_df = logomaker.transform_matrix(counts_df,
                                    from_type='counts',
                                    to_type='information',
                                    background=[0.25] * 4)    # uniform; pass real background for corrected IC

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6, 2))
logo = logomaker.Logo(ic_df,
                      color_scheme='classic',           # 'NajafabadiEtAl2017' for protein
                      shade_below=0.5,
                      fade_below=0.5,
                      font_name='Arial Rounded MT Bold')
logo.style_xticks(rotation=0)
logo.ax.set_ylabel('Bits')
```

```python
# Weight matrix (signed) -- enrichment vs depletion
weight_df = logomaker.transform_matrix(counts_df,
                                        from_type='counts',
                                        to_type='weight',
                                        background=genome_composition)
logo = logomaker.Logo(weight_df, color_scheme='classic',
                       flip_below=True)                  # depleted letters below axis
```

## WebLogo (CLI / web)

```bash
weblogo --format pdf --sequence-type dna \
        --color-scheme classic --units bits \
        --composition equiprobable \
        --fineprint '' \
        --size large \
        < aligned.fasta > logo.pdf
```

WebLogo (Crooks 2004) is the original; supports many formats and is scriptable. For reproducible figures, prefer ggseqlogo or Logomaker (programmatic, easier to integrate with multi-panel figures).

## Background Composition Correction

The bits encoding assumes a uniform background by default. For genome-derived motifs, the background should match the genome:

- Human genome: A=0.29, C=0.21, G=0.21, T=0.29 (approx)
- GC-rich genomes (Streptomyces): A=0.18, C=0.32, G=0.32, T=0.18

Without correction, a motif preferring GC in a genome where GC is rare overestimates information; conversely, an A-rich motif in an AT-rich genome underestimates.

```r
# ggseqlogo: pass `bg_freq`
ggseqlogo(pwm, method = 'bits',
          bg_freq = c(A = 0.29, C = 0.21, G = 0.21, T = 0.29))
```

```python
logomaker.transform_matrix(counts_df, from_type='counts', to_type='information',
                            background=[0.29, 0.21, 0.21, 0.29])
```

## Per-Method Failure Modes

### Probability encoding mistaken for bits

**Trigger:** Default `method = 'probability'` in some implementations.

**Mechanism:** Every position has total height 1; visually flat with all letters same total.

**Symptom:** Reviewer asks "why doesn't the logo show conservation gradient?"

**Fix:** Use `method = 'bits'` for the standard motif encoding.

### Background uniform when genome composition matters

**Trigger:** Uniform `bg_freq = c(0.25, 0.25, 0.25, 0.25)` for a non-uniform genome.

**Mechanism:** Information content overestimates conservation for preferred bases.

**Symptom:** Reported motif looks more conserved than it actually is.

**Fix:** Pass genome composition to `bg_freq` / `background` parameter.

### PWM rows/columns reversed

**Trigger:** Input matrix in samples-as-rows convention; logomaker expects positions-as-rows.

**Mechanism:** Logo renders the wrong dimension as "position."

**Symptom:** Logo has letter count = number of input rows, not motif length.

**Fix:** Transpose the matrix; verify with `print(matrix.shape)` before plotting.

### Custom alphabet not recognized

**Trigger:** RNA logo with `U` instead of `T`; protein logo with `J` or `Z`.

**Mechanism:** ggseqlogo and logomaker auto-detect alphabet from input; unusual characters may fail.

**Symptom:** Letters render as boxes or missing entirely.

**Fix:** Explicit `seq_type = 'rna'` (ggseqlogo) or pass custom color scheme (Logomaker).

### Aligned sequences of unequal length

**Trigger:** Vector of motif instances with different lengths.

**Mechanism:** Most tools require equal-length input.

**Symptom:** Error or only first N positions plotted.

**Fix:** Pre-align (MEME / TOMTOM) or trim to a common length.

### Logo for too few input sequences

**Trigger:** PWM from N=5 sequences plotted as if N=500.

**Mechanism:** Information content has small-N bias; even random sequences look "conserved" at N=5.

**Symptom:** Logo appears more meaningful than the input warrants.

**Fix:** Compute small-sample correction (Schneider 1986; standard in MEME); annotate N in caption; require N ≥ 20 for credible motif.

### Stacked logos with different alphabets compared

**Trigger:** Stacking a DNA logo above a protein logo for visual comparison.

**Mechanism:** Maximum information content differs (2 bits DNA vs 4.3 bits protein); y-axes are not comparable.

**Symptom:** Apparent "weaker" protein logo because of higher possible max.

**Fix:** Normalize both to fractional information (0-1) OR present separately.

## Reconciliation: When Logos Differ

| Pattern | Cause | Action |
|---------|-------|--------|
| ggseqlogo vs Logomaker show different heights | Different default backgrounds (uniform vs explicit) | Standardize background; recompute |
| WebLogo vs Logomaker differ at low-conservation positions | Small-sample correction differs | Use consistent N; report sample-corrected IC |
| JASPAR vs MEME PWM look different | JASPAR uses observed counts; MEME has Dirichlet prior | Document source; cite version |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Max IC per DNA position | 2 bits | Schneider-Stephens 1990 |
| Max IC per protein position | 4.32 bits (log2(20)) | Schneider-Stephens 1990 |
| Min N for credible motif | ≥ 20 instances; ≥ 100 ideal | Common practice |
| Small-sample correction | Schneider 1986 entropy correction | MEME default; ggseqlogo via small-N tools |
| TF binding-site length typical | 6-20 bp | Biology |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Logo flat with all positions = 1 | Probability mode | Switch to bits |
| Motif looks too conserved | Uniform bg in non-uniform genome | Pass `bg_freq` |
| Letter count = N samples not motif length | Matrix transposed | Verify shape |
| RNA U renders as box | Alphabet not recognized | `seq_type = 'rna'` |
| Logos at different scales overlaid | Different alphabets | Normalize OR separate |
| Logo from N=5 looks meaningful | Small-sample bias | Require N>=20; annotate |

## References

- Crooks GE, Hon G, Chandonia JM, Brenner SE. 2004. WebLogo: a sequence logo generator. *Genome Res* 14(6):1188-1190.
- Dey KK, Xie D, Stephens M. 2018. A new sequence logo plot to highlight enrichment and depletion. *bioRxiv*.
- Schneider TD. 1986. Information content of binding sites on nucleotide sequences. *J Mol Biol* 188(3):415-431.
- Schneider TD, Stephens RM. 1990. Sequence logos: a new way to display consensus sequences. *Nucleic Acids Res* 18(20):6097-6100.
- Tareen A, Kinney JB. 2020. Logomaker: beautiful sequence logos in Python. *Bioinformatics* 36(7):2272-2274.
- Wagih O. 2017. ggseqlogo: a versatile R package for drawing sequence logos. *Bioinformatics* 33(22):3645-3647.

## Related Skills

- chip-seq/motif-analysis - Discover the PWM that becomes the logo
- atac-seq/footprinting - Footprinting motifs to visualize
- clip-seq/clip-motif-analysis - CLIP-derived motifs
- alignment/multiple-alignment - Aligned sequences as logo input
- data-visualization/color-palettes - Custom alphabet color schemes
