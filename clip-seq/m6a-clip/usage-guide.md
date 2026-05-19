# m6A CLIP - Usage Guide

## Overview

Map N6-methyladenosine at single-nucleotide resolution. Five method families: antibody-based UV-CL (miCLIP/miCLIP2), antibody-free chemical (GLORI 2023, the new gold standard for stoichiometric m6A), enzyme-fusion editing (DART-seq APOBEC1-YTH), nanopore direct RNA (m6Anet, nanocompore, EpiNano), and peak-level MeRIP-seq. DRACH motif constrains 70-90% of sites but is not exclusive. Cross-method discordance is real: DART 44% in DRACH; GLORI vs miCLIP2 ~60% concordance. Triangulate across orthogonal methods for high-confidence sites.

## Prerequisites

```bash
conda install -c bioconda samtools bedtools
pip install m6anet nanopolish-eventalign
# miCLIP2 m6Aboost: https://github.com/ZarnackGroup/m6Aboost
# GLORI: https://github.com/liucongcas/GLORI-tools
# Bullseye DART: https://github.com/mekoulnik/Bullseye
```

## Quick Start

Tell your AI agent:
- "Map m6A with miCLIP2 + m6Aboost ML on Mettl3-KO calibrated data"
- "Run GLORI for stoichiometric m6A fraction per site"
- "DART-seq with Bullseye, subtracting APOBEC1-only control"
- "m6Anet on nanopore direct RNA at >= 20 coverage"
- "MeRIP-seq peaks with exomePeak2"
- "Triangulate miCLIP2, GLORI, m6Anet for high-confidence sites"
- "Why is DART editing in non-DRACH? Off-target APOBEC1"

## Example Prompts

### miCLIP2

> "Run miCLIP2 pipeline: eCLIP-style preprocessing + PureCLIP single-nt + m6Aboost ML"

> "Filter at m6Aboost score >= 0.5; expect 10-50k high-confidence sites"

### GLORI

> "GLORI-tools on pre/post conversion BAMs; per-A stoichiometric m6A fraction"

> "Filter at coverage >= 20 and m6A_fraction >= 0.1"

### DART-seq

> "Bullseye for DART editing; subtract APOBEC1-only control; filter for DRACH"

### m6Anet (Nanopore)

> "nanopolish eventalign on aligned BAM, then m6anet dataprep + inference"

### Reconciliation

> "Cross-reference miCLIP2 sites with GLORI - if both call, high confidence"

> "miCLIP2 calls 100k sites - far too many; apply m6Aboost"

> "DART 44% in DRACH - subtract APOBEC1 control, filter for DRACH"

## What the Agent Will Do

1. Pick method based on requirements: GLORI (stoichiometric), miCLIP2 (eCLIP-compatible), m6Anet (nanopore), MeRIP (cheap initial)
2. For miCLIP2: eCLIP preprocessing -> alignment -> single-nt CL -> m6Aboost ML
3. For GLORI: pre/post conversion BAMs -> GLORI-tools -> per-A stoichiometric fraction
4. For DART-seq: Bullseye with APOBEC1-only control; filter DRACH
5. For m6Anet: nanopolish eventalign -> m6anet inference
6. Triangulate across 2+ orthogonal methods
7. Report DRACH-context fraction separately from non-DRACH
8. Document method-specific limitations

## Tips

- **GLORI is the new (2023) gold standard.** Stoichiometric, antibody-free.
- **miCLIP2 + m6Aboost is the eCLIP-compatible method.** Trained on Mettl3 KO.
- **m6Anet for nanopore / isoform context.** AUC 0.83; needs >= 20 cov per DRACH.
- **DART-seq needs APOBEC1-only control.** 30-50% off-target without it.
- **DRACH covers 70-90% of m6A.** Filter selectively; some real sites are non-DRACH.
- **Method discordance is real.** 60-75% concordance across methods.
- **MeRIP-seq is peak-level.** Cannot pinpoint single A.
- **Mettl3 KO is the negative control.** Use it for any antibody-based method.
- **GLORI degrades long RNAs.** Use shorter conversion or other method.
- **Triangulate 2+ methods for publication.** Single method = lower confidence.

## Related Skills

- clip-seq/clip-preprocessing - miCLIP2 preprocessing
- clip-seq/clip-alignment - miCLIP2 alignment
- clip-seq/crosslink-site-detection - miCLIP2 CL detection
- clip-seq/clip-peak-calling - MeRIP-seq exomePeak2
- clip-seq/stamp-antibody-free - STAMP / DART antibody-free
- long-read-sequencing/nanopore-methylation - Native nanopore m6A
- long-read-sequencing/basecalling - dRNA basecalling
- epitranscriptomics/m6a-peak-calling - MeRIP peaks
- epitranscriptomics/m6a-differential - Differential m6A
- epitranscriptomics/m6anet-analysis - m6Anet workflow
