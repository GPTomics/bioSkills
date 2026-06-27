# miRDeep2 Analysis - Usage Guide

## Overview

Discover novel miRNAs and quantify known miRNAs with miRDeep2, which scores genome-mapped read stacks against the Dicer/Drosha biogenesis signature (a sharp mature arm, a lower-abundance star arm with ~2-nt 3' overhangs, a stable hairpin fold, and a randfold p-value). The key judgment is that a miRDeep2 score is a structural and expression hypothesis, not a validated miRNA: novel discovery is intrinsically high false-positive, and the classic false positives are tRNA and rRNA fragments whose hairpins throw scoring read stacks. There is no universal score cutoff - survey.pl reports signal-to-noise and an estimated FDR across cutoffs, and the analyst picks and reports one. Reach for miRDeep2 only when de novo discovery is actually needed; otherwise quantify known miRNAs with a faster tool.

## Prerequisites

```bash
conda install -c bioconda mirdeep2 bowtie viennarna
# miRDeep2 uses bowtie 1 (not bowtie2)

# miRBase references (pin the version)
wget https://www.mirbase.org/download/mature.fa
wget https://www.mirbase.org/download/hairpin.fa
```

## Quick Start

Tell your AI agent:
- "Discover novel miRNAs in my collapsed small RNA-seq reads with miRDeep2"
- "Set up the human genome bowtie 1 index and run mapper.pl"
- "Quantify known human miRNAs only, skipping discovery"
- "Filter my novel candidates against tRNA and rRNA loci"
- "Use survey.pl signal-to-noise to choose a score cutoff"

## Example Prompts

### Novel miRNA Discovery

> "Run miRDeep2 discovery on my collapsed FASTA with human and mouse miRBase references"

> "My novel candidates overlap tRNA loci - how do I reject those false positives?"

> "Which score cutoff should I use, and how do I justify it from the survey output?"

### Quantification

> "Quantify known human miRNAs with quantifier.pl without running discovery"

> "Why do quantifier.pl and miRDeep2.pl give different counts for the same miRNA?"

### Results Analysis

> "Parse the result CSV and keep candidates above my chosen signal-to-noise cutoff"

> "Check star-arm read support and randfold p-value for my top novel calls"

## What the Agent Will Do

1. Build the bowtie 1 genome index and confirm reads are adapter-trimmed and collapsed
2. Run mapper.pl to align collapsed reads and emit the ARF alignment file
3. Run miRDeep2.pl with same-species and related-species miRBase references for discovery and quantification (or quantifier.pl for known-only)
4. Choose a score cutoff from survey.pl signal-to-noise / FDR and report it
5. Filter novel candidates against tRNA/rRNA/snoRNA loci and require star-arm support and reproducibility before trusting them

## Tips

- A miRDeep2 score is a hypothesis, not a validated miRNA - novel calls need orthogonal validation (Northern/qPCR, Dicer/Drosha-knockdown loss of signal, conservation)
- There is no universal cutoff; the "score > 10 = high confidence" rule is folklore - use survey.pl signal-to-noise (Friedländer used the lowest cutoff giving SNR >= 5) and report your choice
- The classic false positive is a tRNA/rRNA fragment hairpin - always intersect novel candidates against structured-RNA annotations
- miRDeep2 needs a reference genome and bowtie 1 (not bowtie2), and pre-collapsed reads with the `_xN` count
- Pin the miRBase version; convert U->T and strip header whitespace in reference FASTAs
- Provide a related-species mature FASTA for conservation evidence; use 'none' only if no relative is available
- If you only need known-miRNA counts, run quantifier.pl (or mirge3) and skip the expensive discovery engine
- Judge novel calls against the community criteria (Ambros 2003), where consistent 5' processing of both arms is the single most discriminating signal - deliver a criteria table, not just a score ranking
- miRDeep2 is genome-anchored; for a non-model animal with no assembly (or single-cell), use a genome-free tool such as Mirnovo

## Related Skills

- smrna-preprocessing - Adapter trimming and read collapsing
- mirge3-analysis - Faster known-miRNA and isomiR quantification
- differential-mirna - Differential expression of the count matrix
- trf-pirna-profiling - tRF/piRNA biology that otherwise appears as false positives
- genome-annotation/ncrna-annotation - Annotate tRNA/rRNA loci to filter candidates
