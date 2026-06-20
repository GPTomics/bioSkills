# Translation Efficiency - Usage Guide

## Overview

Quantify translation efficiency (TE = ribosome occupancy / mRNA abundance) and test for differential TE between conditions. The key decisions are counting both assays over the CDS, using count-based GLMs rather than ratio tests, and distinguishing genuine translational control from buffering.

## Prerequisites

```r
# R (differential TE)
BiocManager::install(c("riborex", "xtail", "anota2seq", "DESeq2"))
```

```bash
# Python (quick ranking screen)
pip install pandas numpy statsmodels
```

## Quick Start

Tell your AI agent what you want to do:
- "Calculate translation efficiency from my Ribo-seq and RNA-seq CDS counts"
- "Find genes with differential TE with riborex"
- "Tell translational control from buffering with anota2seq"
- "Run a DESeq2 interaction model for TE"

## Example Prompts

### Calculate TE

> "Compute per-gene log2 TE over the CDS for ranking"

> "Why do my long-UTR genes have systematically low TE?"

### Differential TE

> "Run riborex to find genes with significantly changed TE"

> "Use anota2seq to classify genes as translation, buffering, or abundance"

> "Set up a DESeq2 condition-by-assay interaction and pick the right coefficient"

### Interpretation

> "Which genes are translationally activated versus buffered?"

> "Is this TE change real translational control or secondary to a uORF?"

## What the Agent Will Do

1. Load matched Ribo-seq and RNA-seq CDS count matrices
2. For ranking, compute per-gene log2 TE with a pseudocount
3. For testing, run a count-based GLM (riborex/Xtail/anota2seq/DESeq2 interaction)
4. With anota2seq, classify each gene's mode of regulation
5. Filter by adjusted p-value and check uORF/isoform confounders

## Tips

- **Count both over the CDS** - full-transcript RNA against CDS RPF biases long-UTR genes
- **Ratios rank, GLMs test** - never t-test log-TE for inference
- **Buffering looks like control** - only anota2seq separates them
- **Do not hardcode the DESeq2 interaction name** - pick it from resultsNames(dds)
- **Per-assay size factors** - assume median TE unchanged; use spike-ins for global shifts
- **Trim codons** - exclude first ~15 / last ~5 codons from the RPF count
- **Replicates** - count GLMs and anota2seq RVM want >=3 per condition per assay; n=2 is illustrative only

## Related Skills

- ribosome-periodicity - Calibrate P-site offsets for CDS counts
- orf-detection - Rule out uORF-driven TE changes
- rna-quantification/featurecounts-counting - Matched RNA-seq CDS counts
- differential-expression/deseq2-basics - Count-based DE foundations
