# Codon Usage - Usage Guide

## Overview

This skill enables AI agents to help you analyze codon usage in coding sequences using Biopython: codon counting, CAI (Codon Adaptation Index) scoring against a host, RSCU and Nc bias metrics, GC at codon positions, and naive max-CAI codon optimization for heterologous expression. It uses the redesigned `Bio.SeqUtils.CodonAdaptationIndex` class (the old `Bio.SeqUtils.CodonUsage` module was removed in Biopython 1.82).

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Calculate codon frequencies and RSCU for this coding sequence"
- "What is the CAI of this gene against my E. coli reference genes?"
- "Build a CAI index from these highly expressed genes and score my construct"
- "Codon-optimize this gene for the host, but warn me about expression risks"
- "Find rare/low-RSCU codons in my sequence"
- "What is the GC content at each codon position?"

## Example Prompts

### CAI Against a Host Reference
> "Build a Codon Adaptation Index from this FASTA of highly expressed E. coli genes and report the CAI of my target gene."

### Codon Optimization
> "Codon-optimize this CDS for the host using max-CAI, confirm the protein is unchanged, and list which tradeoffs (ramp, 5' structure, cryptic sites) I should screen before ordering it."

### RSCU Analysis
> "Calculate RSCU values for my coding sequence and flag over- and under-used synonymous codons."

### Rare Codon Detection
> "Find the codons with RSCU below 0.5 in this gene."

### Overall Bias
> "Compute the effective number of codons (Nc) for this sequence."

### Wobble Position Bias
> "What is the GC content at the first, second, and third codon positions?"

### Mitochondrial or Bacterial Genes
> "Score this mitochondrial gene's CAI using the vertebrate mitochondrial codon table."

## What the Agent Will Do

1. Import `CodonAdaptationIndex` from `Bio.SeqUtils` (not the removed `CodonUsage` module)
2. For CAI: build the index from a reference set of highly expressed host genes, then call `calculate`
3. Confirm the query is in frame before scoring (out-of-frame input is silently mis-scored)
4. For optimization: call `optimize`, verify the protein is preserved, and flag expression-risk tradeoffs
5. Compute requested bias metrics (RSCU, Nc, GC123) and return formatted results

## Key Concepts

### CAI (Codon Adaptation Index)
- Geometric mean of per-codon relative adaptiveness (w), range 0-1, higher = better adapted
- w must come from highly expressed genes of the TARGET organism; CAI against a whole-genome or wrong-organism reference is meaningless
- ATG, TGG, and stop codons are excluded; unobserved codons get w = 0.5

### RSCU (Relative Synonymous Codon Usage)
- Observed count / expected-if-uniform within a synonymous family
- 1.0 = no bias, >1 over-used, <1 under-used; the w weights in CAI are built from RSCU ratios

### Nc (Effective Number of Codons)
- Reference-free overall bias measure, range ~20 (fully biased) to 61 (unbiased)

### tAI (tRNA Adaptation Index)
- Supply-side alternative weighting codons by tRNA gene copy number and wobble efficiency
- Often tracks expression better than CAI, but is not in Biopython (use the R `tAI` package)

### Why naive max-CAI optimization can hurt expression
- `optimize()` only maximizes single-codon CAI; it ignores the translation ramp, 5' mRNA structure, GC extremes, cryptic regulatory elements, and codon-pair bias
- Failure is silent: correct protein, high CAI, poor expression. Treat the output as a draft to screen.

## Tips

- The biggest trap: `from Bio.SeqUtils import CodonUsage` raises ImportError on Biopython >= 1.82. Use `from Bio.SeqUtils import CodonAdaptationIndex`.
- There is no bundled Sharp E. coli index. Always build the index from your own reference CDS.
- Keep queries in frame (length divisible by 3, first base = codon 1); out-of-frame input does not raise, it is silently corrupted.
- Pass the matching `table=` for mitochondrial or bacterial genes.
- `optimize(strict=True)` raises on a tie; use `strict=False` to let Biopython pick one.
- GC123 returns percentages (0-100); `gc_fraction` returns a 0-1 fraction. GC3 correlates with genome GC.

## Related Skills

- transcription-translation - Translate CDS and select the correct codon table
- sequence-properties - GC123 and per-position GC content
- sequence-io/read-sequences - Parse reference CDS from FASTA/GenBank for CAI training
- database-access/entrez-fetch - Fetch highly expressed gene sets from NCBI for CAI references
