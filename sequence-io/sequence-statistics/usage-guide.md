# Sequence Statistics - Usage Guide

## Overview

This skill enables AI agents to help you calculate decision-grade statistics for sequence datasets and genome assemblies using Biopython: N50/L50, auN, NG50, length distributions, and GC content with explicit ambiguity handling.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Calculate N50, L50, and auN for my assembly"
- "Compute NG50 against a genome size of 3.1 Gb"
- "Generate a summary report for these sequences"
- "Compare contiguity across multiple FASTA files"
- "Report GC content but exclude ambiguous bases"

## Example Prompts

### Assembly QC
> "Calculate N50, L50, auN, and total bases for my genome assembly"

### Cross-Assembly Comparison
> "Compare NG50 and auN for all FASTA files in the assemblies folder using a 3.1 Gb genome size"

### GC Content
> "Report mean GC content, excluding N bases, for my contigs"

### Distribution
> "Show me the contig length distribution histogram for my sequences"

## Key Metrics

| Metric | What It Tells You |
|--------|-------------------|
| N50 | Contiguity using assembly size (not comparable across assemblies) |
| L50 | How many contigs contain half the assembled bases |
| auN | Smooth, threshold-free contiguity; preferred for ranking |
| NG50 | Contiguity using genome size; comparable across assemblies |
| NGA50 | Like NG50 but broken at misassembly points (NGA50 << NG50 = misassemblies) |
| GC% | G+C proportion; depends on the chosen ambiguity mode |

## What the Agent Will Do
1. Load sequences from the input file (gzip-aware)
2. Compute length, GC, and N-content arrays in one pass
3. Derive N50/L50, auN, and (with a genome size) NG50
4. Report a summary or a side-by-side comparison table

## Tips

- N50 measures contiguity, not correctness: a misassembly can post a large N50, so pair it with BUSCO completeness and read-backed validation.
- Prefer auN for comparing assemblies: it is smooth and every join raises it, whereas N50 can ignore real joins.
- Use NG50 (genome-size denominator), not N50, to compare assemblies of the same genome on a common baseline.
- Set `gc_fraction(seq, ambiguous=...)` explicitly: `remove` excludes Ns, `ignore` counts them in the denominator, `weighted` apportions expected GC; the modes give different answers on the same sequence.
- `gc_fraction` returns a fraction 0-1; multiply by 100 only for display (the removed `GC()` returned a percent).
- Median contig length is near-useless for assemblies (dominated by tiny contigs); rely on N50/auN instead.
- GC content varies by organism (humans ~41%, E. coli ~51%).

## Related Skills

- read-sequences - Parse sequences for statistics calculation
- batch-processing - Calculate stats across multiple files
- fastq-quality - Quality score statistics for FASTQ files
- sequence-manipulation/sequence-properties - Per-sequence GC content and properties
- alignment-files/bam-statistics - samtools stats/flagstat for alignment statistics
