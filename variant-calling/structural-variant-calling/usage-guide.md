# Structural Variant Calling Usage Guide

## Overview

Structural variants (SVs, typically >=50 bp: deletions, insertions, inversions, duplications, translocations) are never observed directly by a read - they are reconstructed from four orthogonal signals (discordant read pairs, split reads, read depth, local assembly). Each caller fuses a subset of those signals, and its blind spots follow directly from which it omits. The practical consequences: the real short-read tradeoff is sensitivity versus breakpoint resolution (not specificity), insertions are a physics limit rather than a tuning failure, genotyping is a different problem from discovery, and every benchmark F1 is meaningless without its matching parameters.

## Prerequisites

```bash
conda install -c bioconda manta delly smoove svaba gridss survivor truvari annotsv bcftools samtools pysam
# Long-read callers
conda install -c bioconda sniffles cutesv pbsv
```

## Quick Start

Tell your AI agent what you want to do:
- "Call structural variants from my BAM with multiple callers and keep the 2/4 consensus"
- "Run somatic SV calling on my tumor-normal pair with Manta and GRIDSS2"
- "Build a genotyped cohort SV matrix, not just a union of discovery calls"
- "Merge my population SV callsets in a sequence-aware way so allele frequencies are not inflated"
- "Benchmark my SV calls against GIAB with explicit Truvari parameters"
- "Decide whether my insertion-heavy analysis needs long reads"

## The four SV signals

| Signal | Detects | Breakpoint resolution | Blind to |
|--------|---------|-----------------------|----------|
| Discordant read pairs (RP) | DEL, INS, INV, DUP, BND | ~+/-300-500 bp (IMPRECISE) | nothing, but never base-precise |
| Split reads (SR, SA tag) | all types | ~+/-0-10 bp (PRECISE) | junctions inside repeats > read length |
| Read depth (RD) | DEL, DUP (dosage) | bin-size limited (100 bp-1 kb) | balanced events (INV, balanced BND) |
| Local assembly (AS) | all types incl. novel INS | base-precise + inserted sequence | large INS/repeats beyond short-read assembly |

## Caller selection guide

| Scenario | Recommended |
|----------|-------------|
| General germline | Manta |
| Somatic SVs (cancer) | Manta tumor-normal, or GRIDSS2 -> GRIPSS -> PURPLE -> LINX |
| Highest precision / single breakends | GRIDSS2 |
| 20-300 bp indel/SV "twilight zone" | SvABA |
| WES/panel | Manta `--exome` |
| Large cohort | DELLY site-list-then-regenotype, or force-genotype with Paragraph/GraphTyper2 |
| Simple germline pipeline | smoove (LUMPY + svtyper + duphold) |
| Insertion-heavy / repeat-mediated / phased | long reads (Sniffles2, cuteSV, pbsv) |

## Example Prompts

### Discovery
> "Call structural variants from my WGS BAM using Manta, DELLY, GRIDSS, and SvABA, then keep only calls supported by at least two callers"

> "Run Manta on my exome BAM - make sure the depth filter does not drop true SVs at high-depth targeted loci"

### Somatic
> "Run somatic SV calling on my tumor-normal pair and interpret complex rearrangements with the GRIDSS2 to GRIPSS to PURPLE to LINX chain"

### Cohort and population
> "I have per-sample SV discovery VCFs for 200 samples - build a properly genotyped population matrix instead of taking their union"

> "Merge my population SV callsets so allele frequencies are not inflated by position-only collapsing"

### Benchmarking
> "Benchmark my short-read SV calls against the GIAB HG002 truth set, report every Truvari parameter, and stratify by Tier 1 versus CMRG regions"

> "My short-read insertion recall is only 40% against a long-read truth set - is that a caller problem or a physics limit?"

## What the Agent Will Do

1. Map the analysis to the caller(s) whose fused signals fit the target SV types, and warn about each caller's blind spots (e.g. smoove cannot represent insertions).
2. Run discovery, applying mode flags that matter (`--exome` for panels, tumor-normal for somatic), and for GRIDSS route the raw breakpoint graph through GRIPSS/PURPLE/LINX rather than treating it as a callset.
3. For cohorts, force-genotype every sample at merged sites rather than unioning discovery VCFs, and choose sequence-aware Truvari over position-only SURVIVOR for allele-frequency work.
4. Filter with `ABS(SVLEN)` (the DEL sign convention) and annotate with AnnotSV.
5. Benchmark with explicit Truvari parameters, stratify by region (Tier 1 vs CMRG), and report event vs breakpoint vs genotype accuracy separately.

## Tips

- The merge parameters ARE the result: position-only SURVIVOR-1000 can inflate allele frequency up to 2.2x versus sequence-aware Truvari - report the merger and its parameters.
- Never filter on raw SVLEN: deletions carry a negative SVLEN, so `SVLEN >= 50` silently drops every DEL - always use `ABS(SVLEN)`.
- CIPOS/CIEND and the IMPRECISE flag encode the RP-vs-SR resolution of each call; use them to judge how much to trust a breakpoint before merging tightly.
- One biological inversion can appear as `<INV>` in Manta and as 2+ paired BND records in GRIDSS - a naive merger triple-counts or drops it.
- Below ~30x coverage split-read evidence thins and callers drift into the low-resolution RP-only regime.
- Run GIAB-CMRG, not just Tier 1: Tier 1 excludes the medically relevant repetitive genes, and reference false-duplications (e.g. CBS, KCNE1) cause reference-specific misses.

## Related Skills

- copy-number/cnvkit-analysis - read-depth CNV detection for dosage changes complementing junction-based SV callers
- long-read-sequencing/structural-variants - full long-read SV pipelines
- variant-calling/consensus-sequences - applying variants to a reference and phasing before haplotype extraction
- variant-calling/vcf-manipulation - view and query SV VCF files
- variant-calling/filtering-best-practices - general variant filtering principles
- variant-calling/variant-annotation - annotate SVs with functional information
