# miRge3 Analysis - Usage Guide

## Overview

Quantify known miRNAs, isomiRs, tRFs, and A-to-I editing quickly with miRge3.0, which Bowtie-aligns collapsed reads to curated miRBase or MirGeneDB libraries rather than doing genome-wide discovery. This makes it the default tool for a routine differential-expression study on one of the six supported species (human, mouse, rat, zebrafish, nematode, fruitfly), but means it cannot help on an unsupported organism. Two judgments matter: choosing miRBase (large, permissive) versus MirGeneDB (small, curated) as the reference, and deciding whether to collapse isomiRs to the parent miRNA or keep seed-shifting 5' isomiRs separate. RPM output is for display only; raw counts go to DESeq2/edgeR for testing.

## Prerequisites

```bash
pip install mirge3
# bowtie (v1), RNAfold (ViennaRNA), and samtools are runtime dependencies

# Libraries are NOT fetched by miRge3 itself - download and extract them from SourceForge:
wget https://sourceforge.net/projects/mirge3/files/miRge3_Lib/human.tar.gz
tar -xzf human.tar.gz
```

## Quick Start

Tell your AI agent:
- "Quantify human miRNAs with miRge3 using MirGeneDB"
- "Emit isomiR results as mirGFF3 and detect A-to-I editing"
- "My organism is not one of the six supported - what are my options?"
- "Give me raw counts for DESeq2, not RPM"
- "Summarize isomiR diversity but keep 5' isomiRs separate"

## Example Prompts

### Basic Quantification

> "Run miRge3 annotate on my trimmed FASTQs with the human miRBase library"

> "Quantify with MirGeneDB instead of miRBase and tell me how the feature count differs"

> "Build a count matrix across all my samples for differential expression"

### IsomiR and Editing Analysis

> "Emit mirGFF3 isomiR output and classify variants as 5' versus 3'"

> "Should I collapse isomiRs to the parent miRNA for my DE analysis?"

> "Detect A-to-I editing sites and report editing frequencies per miRNA"

### Reference and Species Choices

> "My samples are rat - is rat a supported species and which database should I use?"

> "How do I build a custom library for an unsupported organism?"

## What the Agent Will Do

1. Confirm the organism is among the six supported species and locate the extracted library
2. Choose miRBase or MirGeneDB based on whether recall or curation matters, and pin the version
3. Run `miRge3.0 annotate` with the adapter (name or sequence), `-gff` for isomiR mirGFF3, and `-ai` for A-to-I editing
4. Load the raw `miR.Counts.csv` (not RPM) and filter near-zero miRNAs
5. Aggregate isomiRs deliberately - collapsing 3' variants but keeping 5' seed-shifting variants separate when isomiR identity is the question

## Tips

- miRge3 is faster than miRDeep2 for known-miRNA quantification because it aligns to small curated libraries, not the genome
- There is no `--isomir` flag - isomiR counts are produced by default and `-gff` emits the mirGFF3 standard format
- There is no `--download-library` command - libraries come from SourceForge and are extracted with `tar`
- miRge3.0 has no documented Python API; orchestrate it with `subprocess`, not `from mirge3.annotate import ...`
- Only six species ship libraries - for anything else build a custom library with miRge3_build, or use miRDeep2/sRNAbench
- miRBase (large, permissive) versus MirGeneDB (small, curated) genuinely changes the result; report which one you used
- Feed RAW counts to DESeq2/edgeR; RPM is for display only
- 5' isomiRs shift the seed and retarget the miRNA - do not silently collapse them into the canonical sequence
- Low-count 3'/internal isomiRs are often sequencing/ligation artifacts - filter hard and require replicate or UMI support; trust 5' isomiRs more
- A germline seed SNP (polymiR) looks like an isomiR/edit; fold genotypes into the reference (OptimiR) when available
- A seed A-to-I edit retargets the miRNA; keep A-to-I detection on, or permissive alignment silently merges edited reads into the canonical count
- Report and quantify by arm (-5p/-3p), never sum the two arms - the dominant arm switches across tissue/condition and the arms have different seeds and targets

## Related Skills

- smrna-preprocessing - Adapter and UMI handling
- mirdeep2-analysis - Novel miRNA discovery when needed
- differential-mirna - DE from the raw count matrix
- trf-pirna-profiling - Deeper tRF and piRNA analysis
