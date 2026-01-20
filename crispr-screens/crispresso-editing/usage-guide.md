# CRISPResso2 Editing Analysis - Usage Guide

## Overview
CRISPResso2 analyzes amplicon sequencing data from CRISPR editing experiments. It quantifies indels, HDR efficiency, and characterizes editing outcomes.

## Prerequisites
```bash
pip install crispresso2
# or via conda
conda install -c bioconda crispresso2
```

## Quick Start
Tell your AI agent what you want to do:
- "Analyze my CRISPR editing experiment to quantify indel rates"
- "Calculate HDR efficiency from my amplicon sequencing"
- "Compare editing outcomes between samples"

## Example Prompts

### Basic Editing Analysis
> "I have amplicon sequencing data from a CRISPR knockout experiment. Run CRISPResso to quantify the indel rate at my target site."

> "Analyze my CRISPR editing FASTQ files. The amplicon is 300bp and my guide sequence is ACGTACGTACGTACGTACGT."

### HDR Analysis
> "I'm doing HDR with a repair template. Run CRISPResso with my expected HDR sequence to calculate precise editing efficiency."

> "Quantify knock-in efficiency from my amplicon sequencing. I have the wild-type and expected HDR amplicon sequences."

### Batch and Comparison
> "I have 20 samples all targeting the same site. Run CRISPRessoBatch to analyze them together."

> "Compare editing efficiency between my control and treated samples using CRISPRessoCompare."

### Base and Prime Editing
> "Analyze my base editing experiment. I'm looking for C>T conversions at position 6 of the guide."

> "Run CRISPResso on my prime editing samples and quantify the precise edit rate."

## What the Agent Will Do
1. Identify appropriate CRISPResso mode (single, batch, pooled)
2. Set up amplicon and guide sequences
3. Run alignment and quantification
4. Calculate editing rates (NHEJ, HDR, other)
5. Generate indel size distribution plots
6. Create comparison reports if multiple samples

## Tips
- Minimum 1000x read depth recommended for reliable quantification
- Check mapping percentage (should be >80%)
- Use --exclude_bp_from_left/right to ignore primer artifacts
- For HDR, always provide --expected_hdr_amplicon_seq
- Review the allele table for detailed editing outcomes

## Analysis Types

| Tool | Use Case |
|------|----------|
| CRISPResso | Single amplicon |
| CRISPRessoBatch | Multiple samples, same amplicon |
| CRISPRessoPooled | Multiple amplicons per sample |
| CRISPRessoWGS | Whole genome off-target |
| CRISPRessoCompare | Compare conditions |

## References
- CRISPResso2: doi:10.1038/s41587-019-0032-3
- Documentation: https://crispresso.pinellolab.partners.org/
