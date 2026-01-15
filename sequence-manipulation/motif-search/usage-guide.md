# Motif Search - Usage Guide

## Overview

This skill enables AI agents to help you find patterns, motifs, and subsequences in biological sequences using Biopython. It covers exact matches, regex patterns, IUPAC ambiguity codes, and probabilistic motif searching.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Find all occurrences of GAATTC in this sequence"
- "Search for EcoRI and BamHI restriction sites"
- "Find TATA box variants in the promoter region"
- "Search both strands for this motif"

## Example Prompts

### Exact Pattern Search
> "Find all positions where GAATTC occurs in my sequence"

### Restriction Sites
> "What restriction enzymes cut this sequence?"

### Both Strands
> "Search for GAATTC on both strands"

### Ambiguous Pattern
> "Find all matches to the pattern GATNNTC where N is any base"

### Regulatory Elements
> "Search for TATA box variants in the first 500 bp"

### Repeats
> "Find tandem repeats of CAG in this sequence"

### Start Codons
> "Find all ATG start codons and their positions"

## What the Agent Will Do

1. Import Bio.Seq and regex if needed
2. Create appropriate search function for the pattern type
3. Search the sequence (and reverse complement if requested)
4. Return positions and matched subsequences
5. Format results clearly

## Pattern Types

- **Exact**: Fixed sequence like `GAATTC`
- **IUPAC**: With ambiguity codes like `GATNNTC` (N = any base)
- **Regex**: Flexible patterns like `TATA[AT]A[AT]`
- **PWM/PSSM**: Probabilistic scoring matrices

## Tips

- Always search both strands for biological relevance
- Use IUPAC codes for patterns with known ambiguity
- For regulatory elements, use regex for variant matching
- Consider frame position for coding sequence motifs
- Large sequences benefit from compiled regex patterns
