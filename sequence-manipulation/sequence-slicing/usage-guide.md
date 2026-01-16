# Sequence Slicing - Usage Guide

## Overview

This skill enables AI agents to help you extract, slice, and concatenate biological sequences using Biopython. It covers indexing, slicing, joining sequences, and working with coordinate systems.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Extract positions 100-200 from this sequence"
- "Get the first 50 bases of each sequence"
- "Join these sequences with an NNN linker"
- "Split this sequence into codons"

## Example Prompts

### Basic Extraction
> "Extract nucleotides 100 to 500 from my sequence"

### With 1-Based Coordinates
> "The gene starts at position 1234 and ends at 5678 (1-based). Extract it."

### Split into Parts
> "Split this sequence into 100 bp chunks"

### Codon Extraction
> "Show me each codon in this coding sequence"

### Concatenation
> "Join these three exon sequences together"

### With Linker
> "Concatenate these sequences with NNNNNN between each"

### Flanking Regions
> "Get 50 bp upstream and downstream of position 1000"

## What the Agent Will Do

1. Import Bio.Seq
2. Create or use existing Seq objects
3. Apply Python slicing syntax with proper indexing
4. Handle coordinate system conversions if needed
5. Return the extracted or combined sequence

## Coordinate System Notes

Python uses 0-based indexing with exclusive end positions:
- `seq[0]` is the first base
- `seq[0:3]` gets bases at positions 0, 1, 2 (three bases total)
- `seq[100:200]` gets 100 bases (positions 100-199)

Biology often uses 1-based coordinates:
- Position 1 is the first base
- Position 100-200 means 101 bases
- Convert: `seq[start-1:end]` for 1-based coordinates

## Tips

- Always clarify if coordinates are 0-based or 1-based
- Slicing returns Seq objects, not strings
- Use `+` to concatenate, or `linker.join(seqs)` for multiple
- For GenBank features, use `feature.extract(record.seq)`
- Check sequence length before slicing to avoid errors
