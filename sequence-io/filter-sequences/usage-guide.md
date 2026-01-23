# Filter Sequences - Usage Guide

## Overview

This skill enables AI agents to help you filter and select sequences based on various criteria using Biopython.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Keep only sequences longer than 500 bp"
- "Extract sequences matching these IDs"
- "Filter out sequences with N's"
- "Select sequences with GC content between 40-60%"

## Example Prompts

### By Length
> "Remove sequences shorter than 200 bp from my FASTA file"

### By ID
> "Extract sequences with IDs listed in wanted_ids.txt"

### By Content
> "Filter out sequences containing more than 5% N bases"

### By GC Content
> "Keep only sequences with GC content above 50%"

### Combined
> "Filter sequences: length >= 100, no N's, GC 40-60%"

## What the Agent Will Do
1. Load sequences from input file
2. Apply specified filter criteria
3. Keep sequences passing all filters
4. Write filtered sequences to output

## Tips

- Use generator expressions for large files (memory efficient)
- For ID-based filtering, use sets for O(1) lookup
- Combine multiple criteria in a single filter function
- Generators can only be used once - recreate if needed
