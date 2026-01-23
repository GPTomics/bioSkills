# Format Conversion - Usage Guide

## Overview

This skill enables AI agents to help you convert between biological sequence file formats using Biopython.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Convert this GenBank file to FASTA"
- "Convert all .gb files in this folder to FASTA"
- "Change my FASTQ to FASTA format"
- "Convert FASTA to GenBank with DNA molecule type"

## Example Prompts

### Simple Conversion
> "Convert sequence.gb to FASTA format"

### Batch Conversion
> "Convert all GenBank files in data/ to FASTA"

### With Modifications
> "Convert to FASTA and uppercase all sequences"

## What the Agent Will Do
1. Parse sequences from input format
2. Convert to target format representation
3. Write sequences in new format
4. Preserve sequence metadata where applicable

## Tips

- GenBank to FASTA is always safe (just loses annotations)
- FASTA to GenBank requires adding molecule_type
- FASTA to FASTQ requires adding quality scores
- Use SeqIO.convert() for simple conversions, parse/write for complex ones
