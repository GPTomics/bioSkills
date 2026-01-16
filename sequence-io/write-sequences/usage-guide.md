# Write Sequences - Usage Guide

## Overview

This skill enables AI agents to help you write biological sequence data to files in various formats using Biopython.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Create a FASTA file with these sequences"
- "Save the modified sequences to a new file"
- "Convert these records to GenBank format"
- "Write FASTQ output with quality scores"

## Example Prompts

### Basic Writing
> "Create a FASTA file with sequences for gene1, gene2, and gene3"

### Modifying and Saving
> "Read input.fasta, uppercase all sequences, and save to output.fasta"

### Format-Specific
> "Write these sequences to GenBank format with organism annotation"

### Appending
> "Add this new sequence to my existing FASTA file"

## What the Agent Will Do

1. Import Bio.SeqIO and related modules
2. Create SeqRecord objects with required attributes
3. Add format-specific annotations if needed (molecule_type for GenBank, quality scores for FASTQ)
4. Write records using SeqIO.write()
5. Verify the output was created successfully

## Tips

- FASTA is the simplest format - just needs id and sequence
- FASTQ requires quality scores for each base
- GenBank/EMBL require molecule_type annotation
- Use file handles with 'a' mode to append to existing files
