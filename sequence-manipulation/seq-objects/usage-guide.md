# Seq Objects - Usage Guide

## Overview

This skill enables AI agents to help you create and manipulate biological sequence objects in Biopython. It covers the core Seq, MutableSeq, and SeqRecord classes that form the foundation of sequence analysis.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Create a Seq object from this DNA string"
- "I need a mutable sequence I can modify in place"
- "Build a SeqRecord with ID and description for writing to FASTA"
- "Copy this SeqRecord and change its ID"

## Example Prompts

### Basic Sequence Creation
> "Create a Seq object from 'ATGCGATCGATCG' and show its length"

### Mutable Sequences
> "I need to replace position 5 with a G in this sequence"

### SeqRecord Creation
> "Create a SeqRecord with ID 'gene1' and description 'Example gene' that I can write to a file"

### Annotated Records
> "Build a SeqRecord with organism annotation set to 'E. coli'"

### Batch Processing
> "Create SeqRecords from this list of sequences with sequential IDs"

## What the Agent Will Do

1. Import the appropriate classes (Seq, MutableSeq, SeqRecord)
2. Create the sequence object with the correct type
3. Set any required attributes (id, description, annotations)
4. Perform the requested operations
5. Return or display the result

## When to Use Each Object

- **Seq**: Most common. Use for reading, analyzing, and transforming sequences.
- **MutableSeq**: When you need to modify individual bases/residues in place.
- **SeqRecord**: When you need metadata (ID, description) or plan to write to files.

## Tips

- If you're writing to a file, you need SeqRecord (not just Seq)
- For GenBank format output, set `record.annotations['molecule_type']`
- When copying SeqRecords, use `deepcopy()` to avoid reference issues
- Seq objects are immutable - use MutableSeq for in-place changes
