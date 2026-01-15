# Read Sequences - Usage Guide

## Overview

This skill enables AI agents to help you read biological sequence data from common file formats using Biopython. It covers FASTA, FASTQ, GenBank, and 40+ other formats.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Read the sequences from my FASTA file"
- "Load this GenBank file and show me the gene features"
- "Count how many sequences are in this FASTQ file"
- "Get the sequence with ID 'NM_001234' from this large FASTA"

## Example Prompts

### Basic Reading
> "Parse sequences.fasta and print each sequence ID and length"

### Working with GenBank
> "Read the GenBank file and extract all CDS product names"

### FASTQ Quality Analysis
> "Load reads.fastq and calculate the average quality score for each read"

### Large File Handling
> "I have a 10GB FASTA file and need to extract just one sequence by ID without loading the whole file"

## What the Agent Will Do

1. Import Bio.SeqIO
2. Choose the appropriate function (parse, read, index, or to_dict)
3. Specify the correct format string for your file type
4. Iterate or access records as needed
5. Extract the requested information from SeqRecord objects

## Supported File Types

The agent can read these common formats:
- **FASTA** (.fasta, .fa, .fna, .faa)
- **FASTQ** (.fastq, .fq)
- **GenBank** (.gb, .gbk)
- **EMBL** (.embl)
- **Swiss-Prot** (.dat)
- **Alignment formats** (PHYLIP, Clustal, Stockholm)

## Tips

- Always provide the file path when asking the agent to read sequences
- Mention if your file is compressed (.gz) so the agent handles it correctly
- For large files, mention the size so the agent chooses memory-efficient methods
- Specify what information you need (IDs, sequences, annotations, features)
