# Local BLAST - Usage Guide

## Overview

This skill enables AI agents to help you run BLAST searches locally using BLAST+ command-line tools, enabling fast unlimited searches against custom or downloaded databases.

## Prerequisites

Install BLAST+:

```bash
# macOS
brew install blast

# Ubuntu/Debian
sudo apt install ncbi-blast+

# conda
conda install -c bioconda blast
```

Verify:
```bash
blastn -version
```

## Quick Start

Tell your AI agent what you want to do:

- "Create a BLAST database from my reference sequences"
- "BLAST my query against the local database"
- "Set up a local copy of the nt database"
- "Run reciprocal best BLAST between two species"

## Example Prompts

### Create Database
> "Make a nucleotide BLAST database from references.fasta"

### Run Search
> "BLAST query.fasta against my_db with 8 threads"

### Custom Output
> "Run BLASTP and output tabular format with query ID, subject ID, percent identity, and E-value"

### Large Scale
> "I have 10,000 sequences to BLAST - set up an efficient workflow"

### Best Hits
> "Get only the best hit for each query sequence"

## What the Agent Will Do

1. Create BLAST databases with makeblastdb
2. Run appropriate BLAST program (blastn, blastp, etc.)
3. Configure output format (tabular, XML, etc.)
4. Set performance options (threads, e-value)
5. Parse and filter results as needed

## When to Use Local vs Remote BLAST

**Use Local When:**
- Many queries to run
- Need fast turnaround
- Custom database
- Offline access needed

**Use Remote When:**
- One-off searches
- Need latest NCBI databases
- Don't want to install software

## Tips

- Use -num_threads for faster searches
- -outfmt 6 gives easy-to-parse tabular output
- -max_target_seqs limits hits per query
- Create database with -parse_seqids to enable sequence extraction
- For large databases, consider downloading pre-built NCBI databases
