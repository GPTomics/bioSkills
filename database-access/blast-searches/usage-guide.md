# BLAST Searches - Usage Guide

## Overview

This skill enables AI agents to help you run BLAST searches against NCBI databases using Biopython, finding similar sequences to identify unknown genes or find homologs.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "BLAST this sequence against the nr database"
- "Find what organism this DNA sequence comes from"
- "Search for proteins similar to my query sequence"
- "BLAST my sequence but only against human sequences"

## Example Prompts

### Identify a Sequence
> "Here's a DNA sequence I found - BLAST it to identify what it is"

### Find Homologs
> "Find protein homologs for this sequence in mammals"

### Specific Database
> "BLAST my sequence against RefSeq mRNA only"

### With Filters
> "BLAST this protein against SwissProt with E-value cutoff 1e-10"

## What the Agent Will Do

1. Determine the appropriate BLAST program (blastn, blastp, blastx, etc.)
2. Choose the right database
3. Submit the query to NCBI
4. Parse the XML results
5. Present the top hits with E-values and identities

## BLAST Programs

- **blastn** - DNA vs DNA
- **blastp** - Protein vs Protein
- **blastx** - DNA query translated, searched against proteins
- **tblastn** - Protein query vs translated DNA

## Tips

- Remote BLAST takes time (30 sec to several minutes)
- Specify organism to narrow results
- Lower E-value = more significant hit
- For many queries, consider local BLAST
- Save results to file if you need them later
