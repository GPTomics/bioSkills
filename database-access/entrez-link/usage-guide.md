# Entrez Link - Usage Guide

## Overview

This skill enables AI agents to help you navigate between NCBI databases, finding related records across different data types (genes, proteins, sequences, publications).

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Find all proteins encoded by gene ID 672"
- "What publications mention this sequence?"
- "Get the gene record for accession NM_007294"
- "Find 3D structures for BRCA1 protein"

## Example Prompts

### Gene to Protein
> "Get the RefSeq proteins for human BRCA1 (gene ID 672)"

### Sequence to Gene
> "What gene does NM_007294 belong to?"

### Related Publications
> "Find PubMed articles related to PMID 35412348"

### Multi-Database Navigation
> "Start with gene TP53, find its proteins, then get any 3D structures"

### Discover Links
> "What other databases can I link to from a gene record?"

## What the Agent Will Do

1. Import Bio.Entrez and set up authentication
2. Use elink() to find cross-references
3. Parse the LinkSetDb results
4. Chain multiple links if navigating several databases
5. Return the linked record IDs for further use

## Common Link Paths

- Gene -> Protein: Find encoded proteins
- Gene -> Nucleotide: Find mRNA sequences
- Nucleotide -> Gene: Find the gene a sequence belongs to
- Any -> PubMed: Find related publications
- Protein -> Structure: Find 3D structures

## Tips

- Provide the record ID when you have it
- Specify which database you're starting from
- Mention the target database or data type you need
- For complex navigations, describe the path you want
- Remember that elink returns IDs - use efetch to get actual records
