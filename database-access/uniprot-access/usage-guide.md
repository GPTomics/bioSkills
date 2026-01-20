# UniProt Access - Usage Guide

## Overview

This skill enables AI agents to help you query UniProt programmatically using the REST API to retrieve protein sequences, annotations, and functional information.

## Prerequisites

```bash
pip install requests
```

## Quick Start

Tell your AI agent what you want to do:

- "Download the FASTA sequence for UniProt P04637"
- "Search UniProt for all human kinases"
- "Map these gene names to UniProt accessions"
- "Get functional annotations for BRCA1"

## Example Prompts

### Sequence Retrieval
> "Download the protein sequence for UniProt accession P53_HUMAN"

### Protein Search
> "Find all reviewed human proteins with 'transcription factor' in their name"

### ID Mapping
> "Convert these gene names to UniProt IDs: TP53, BRCA1, EGFR"

### Bulk Download
> "Download all Swiss-Prot entries for E. coli K-12 as FASTA"

### Specific Fields
> "Get accession, gene name, and subcellular location for all human membrane proteins"

## What the Agent Will Do

1. Construct the appropriate REST API query
2. Choose the right endpoint (search, stream, idmapping)
3. Set output format (FASTA, JSON, TSV, XML)
4. Handle pagination for large result sets
5. Parse and present the results

## Output Formats

| Format | Use Case |
|--------|----------|
| `fasta` | Sequences for alignment/BLAST |
| `json` | Programmatic parsing |
| `tsv` | Tabular analysis, pandas |
| `xml` | Full structured data |
| `txt` | Human-readable flat file |
| `gff` | Feature coordinates |

## Tips

- Use `reviewed:true` for Swiss-Prot (curated) entries only
- Use the stream endpoint for >500 results
- Specify fields in TSV format to reduce data transfer
- Be polite with requests - add delays for batch queries
- Cache results when possible - UniProt updates frequently but not constantly
