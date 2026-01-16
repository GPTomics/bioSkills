# Entrez Fetch - Usage Guide

## Overview

This skill enables AI agents to help you retrieve records from NCBI databases. It covers fetching sequences, GenBank records, PubMed abstracts, and document summaries.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Download the FASTA sequence for NM_007294"
- "Fetch the GenBank record for this accession and show me the CDS features"
- "Get the PubMed abstract for PMID 35412348"
- "Get summary information for these 10 nucleotide IDs"

## Example Prompts

### Download Sequences
> "Download the FASTA sequences for these accessions: NM_007294, NM_000059"

### GenBank Records
> "Fetch the GenBank record for NC_000001 and list all gene features"

### PubMed
> "Get the abstracts for these PubMed IDs: 35412348, 34502548"

### Metadata Only
> "Get the lengths and organisms for these 50 nucleotide records without downloading full sequences"

### Search and Fetch
> "Find all human insulin mRNA sequences and download them as FASTA"

## What the Agent Will Do

1. Import Bio.Entrez and set up authentication
2. Determine the appropriate return type (fasta, gb, abstract, etc.)
3. Call efetch() or esummary() with correct parameters
4. Parse results with SeqIO if needed
5. Save to file or display as requested

## Output Formats

- **FASTA** - Plain sequence with header
- **GenBank** - Full annotations and features
- **Abstract** - PubMed article abstracts
- **XML** - Structured data for complex parsing

## Tips

- Provide accession numbers when you have them
- Specify the format you need (FASTA, GenBank, etc.)
- For metadata only, ask for "summaries" - it's faster
- Mention if you want to save to a file
- For large downloads, see the batch-downloads skill
