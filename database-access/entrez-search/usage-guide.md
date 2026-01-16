# Entrez Search - Usage Guide

## Overview

This skill enables AI agents to help you search NCBI databases using Biopython's Entrez module. It covers keyword searches, database exploration, and query building.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Search PubMed for articles about CRISPR from 2024"
- "Find all human BRCA1 mRNA sequences in GenBank"
- "How many records mention 'machine learning' across all NCBI databases?"
- "What fields can I search in the protein database?"

## Example Prompts

### Basic Searches
> "Search the nucleotide database for human insulin sequences"

### Date-Limited Searches
> "Find PubMed articles about COVID-19 published in the last 30 days"

### Complex Queries
> "Search for mouse OR rat sequences of the TP53 gene longer than 1000 bp"

### Database Discovery
> "List all available NCBI databases and their record counts"

### Field Discovery
> "What search fields are available in the SRA database?"

## What the Agent Will Do

1. Import Bio.Entrez and set up email authentication
2. Build the appropriate search query with field tags
3. Call esearch(), einfo(), or egquery() as needed
4. Parse and return the results
5. Explain query translations if the search returns unexpected results

## Common Databases

- **pubmed** - Biomedical literature
- **nucleotide** - DNA/RNA sequences
- **protein** - Protein sequences
- **gene** - Gene records
- **sra** - Sequencing data archives
- **taxonomy** - Organism taxonomy

## Tips

- Be specific about which database to search
- Mention if you need all results or just a count
- Specify date ranges when searching for recent publications
- Ask about available fields if you're unsure how to narrow your search
- For very large result sets, mention you need batch downloading
