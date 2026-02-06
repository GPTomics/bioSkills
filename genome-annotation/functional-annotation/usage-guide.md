# Functional Annotation - Usage Guide

## Overview

Assign GO terms, KEGG orthologs, Pfam domains, and EC numbers to predicted protein sequences using eggNOG-mapper (orthology-based) and InterProScan (domain-based). Combining both tools maximizes functional annotation coverage.

## Prerequisites

```bash
# eggNOG-mapper
conda install -c bioconda eggnog-mapper

# Download eggNOG database (~44 GB for full, ~9 GB for DIAMOND only)
download_eggnog_data.py --data_dir /path/to/eggnog_db -y

# InterProScan
conda install -c bioconda interproscan

# Python utilities
pip install pandas biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Add functional annotations to my predicted proteins"
- "Run eggNOG-mapper on my protein FASTA"
- "Find GO terms and KEGG pathways for my gene predictions"

## Example Prompts

### eggNOG-mapper

> "Run eggNOG-mapper on my predicted proteins with the bacterial taxonomic scope"

> "Annotate my proteins with GO terms, KEGG, and Pfam using eggNOG-mapper"

### InterProScan

> "Run InterProScan on my proteins to find Pfam domains and signal peptides"

> "Search my proteins against Pfam, TIGRFAM, and CDD databases"

### Combined Analysis

> "Run both eggNOG-mapper and InterProScan and merge the results"

> "How many of my proteins have functional annotations?"

## What the Agent Will Do

1. Verify input protein FASTA quality
2. Run eggNOG-mapper with appropriate taxonomic scope
3. Run InterProScan for domain-based annotation
4. Merge results from both tools
5. Report annotation coverage statistics
6. Export combined annotation table

## Tips

- **Run both tools** - eggNOG-mapper and InterProScan are complementary; merging improves coverage by 10-20%
- **Taxonomic scope** - Restricting to the correct taxon improves specificity
- **GO evidence** - Use `--go_evidence non-electronic` for higher quality GO terms
- **Database version** - Document the eggNOG database version (v5.0.x) for reproducibility
- **Memory** - InterProScan is memory-intensive; allocate 8+ GB RAM
- **Annotation rate** - Expect 60-80% coverage for well-studied organisms, lower for novel taxa

## Related Skills

- genome-annotation/prokaryotic-annotation - Bakta includes basic functional annotation
- genome-annotation/eukaryotic-gene-prediction - Produces proteins for annotation
- pathway-analysis/go-enrichment - Downstream GO enrichment analysis
- pathway-analysis/kegg-pathways - KEGG pathway mapping
