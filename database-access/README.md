# database-access

## Overview

Access NCBI and UniProt databases, download sequences, query SRA/GEO, and run BLAST searches. Often the starting point of bioinformatics workflows for fetching data before local processing.

**Tool type:** mixed | **Primary tools:** Bio.Entrez, Bio.Blast, SRA toolkit, BLAST+, UniProt REST API

## Skills

| Skill | Description |
|-------|-------------|
| entrez-search | Search NCBI databases with ESearch, EInfo, EGQuery |
| entrez-fetch | Retrieve records by ID with EFetch, ESummary |
| entrez-link | Cross-database references with ELink |
| batch-downloads | Large-scale downloads using history server and batching |
| sra-data | Download SRA sequencing reads with fasterq-dump, prefetch |
| geo-data | Query GEO expression datasets, link GEO to SRA |
| blast-searches | Remote BLAST searches via NCBI |
| local-blast | Local BLAST databases and searches with BLAST+ |
| sequence-similarity | PSI-BLAST, HMMER, reciprocal best hits for remote homologs |
| uniprot-access | Query UniProt protein database, ID mapping, batch retrieval |

## Example Prompts

- "Find PubMed articles about CRISPR published in 2024"
- "Download the GenBank record for NM_001234"
- "Download all RefSeq mRNA sequences for human TP53"
- "Find all proteins linked to this gene ID"
- "Download the FASTQ files for SRA accession SRR12345678"
- "Find GEO datasets for breast cancer RNA-seq"
- "BLAST this sequence against the nr database"
- "Set up a local BLAST database from my sequences"
- "Run blastp on my protein against the Swiss-Prot database"
- "Find remote homologs with PSI-BLAST"
- "Search Pfam with my protein sequence"
- "Find orthologs between two species"
- "Get the taxonomy ID for Escherichia coli K-12"
- "Fetch gene information for human BRCA1"
- "Find all human kinases in UniProt"
- "Download the UniProt entry for P53_HUMAN"
- "Map these gene symbols to UniProt IDs"

## Requirements

```bash
# Biopython
pip install biopython

# SRA Toolkit
conda install -c bioconda sra-tools

# BLAST+
conda install -c bioconda blast

# HMMER
conda install -c bioconda hmmer
```

## Related Skills

- **sequence-io** - Read/write sequence files after downloading
- **sequence-manipulation** - Work with downloaded sequences
- **alignment-files** - Process alignments after downloading from SRA
- **variant-calling** - Analyze variants from downloaded data
