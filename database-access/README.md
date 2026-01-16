# database-access

## Overview

Access NCBI databases, download sequences, query SRA/GEO, and run BLAST searches. This category is often the starting point of bioinformatics workflows, fetching data from NCBI before local processing.

**Tool type:** mixed (Python + CLI)
**Primary tools:** Biopython Bio.Entrez, Bio.Blast, SRA toolkit, BLAST+

## Workflow

```
[database-access] <--- YOU ARE HERE
    |                  - Search NCBI databases
    |                  - Fetch sequences/records
    |                  - Download SRA reads
    |                  - Run BLAST searches
    v
Raw Reads (FASTQ) / Sequences (FASTA/GenBank)
    |
    v
[sequence-io] -----> Local file processing
    |
    v
[sequence-manipulation] --> Analysis
    |
    v
[Aligner: bwa/bowtie2/STAR]
    |
    v
[alignment-files] --> SAM/BAM processing
```

## Skills

| Skill | Description | Key Functions |
|-------|-------------|---------------|
| [entrez-search](entrez-search/) | Search NCBI databases | `Entrez.esearch()`, `Entrez.einfo()`, `Entrez.egquery()` |
| [entrez-fetch](entrez-fetch/) | Retrieve records from NCBI | `Entrez.efetch()`, `Entrez.esummary()` |
| [entrez-link](entrez-link/) | Cross-database references | `Entrez.elink()`, `Entrez.read()` |
| [batch-downloads](batch-downloads/) | Large-scale downloads | `Entrez.epost()`, `usehistory='y'`, batched `efetch()` |
| [sra-data](sra-data/) | Download SRA sequencing data | `fasterq-dump`, `prefetch`, `vdb-validate` |
| [geo-data](geo-data/) | Query GEO expression data | `Entrez.esearch(db='gds')`, `Entrez.esummary()` |
| [blast-searches](blast-searches/) | Remote BLAST searches | `NCBIWWW.qblast()`, `NCBIXML.read()` |
| [local-blast](local-blast/) | Local BLAST with BLAST+ | `blastn`, `blastp`, `makeblastdb` |

## Entrez Utilities

The E-utilities are web services providing access to NCBI databases. Common workflows:

```
Search -> Fetch (simple)
┌─────────────────┐     ┌─────────────────┐
│  Entrez.esearch │ --> │  Entrez.efetch  │ --> Records
│  (get IDs)      │     │  (get data)     │
└─────────────────┘     └─────────────────┘

Search -> Link -> Fetch (cross-database)
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│  Entrez.esearch │ --> │  Entrez.elink   │ --> │  Entrez.efetch  │ --> Records
│  (gene IDs)     │     │  (protein IDs)  │     │  (proteins)     │
└─────────────────┘     └─────────────────┘     └─────────────────┘

Search -> History -> Batch Fetch (large scale)
┌─────────────────┐     ┌─────────────────┐
│  Entrez.esearch │ --> │  Entrez.efetch  │ --> Records
│  usehistory='y' │     │  (in batches)   │     (streamed to file)
│  WebEnv/query_key     └─────────────────┘
└─────────────────┘
```

## NCBI Databases

Common databases accessible via Entrez:

| Database | ID | Description |
|----------|-----|-------------|
| PubMed | `pubmed` | Biomedical literature |
| Nucleotide | `nucleotide` | DNA/RNA sequences (GenBank/RefSeq) |
| Protein | `protein` | Protein sequences |
| Gene | `gene` | Gene records |
| SRA | `sra` | Sequence Read Archive |
| GEO | `gds` | Gene Expression Omnibus datasets |
| GEO Profiles | `geoprofiles` | Individual expression profiles |
| Taxonomy | `taxonomy` | Organism taxonomy |
| Assembly | `assembly` | Genome assemblies |
| BioProject | `bioproject` | Biological projects |
| BioSample | `biosample` | Biological samples |
| dbSNP | `snp` | Single nucleotide polymorphisms |
| ClinVar | `clinvar` | Clinical variants |
| Structure | `structure` | 3D macromolecular structures |

Use `Entrez.einfo()` to get the full list of 30+ databases.

## Authentication

NCBI requires identification for Entrez access. Set `Entrez.email` (required) and optionally `Entrez.api_key` to increase rate limits.

| Configuration | Rate Limit |
|---------------|------------|
| Email only | 3 requests/second |
| Email + API key | 10 requests/second |

Get an API key at: https://www.ncbi.nlm.nih.gov/account/settings/

## Example Prompts

Ask your AI agent naturally:

| Task | Example Prompt |
|------|----------------|
| Search PubMed | "Find PubMed articles about CRISPR published in 2024" |
| Fetch sequences | "Download the GenBank record for NM_001234" |
| Batch download | "Download all RefSeq mRNA sequences for human TP53" |
| Cross-reference | "Find all proteins linked to this gene ID" |
| SRA download | "Download the FASTQ files for SRA accession SRR12345678" |
| GEO query | "Find GEO datasets for breast cancer RNA-seq" |
| Remote BLAST | "BLAST this sequence against the nr database" |
| Local BLAST | "Set up a local BLAST database from my sequences" |
| Taxonomy lookup | "Get the taxonomy ID for Escherichia coli K-12" |
| Gene info | "Fetch gene information for human BRCA1" |

## Requirements

### Python (Biopython)
```bash
pip install biopython
```

### SRA Toolkit (for sra-data skill)
```bash
# macOS
brew install sratoolkit

# Ubuntu/Debian
sudo apt install sra-toolkit

# conda
conda install -c bioconda sra-tools
```

### BLAST+ (for local-blast skill)
```bash
# macOS
brew install blast

# Ubuntu/Debian
sudo apt install ncbi-blast+

# conda
conda install -c bioconda blast
```

## Related Skills

- **sequence-io** - Read/write sequence files after downloading from NCBI
- **sequence-manipulation** - Work with downloaded sequences (translation, GC content)
- **alignment-files** - Process alignments after downloading from SRA
- **variant-calling** - Analyze variants from downloaded data
