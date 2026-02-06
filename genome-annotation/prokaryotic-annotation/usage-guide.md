# Prokaryotic Annotation - Usage Guide

## Overview

Annotate bacterial and archaeal genomes with structural and functional gene predictions using Bakta (preferred) or Prokka (legacy). Produces GFF3, GenBank, protein FASTA, and other standard output formats ready for downstream analysis or NCBI submission.

## Prerequisites

```bash
# Bakta (recommended)
conda install -c bioconda bakta

# Download database (~30 GB for full, ~1.5 GB for light)
bakta_db download --output /path/to/bakta_db --type full

# Prokka (legacy alternative)
conda install -c bioconda prokka

# Python parsing
pip install gffutils biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate my bacterial genome assembly with Bakta"
- "Run Prokka on my E. coli assembly"
- "Prepare my genome annotation for NCBI submission"

## Example Prompts

### Basic Annotation

> "Annotate this bacterial assembly with Bakta using the full database"

> "Run Bakta on my Staphylococcus genome with locus tag prefix SAUR"

### With Organism Info

> "Annotate my E. coli K-12 genome with Bakta including genus, species, and strain info"

> "Run Prokka on my Pseudomonas assembly with the genus-specific database"

### QC and Submission

> "Check the annotation quality - how many genes, coding density, hypothetical proteins"

> "Prepare my Bakta annotation for NCBI GenBank submission"

## What the Agent Will Do

1. Verify assembly quality and check for appropriate database
2. Run Bakta (or Prokka) with organism-specific parameters
3. Generate GFF3, GenBank, protein FASTA, and summary outputs
4. Report annotation statistics (gene count, coding density, tRNA/rRNA counts)
5. Flag potential issues (low coding density, many hypothetical proteins)
6. Suggest functional annotation follow-up if needed

## Tips

- **Database choice** - Use Bakta full database for comprehensive results; light database for quick screening
- **Translation table** - Default table 11 works for most bacteria; use table 4 for Mycoplasma/Spiroplasma
- **Complete genomes** - Add `--complete` flag for circularized assemblies to enable origin detection
- **Gram stain** - Specify `--gram +` or `--gram -` for improved signal peptide prediction
- **Locus tags** - Use 3-12 uppercase alphanumeric characters for NCBI compliance
- **Hypothetical proteins** - 30-50% hypothetical is normal for novel organisms; run eggNOG-mapper for more

## Related Skills

- genome-annotation/functional-annotation - Add GO/KEGG/Pfam to predicted proteins
- genome-annotation/ncrna-annotation - Detailed ncRNA identification
- genome-assembly/assembly-qc - Check assembly quality before annotation
- genome-intervals/gtf-gff-handling - Parse GFF3 output files
