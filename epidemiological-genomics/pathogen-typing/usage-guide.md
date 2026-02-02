# Pathogen Typing - Usage Guide

## Overview

Identify bacterial strain types using MLST, cgMLST, and SNP-based approaches for outbreak investigation and surveillance.

## Prerequisites

```bash
conda install -c bioconda mlst
pip install chewbbaca pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Type my Salmonella isolates with MLST"
- "Run cgMLST on these E. coli genomes"
- "Find which isolates cluster together"

## Example Prompts

### Basic MLST

> "Determine the sequence type of my bacterial genome"

> "Run MLST on all FASTAs in this directory"

### Core Genome MLST

> "Run cgMLST and identify outbreak clusters"

> "Calculate allelic distances between my isolates"

### Outbreak Investigation

> "Which isolates are within 5 allele differences?"

> "Identify transmission clusters from my typing results"

## What the Agent Will Do

1. Identify organism and select appropriate scheme
2. Run MLST/cgMLST typing
3. Parse results into structured format
4. Calculate distances between isolates
5. Identify clusters based on thresholds
6. Report sequence types and cluster assignments

## Tips

- **Scheme selection** - mlst auto-detects; specify if known
- **cgMLST thresholds** - Pathogen-specific (E. coli: 10, Salmonella: 7)
- **Novel STs** - Submit new types to PubMLST for official designation
- **Assembly quality** - Poor assemblies may give incomplete profiles
- **Mixed cultures** - MLST assumes clonal isolates
