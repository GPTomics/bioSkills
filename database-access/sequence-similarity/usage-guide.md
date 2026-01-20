# Sequence Similarity - Usage Guide

## Overview

This skill enables AI agents to help you find homologous sequences using advanced methods including PSI-BLAST, HMMER, and reciprocal best hits.

## Prerequisites

```bash
# BLAST+
conda install -c bioconda blast

# HMMER
conda install -c bioconda hmmer

# OrthoFinder (for ortholog identification)
conda install -c bioconda orthofinder
```

## Quick Start

Tell your AI agent what you want to do:

- "Run PSI-BLAST to find distant homologs of my protein"
- "Search Pfam domains in my sequence using HMMER"
- "Find orthologs between mouse and human proteomes"
- "Build a profile HMM from my alignment"

## Example Prompts

### PSI-BLAST Searches
> "Run PSI-BLAST for 5 iterations against nr to find remote homologs of this protein"

### HMMER Domain Search
> "Search my protein against Pfam to identify conserved domains"

### Reciprocal Best Hits
> "Find 1:1 orthologs between human and mouse using reciprocal best BLAST"

### Profile Building
> "Build an HMM profile from my multiple sequence alignment and search against a proteome"

## What the Agent Will Do

1. Select the appropriate tool (BLAST, PSI-BLAST, HMMER) based on sensitivity needs
2. Configure search parameters (iterations, E-values, databases)
3. Run the similarity search
4. Parse and filter results by significance
5. Interpret homology relationships (orthologs, paralogs, domains)

## Method Comparison

| Method | Speed | Sensitivity | Use Case |
|--------|-------|-------------|----------|
| BLASTP | Fast | Moderate | Close homologs |
| PSI-BLAST | Medium | High | Remote homologs |
| HMMER | Slow | Highest | Protein families |
| Delta-BLAST | Medium | High | Domain-aware |

## E-value Guidelines

| E-value | Conclusion |
|---------|------------|
| < 1e-50 | Strong homology |
| 1e-50 to 1e-10 | Significant |
| 1e-10 to 1e-3 | Marginal |
| > 0.01 | Not significant |

## Tips

- Run PSI-BLAST for 3-5 iterations for best sensitivity
- Use lower inclusion E-value (0.001) for cleaner profiles
- HMMER is better for very distant relationships
- Always validate with reciprocal searches for ortholog calls
- Use hmmscan for domain identification, hmmsearch for protein family members
