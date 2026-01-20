# Protein Inference - Usage Guide

## Overview
Determine which proteins are present in a sample from identified peptides, handling shared peptides between homologous proteins through grouping and parsimony.

## Prerequisites
```bash
pip install pyopenms pandas
# CLI: ProteinProphet (TPP), EPIFANY (OpenMS)
# R alternative: BiocManager::install("MSnbase")
```

## Quick Start
Tell your AI agent what you want to do:
- "Group proteins by shared peptide evidence from my search results"
- "Apply parsimony to find the minimum protein set"
- "Filter to proteins with at least 2 unique peptides"

## Example Prompts

### Protein Grouping
> "Parse MaxQuant proteinGroups.txt and explain the protein group structure"

> "Identify protein groups that share all peptide evidence (indistinguishable)"

> "List proteins with unique peptides vs those identified only by shared peptides"

### Inference Methods
> "Apply parsimony principle to report minimum protein set explaining all peptides"

> "Run EPIFANY for probabilistic protein inference with FDR control"

> "Use Occam's razor to assign shared peptides to the most likely protein"

### Quality Filtering
> "Filter to proteins with at least 2 unique peptides for confident identification"

> "Apply 1% protein-level FDR and report the number of protein groups"

> "Identify single-peptide hits and flag them for review"

### Reporting
> "Create a protein list with gene names, unique peptides, and coverage"

> "Export protein groups in mzIdentML format"

> "Summarize how many proteins are identified at each evidence level"

## What the Agent Will Do
1. Load peptide-protein mappings from search results
2. Build protein groups based on shared peptides
3. Apply inference method (parsimony, probabilistic)
4. Calculate protein-level FDR
5. Filter by unique peptides and FDR threshold
6. Generate protein list with evidence summary

## Inference Strategies

| Method | Description |
|--------|-------------|
| Parsimony | Minimum protein set explaining all peptides |
| Occam's Razor | Assign shared peptides to protein with most evidence |
| Probabilistic | ProteinProphet, EPIFANY - probability-based |
| All peptides | Report all possible proteins (most inclusive) |

## Key Concepts
- **Unique peptides**: Map to only one protein (strongest evidence)
- **Razor peptides**: Shared peptides assigned to winning protein
- **Protein groups**: Proteins with indistinguishable evidence

## Tips
- Require >= 2 unique peptides for confident identification
- Protein FDR is separate from peptide FDR (typically 1-5%)
- Report protein groups, not just lead proteins
- Be cautious with single-peptide identifications
- Document which inference method was used
