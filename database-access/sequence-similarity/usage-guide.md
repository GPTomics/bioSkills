# Sequence Similarity Usage Guide

This guide covers advanced methods for finding homologous sequences.

## When to Use Each Tool

### Standard BLAST
- First-pass search
- Close homologs
- Fast results

### PSI-BLAST
- Remote homologs
- Protein family members
- When BLAST finds few hits

### HMMER
- Highest sensitivity
- Profile-based search
- Domain identification

### Reciprocal Best Hit
- Ortholog identification
- Comparative genomics
- Protein function transfer

## Method Comparison

| Method | Speed | Sensitivity | Use Case |
|--------|-------|-------------|----------|
| BLASTP | Fast | Moderate | Close homologs |
| PSI-BLAST | Medium | High | Remote homologs |
| HMMER | Slow | Highest | Protein families |
| Delta-BLAST | Medium | High | Domain-aware |

## Requirements

```bash
# BLAST+
conda install -c bioconda blast

# HMMER
conda install -c bioconda hmmer

# OrthoFinder
conda install -c bioconda orthofinder
```

## E-value Guidelines

| E-value | Conclusion |
|---------|------------|
| < 1e-50 | Strong homology |
| 1e-50 to 1e-10 | Significant |
| 1e-10 to 1e-3 | Marginal |
| > 0.01 | Not significant |

## Tips

- Run PSI-BLAST for 3-5 iterations
- Use lower inclusion E-value (0.001) for cleaner profiles
- HMMER is better for very distant relationships
- Always validate with reciprocal searches
