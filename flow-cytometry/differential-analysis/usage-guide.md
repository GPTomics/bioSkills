# Differential Analysis Usage Guide

## Overview

Differential analysis identifies cell populations that differ in abundance or marker expression between experimental conditions.

## Analysis Types

### Differential Abundance (DA)
- Are cluster proportions different between conditions?
- Uses edgeR or generalized linear models
- Reports fold change in frequency

### Differential State (DS)
- Is marker expression different within clusters?
- Uses limma
- Tests each marker in each cluster

## Statistical Considerations

### Sample Size
- Need biological replicates (n ≥ 3 per group)
- Technical replicates don't count

### Normalization
- For DA: proportions or counts with TMM
- For DS: per-sample normalization

### Multiple Testing
- Correct for # clusters tested (DA)
- Correct for # clusters × # markers (DS)

## Methods

| Method | Use Case |
|--------|----------|
| edgeR | DA, count-based |
| limma | DS, continuous |
| Mixed models | Paired designs |
| CITRUS | Automated discovery |

## Interpretation

### DA Results
- logFC > 0: Higher in treatment
- p_adj < 0.05: Significant

### DS Results
- Per marker per cluster
- Consider biological relevance

## Experimental Design

- Include biological replicates
- Match batches to conditions
- Consider paired designs

## References

- diffcyt: doi:10.1038/s41467-017-00707-4
- CITRUS: doi:10.1073/pnas.1408792111
