# Batch Correction - Usage Guide

## Overview
Batch effects in CRISPR screens arise from technical variation between screening rounds, different plasmid preps, or varying experimental conditions. Correction is essential for combining data across batches.

## Prerequisites
```bash
pip install pandas numpy scipy
pip install combat  # or pycombat
# R packages
# BiocManager::install(c('sva', 'limma'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Correct batch effects between my screening rounds"
- "Normalize counts across different batches before analysis"
- "Check if my batches need correction using PCA"

## Example Prompts

### Batch Assessment
> "Run PCA on my CRISPR screen data colored by batch to see if there are batch effects."

> "Calculate the correlation between replicates within and across batches. Is batch correction needed?"

> "Visualize the distribution of control guides across batches to assess technical variation."

### Correction Methods
> "Apply ComBat batch correction to my count matrix before running MAGeCK."

> "Use median normalization to correct for sequencing depth differences between batches."

> "Normalize my screen data using non-targeting controls as reference."

### Multi-Batch Analysis
> "I have screens from 3 different weeks. Correct batch effects and combine for joint analysis."

> "Integrate data from two different screening facilities accounting for site-specific effects."

### Validation
> "After batch correction, check that essential gene dropout is preserved."

> "Verify batch correction worked by re-running PCA. Batches should now overlap."

## What the Agent Will Do
1. Assess batch effects using PCA and correlation
2. Select appropriate correction method
3. Apply normalization or batch correction
4. Verify correction preserved biological signal
5. Generate before/after comparison plots
6. Output corrected count matrix

## Tips
- Always visualize batch effects before correction (PCA, correlation)
- Include same controls in all batches for validation
- Median normalization is simplest, ComBat is most powerful
- Check that essential gene dropout is preserved after correction
- Over-correction can remove biological signal - validate carefully

## Method Selection

| Situation | Recommended Method |
|-----------|-------------------|
| Minor depth variation | Median normalization |
| Moderate batch effects | Size factor normalization |
| Strong batch effects | ComBat |
| Biological signal priority | Control-based normalization |

## QC Metrics

| Metric | Excellent | Good | Investigate |
|--------|-----------|------|-------------|
| Within-batch correlation | > 0.95 | 0.9-0.95 | < 0.9 |
| Cross-batch correlation | > 0.90 | 0.8-0.90 | < 0.8 |
| Batch effect (PCA) | Not visible | Minor | Dominant |

## Validation Checklist
- PCA no longer separates by batch
- Replicate correlations improved
- Control guide distributions aligned
- Essential gene depletion preserved
- Non-targeting controls stable

## References
- MAGeCK batch correction: doi:10.1186/s13059-015-0843-6
- ComBat: doi:10.1093/biostatistics/kxj037
- Screen normalization: doi:10.1038/nmeth.3935
