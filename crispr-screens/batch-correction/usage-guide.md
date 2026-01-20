# Batch Correction Usage Guide

## Overview

Batch effects in CRISPR screens arise from technical variation between screening rounds, different plasmid preps, or varying experimental conditions. Correction is essential for combining data across batches.

## When to Use This Skill

- Combining multiple screening rounds
- Different cell line passages
- Different plasmid preparations
- Multi-site screens
- Screens performed over time

## Types of Batch Effects

### Systematic Effects
- Library representation differences
- Sequencing depth variation
- Cell culture conditions

### Replicate-Level Effects
- Infection efficiency
- Selection timing
- DNA extraction quality

## Normalization Methods

### Median Normalization
- Simple and robust
- Best for minor batch effects
- Preserves relative relationships

### Size Factor Normalization
- DESeq2-style approach
- Accounts for library composition
- Good for count data

### Quantile Normalization
- Forces same distribution
- Aggressive correction
- May over-correct biological signal

### Control-Based Normalization
- Uses non-targeting controls
- Preserves biological variation
- Requires sufficient controls

## When to Use Each Method

| Situation | Recommended Method |
|-----------|-------------------|
| Minor technical variation | Median normalization |
| Sequencing depth variation | Size factors |
| Strong batch effects | ComBat |
| Biological signal priority | Control-based |
| Multi-batch comparison | Quantile normalization |

## QC Metrics

### Replicate Correlation
| Pearson r | Interpretation |
|-----------|----------------|
| > 0.95 | Excellent |
| 0.9-0.95 | Good |
| 0.8-0.9 | Acceptable |
| < 0.8 | Poor - investigate |

### Batch Separation
- PCA should not separate by batch
- Batch effect ratio < 1 is good
- Consider excluding high-effect batches

## Common Issues

### Over-correction
- Removes true biological signal
- Use lighter normalization
- Validate with controls

### Under-correction
- Batch effects persist
- Try stronger methods
- Consider excluding batches

### Different baseline
- Control guides differ between batches
- Normalize to controls first
- Consider stratified analysis

## Best Practices

1. Include same controls in all batches
2. Process replicates together
3. Check batch effects before correction
4. Validate correction with controls
5. Document normalization method
6. Compare results with/without correction

## Validation

### Check Before Correction
- PCA colored by batch
- Replicate correlations
- Control guide distributions

### Check After Correction
- Batch effect reduced in PCA
- Essential gene depletion preserved
- Control guide distributions aligned

## References

- MAGeCK batch correction: doi:10.1186/s13059-015-0843-6
- ComBat: doi:10.1093/biostatistics/kxj037
- Screen normalization: doi:10.1038/nmeth.3935
