# Hit Calling Usage Guide

## Overview

Multiple methods exist for calling hits in CRISPR screens. The choice depends on screen design, reference data availability, and desired stringency.

## Method Comparison

| Method | Approach | Strengths | Best For |
|--------|----------|-----------|----------|
| MAGeCK RRA | Rank aggregation | No training data needed | General screens |
| MAGeCK MLE | Maximum likelihood | Multi-condition | Complex designs |
| BAGEL2 | Bayesian | Uses reference genes | Well-characterized systems |
| drugZ | Z-score | Simple, interpretable | Drug screens |

## Choosing Thresholds

### MAGeCK
- FDR < 0.1 for discovery
- FDR < 0.05 for high confidence

### BAGEL2
- BF > 3 suggestive
- BF > 5 strong evidence
- BF > 10 very strong

### Z-scores
- |Z| > 2 (~p < 0.05)
- |Z| > 3 (~p < 0.003)

## Consensus Approach

Combine multiple methods:
1. Run MAGeCK, BAGEL2, drugZ
2. Take genes called by 2+ methods
3. Reduces method-specific biases

## Negative vs Positive Selection

### Negative (Dropout)
- Essential genes
- Drug sensitivity
- Look at neg|score, synth

### Positive (Enrichment)
- Resistance genes
- Growth advantage
- Look at pos|score, supp

## References

- BAGEL2: doi:10.1186/s13059-019-1749-z
- drugZ: doi:10.1186/s13073-019-0665-3
