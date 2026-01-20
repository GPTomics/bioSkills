# CRISPR Screen QC Usage Guide

## Overview

Quality control is critical for CRISPR screens. Poor library representation or technical issues can lead to false positives/negatives.

## Key QC Metrics

### Library Representation
- **Zero-count sgRNAs**: <1% ideal, <5% acceptable
- **Low-count sgRNAs**: <10% with <30 reads

### Read Distribution
- **Gini index**: <0.2 ideal, <0.3 acceptable
- Even distribution preferred

### Replicate Correlation
- **Pearson r**: >0.8 ideal, >0.6 minimum
- Log-transformed counts

### Essential Gene Recovery
- **AUC**: >0.85 ideal for gold-standard essentials
- Validates screen worked

## Common Issues

### Low Library Representation
- Causes: Low MOI, bottleneck, poor amplification
- Fix: Increase cell number, optimize PCR

### Skewed Distribution
- Causes: PCR bias, library quality
- Fix: Reduce PCR cycles, new library prep

### Poor Replicate Correlation
- Causes: Technical variation, batch effects
- Fix: Process replicates together

### No Essential Gene Dropout
- Causes: Insufficient selection time, wrong cell line
- Fix: Extend selection, verify model system

## Reference Gene Sets

- **Essential**: DepMap common essentials, Hart et al.
- **Non-essential**: Olfactory receptors, non-expressed genes

## References

- Hart et al. essentials: doi:10.1016/j.molcel.2017.06.014
- DepMap: https://depmap.org/
