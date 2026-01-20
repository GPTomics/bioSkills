# Cytometry QC Usage Guide

## Overview

Comprehensive quality control ensures reliable flow cytometry and CyTOF data by detecting acquisition problems, removing problematic events, and flagging outlier samples.

## When to Use This Skill

- After data acquisition, before analysis
- Batch processing multiple samples
- Troubleshooting inconsistent results
- Multi-site study harmonization
- Longitudinal data comparison

## QC Checks

### Flow Rate Stability
- Events per time unit should be consistent
- Clogs cause rate drops
- Air bubbles cause rate spikes
- CV < 20% is typically acceptable

### Signal Drift
- Detector sensitivity can change over time
- Temperature effects on optics
- Check before/after normalization

### Margin Events
- Events at detector saturation limits
- Cannot be accurately quantified
- Remove before analysis

### Dead Cells
- Take up viability dye
- Should be excluded for most analyses
- Percentage indicates sample quality

## Expected QC Metrics

| Metric | Acceptable | Warning | Fail |
|--------|------------|---------|------|
| Flow rate CV | < 15% | 15-25% | > 25% |
| Signal drift | < 5% | 5-15% | > 15% |
| Margin events | < 1% | 1-5% | > 5% |
| Dead cells | < 10% | 10-30% | > 30% |

## CyTOF-Specific Checks

### Event Length
- Cell size proxy
- Typical: 15-45 for single cells
- Low: debris/dying cells
- High: doublets/aggregates

### DNA Content
- Confirms nucleated cells
- Dead/debris have low DNA
- Use Ir191/Ir193 channels

### Gaussian Parameters
- Push quality metrics
- Center, Width, Residual
- Identify poor acquisitions

## Batch QC

### Outlier Detection
- Event count significantly different
- Flow rate unstable
- Median signals shifted
- Consider excluding or re-acquiring

### Cross-Sample Normalization
- After individual QC passes
- Use bead-normalization skill
- Document normalization parameters

## Common Issues

### High dead cell percentage
- Sample handling issue
- Staining too long
- Temperature problems
- Transport/storage issues

### Flow rate instability
- Clogged sample line
- Air in system
- Sample viscosity
- Low sample volume

### Signal drift
- Laser warmup incomplete
- Temperature fluctuation
- Detector fatigue
- Long acquisition time

## Best Practices

1. Run QC controls at start of each day
2. Include technical replicates
3. Document instrument settings
4. Review QC plots before analysis
5. Set pass/fail thresholds before analysis
6. Archive QC reports with data

## References

- flowAI: doi:10.1093/bioinformatics/btw191
- PeacoQC: doi:10.1002/cyto.a.24501
- CyTOF QC: doi:10.1002/cyto.a.22624
