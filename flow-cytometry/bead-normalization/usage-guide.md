# Bead Normalization Usage Guide

## Overview

Bead-based normalization corrects for signal drift and batch effects in mass cytometry (CyTOF) and high-parameter flow cytometry by using reference beads with known properties.

## When to Use This Skill

- CyTOF experiments (EQ beads standard)
- Multi-batch experiments
- Longitudinal studies
- Multi-site collaborations
- Observed signal drift

## Types of Normalization

### EQ Bead Normalization (CyTOF)
- Standard for CyTOF
- Uses metal-containing beads
- Corrects instrument drift
- Typically done with Fluidigm normalizer

### Reference Sample Normalization
- Use same biological sample across batches
- Build transformation model
- Apply to all samples
- CytoNorm, Harmony methods

### Quantile Normalization
- Force same distribution across samples
- Simple but aggressive
- Use cautiously with biological heterogeneity

## EQ Beads

Fluidigm EQ Four Element Calibration Beads contain:
- Cerium-140
- Europium-151
- Europium-153
- Holmium-165
- Lutetium-175

## Workflow

1. **Identify beads** - Gate on bead-positive events
2. **Calculate reference** - Median intensity per channel
3. **Compute factors** - Ratio to reference
4. **Apply correction** - Multiply or LOESS smooth
5. **Remove beads** - Clean data for analysis

## Drift Patterns

| Type | Cause | Solution |
|------|-------|----------|
| Linear drift | Ion source degradation | Linear correction |
| Step change | Tuning adjustment | Segment normalization |
| Random fluctuation | Unstable conditions | LOESS smoothing |

## Quality Metrics

### Good Normalization
- CV of bead channels < 10%
- No trend in residuals
- Biological patterns preserved

### Warning Signs
- Persistent drift after correction
- Over-correction artifacts
- Loss of biological signal

## Common Issues

### Bead identification failure
- Adjust gating threshold
- Check bead concentration
- Verify channel names

### Over-normalization
- Use appropriate span
- Preserve biological variation
- Validate with known markers

### Batch effects remain
- Consider CytoNorm
- Add reference samples
- Check other confounders

## Best Practices

1. Include EQ beads in every run
2. Run reference sample per batch
3. Normalize before analysis
4. Validate with known biology
5. Document normalization parameters

## References

- CyTOF normalization: doi:10.1002/cyto.a.22271
- CytoNorm: doi:10.1002/cyto.a.24158
- Bead-based QC: doi:10.1002/cyto.a.22624
