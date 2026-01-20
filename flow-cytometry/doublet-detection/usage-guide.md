# Doublet Detection Usage Guide

## Overview

Doublets are events where two or more cells pass through the laser simultaneously. They can confuse gating and clustering by appearing as intermediate populations.

## Why Remove Doublets

- Doublets appear as false intermediate populations
- Affect clustering accuracy
- Skew population frequencies
- Corrupt rare event detection

## Detection Methods

### FSC-A vs FSC-H
- Most common method
- Works on pulse geometry
- FSC-A = Area under pulse
- FSC-H = Height of pulse
- Singlets: linear relationship
- Doublets: higher A for given H

### FSC-A vs FSC-W
- Some instruments provide Width
- A = H × W
- Doublets have increased width

### Ratio Method
- Calculate FSC-A / FSC-H ratio
- Singlets have consistent ratio
- Doublets have elevated ratio

### DNA Content (CyTOF)
- Use DNA intercalator channels
- Doublets have ~2× DNA signal

## Expected Doublet Rates

| Sample Type | Expected Rate |
|-------------|---------------|
| PBMCs | 1-5% |
| Cell lines | 2-10% |
| Tissue digest | 5-15% |
| Sorted cells | <1% |

## Method Selection

| Situation | Recommended Method |
|-----------|-------------------|
| Standard flow | FSC-A vs FSC-H |
| High doublet rate | Combined FSC + SSC |
| CyTOF | DNA or Event_length |
| No FSC-H | Ratio or DNA |

## Quality Checks

### Good Singlet Gate
- Clear diagonal pattern
- ~95% of events retained
- Clean separation from doublets

### Warning Signs
- Bimodal distribution on diagonal
- Very high doublet rate (>15%)
- Gate cuts into true singlets

## Common Issues

### No FSC-H channel
- Use FSC-W if available
- Calculate derived ratio
- Consider DNA staining

### High doublet rate
- Reduce cell concentration
- Check sample preparation
- Ensure proper disaggregation

### Poor discrimination
- Check pulse settings
- Verify threshold settings
- Consider different gating strategy

## References

- flowAI: doi:10.1093/bioinformatics/btw191
- flowDensity: doi:10.1093/bioinformatics/btu677
