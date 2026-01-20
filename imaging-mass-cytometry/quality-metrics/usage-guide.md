# Quality Metrics Usage Guide

## Overview

Quality metrics assess IMC data quality at multiple levels: acquisition, preprocessing, and segmentation. Early QC catches issues before downstream analysis.

## When to Use This Skill

- After data acquisition
- Before analysis begins
- Batch processing QC
- Troubleshooting unexpected results
- Multi-site study harmonization

## Key Metrics

### Signal-to-Noise Ratio
| SNR | Interpretation |
|-----|----------------|
| > 5 | Excellent |
| 3-5 | Good |
| 1.5-3 | Acceptable |
| < 1.5 | Poor - review |

### Tissue Coverage
| Coverage | Status |
|----------|--------|
| > 50% | Good |
| 30-50% | Acceptable |
| < 30% | Poor - may affect analysis |

### Channel Correlation
- Expected: Related markers (CD3/CD45)
- Unexpected high: Possible spillover
- Check against panel design

## Common Issues

### Low SNR
- Antibody concentration issue
- Poor tissue fixation
- Acquisition parameters
- Consider normalization

### Unexpected correlations
- Spillover from adjacent channels
- Apply spillover correction
- Check panel mass spacing

### Acquisition artifacts
- Hot pixels: Dead detector elements
- Striping: Ion beam instability
- Saturation: Overloaded detector

### Tissue quality
- Low coverage: Tissue loss
- High fragmentation: Poor sectioning
- Dead regions: Ablation issues

## QC Thresholds

### Pass Criteria (typical)
- Mean SNR > 2
- All channels SNR > 1.5
- Tissue coverage > 30%
- No saturation > 1%
- No striping artifacts

### Flag for Review
- Any channel SNR < 1.5
- Unexpected correlations > 0.7
- Tissue fragmentation > 0.5
- Hot pixels > 0.1%

### Fail
- Mean SNR < 1
- Coverage < 10%
- Major striping/banding
- Extensive saturation

## Batch QC

### Outlier Detection
- Compare across samples
- Flag > 2 SD from median
- Review before exclusion

### Documentation
- Record QC decisions
- Note excluded samples
- Track batch effects

## Best Practices

1. Run QC before analysis
2. Define thresholds before looking at data
3. Document all QC decisions
4. Keep QC reports with data
5. Review failed samples individually
6. Consider reacquisition for failures

## References

- IMC quality control: doi:10.1038/s41596-021-00508-0
- steinbock QC: doi:10.1038/s41592-021-01308-y
