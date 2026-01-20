# Library Design Usage Guide

## Overview

CRISPR library design involves selecting optimal sgRNAs for genetic screens, balancing efficiency, specificity, and practical cloning considerations.

## When to Use This Skill

- Designing new knockout/activation libraries
- Optimizing existing library designs
- Creating focused sublibaries
- Validating library composition

## Design Considerations

### sgRNA Selection Criteria
| Criterion | Optimal | Acceptable |
|-----------|---------|------------|
| GC content | 40-60% | 30-70% |
| Length | 20 nt | 18-22 nt |
| Poly-T runs | None | < 4 Ts |
| 5' nucleotide | G | A, C |

### Off-Target Risk
- Use validated algorithms (CRISPOR, Azimuth)
- Filter guides with >3 off-targets in exons
- Consider off-target in essential genes

### Library Size
| Screen Type | Guides/Gene | Minimum Cells |
|-------------|-------------|---------------|
| Dropout | 4-6 | 500x coverage |
| Activation | 3-4 | 300x coverage |
| Focused | 6-10 | 1000x coverage |

## Control Design

### Essential Controls
- Include 20-50 essential gene guides
- Should drop out in all conditions
- Validates library/screen quality

### Non-Targeting Controls
- 100-500 non-targeting guides
- Use for background estimation
- Spread across library synthesis

### Safe Harbor Controls
- AAVS1, ROSA26 targeting
- Validate Cas9 activity
- Should show no phenotype

## Cloning Strategy

### Arrayed Oligo Synthesis
- For libraries < 1,000 guides
- Higher per-guide quality
- Individual validation possible

### Pooled Array Synthesis
- For libraries > 1,000 guides
- More cost-effective
- Requires representation QC

### Subpool Strategy
- Divide large libraries into subpools
- Enables focused validation
- Simplifies coverage calculations

## Common Issues

### Low guide scores
- Gene lacks good PAM sites
- Try extended search window
- Consider alternative PAMs (NAG)

### Too many off-targets
- Use stricter filtering
- Prioritize exon-targeting
- Check against recent genome build

### Uneven representation
- Array synthesis variance
- Amplification bias
- May need rebalancing

## Quality Metrics

### Good Library
- 4+ guides per gene
- >90% guides with GC 30-70%
- <5% guides with off-targets
- Controls: 5-10% of library

## Best Practices

1. Use latest genome annotation
2. Validate computationally predicted scores
3. Include comprehensive controls
4. Document design parameters
5. Verify representation after cloning
6. Bank library aliquots

## References

- CRISPOR: doi:10.1093/nar/gkw441
- Azimuth 2.0: doi:10.1038/nbt.3437
- Library design: doi:10.1038/nrg.2017.97
