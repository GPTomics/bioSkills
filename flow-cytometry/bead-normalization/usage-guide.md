# Bead Normalization - Usage Guide

## Overview
Bead-based normalization corrects for signal drift and batch effects in CyTOF and high-parameter flow cytometry by using reference beads with known properties.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('CATALYST', 'flowCore'))

# Or use Fluidigm normalizer for EQ beads
```

## Quick Start
Tell your AI agent what you want to do:
- "Normalize my CyTOF data using EQ beads"
- "Correct for signal drift across my acquisition"
- "Apply batch correction using reference samples"

## Example Prompts
### EQ Bead Normalization
> "Identify EQ beads and normalize my CyTOF run"
> "Correct for instrument drift using bead channels"
> "Remove bead events after normalization"

### Batch Correction
> "Apply CytoNorm to harmonize data across batches"
> "Use reference samples to build a batch correction model"
> "Normalize my multi-site CyTOF study"

### Quality Assessment
> "Plot bead signal over time to visualize drift"
> "Show before/after normalization comparison"
> "Calculate CV of bead channels across the run"

## What the Agent Will Do
1. Identify bead events using bead-specific channels
2. Calculate reference intensity (median per channel)
3. Compute normalization factors over time
4. Apply correction (linear or LOESS smoothing)
5. Remove bead events from final dataset

## Tips
- EQ beads contain Ce-140, Eu-151, Eu-153, Ho-165, Lu-175
- Always include beads in every CyTOF run
- CV of bead channels should be <10% after normalization
- For multi-batch studies, include a reference sample per batch
- CytoNorm provides batch-to-batch normalization using reference samples

## Normalization Types

| Type | Use Case | Method |
|------|----------|--------|
| EQ Bead | Within-run drift | Fluidigm normalizer or CATALYST |
| Reference sample | Between-batch | CytoNorm, Harmony |
| Quantile | Simple normalization | Force same distribution (use cautiously) |

## Drift Patterns

| Pattern | Cause | Solution |
|---------|-------|----------|
| Linear drift | Ion source degradation | Linear correction |
| Step change | Tuning adjustment | Segment normalization |
| Random fluctuation | Unstable conditions | LOESS smoothing |

## Quality Metrics

| Metric | Good | Warning |
|--------|------|---------|
| Bead CV | <10% | >15% |
| Residual trend | None | Persistent drift |
| Biological signal | Preserved | Lost/distorted |

## References
- CyTOF normalization: doi:10.1002/cyto.a.22271
- CytoNorm: doi:10.1002/cyto.a.24158
- Bead-based QC: doi:10.1002/cyto.a.22624
