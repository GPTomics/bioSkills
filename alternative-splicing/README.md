# alternative-splicing

## Overview
Alternative splicing analysis for RNA-seq data covering event quantification, differential splicing detection, isoform switching analysis, and visualization.

**Tool type:** mixed | **Primary tools:** rMATS-turbo, SUPPA2, leafcutter, IsoformSwitchAnalyzeR

## Skills
| Skill | Description |
|-------|-------------|
| splicing-quantification | Quantify PSI values for splicing events from RNA-seq |
| differential-splicing | Detect differential splicing between conditions |
| isoform-switching | Analyze isoform switches and functional consequences |
| sashimi-plots | Visualize splicing events with junction coverage |
| splicing-qc | Quality control metrics for splicing analysis |
| single-cell-splicing | Single-cell resolution splicing analysis |

## Example Prompts
- "Quantify exon skipping events from my RNA-seq BAMs"
- "Find differential splicing between tumor and normal"
- "Create sashimi plots for significant splicing events"
- "Check if my data has sufficient junction coverage for splicing analysis"
- "Analyze splicing heterogeneity in my scRNA-seq data"

## Requirements
```bash
# Python
pip install suppa rmats-turbo ggsashimi rseqc brie2

# R
BiocManager::install(c('IsoformSwitchAnalyzeR', 'leafcutter'))
```

## Related Skills

- **read-alignment** - STAR 2-pass mode for junction detection
- **differential-expression** - Compare gene-level vs splicing changes
- **rna-quantification** - Transcript TPM for SUPPA2
