# Cytometry QC - Usage Guide

## Overview
Comprehensive quality control ensures reliable flow cytometry and CyTOF data by detecting acquisition problems, removing problematic events, and flagging outlier samples.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowCore', 'flowAI', 'CATALYST'))

# Optional for advanced QC
BiocManager::install('PeacoQC')
```

## Quick Start
Tell your AI agent what you want to do:
- "Run quality control on my flow cytometry data"
- "Check for flow rate anomalies and signal drift"
- "Identify and flag low-quality samples"

## Example Prompts
### Event-Level QC
> "Remove margin events and anomalous time segments"
> "Run flowAI to automatically clean my FCS files"
> "Apply PeacoQC to detect acquisition problems"

### Sample-Level QC
> "Check flow rate stability for each sample"
> "Flag samples with high dead cell percentages"
> "Identify outlier samples based on event counts"

### CyTOF-Specific QC
> "Filter events by Event_length and DNA content"
> "Check Gaussian parameters for acquisition quality"
> "Generate a QC report for my CyTOF batch"

### Batch QC
> "Compare QC metrics across all samples in my experiment"
> "Create a QC summary table for my dataset"
> "Flag samples that fail QC thresholds"

## What the Agent Will Do
1. Load QC packages (flowAI, PeacoQC, or CATALYST)
2. Check flow rate stability over time
3. Identify margin events and signal anomalies
4. Remove or flag problematic events
5. Generate QC report with pass/fail metrics

## Tips
- Run QC before any downstream analysis
- flowAI handles flow rate, signal drift, and margin events
- For CyTOF, filter by Event_length (15-45) and DNA content
- Dead cells (>10%) indicate sample handling issues
- Document QC thresholds before analyzing data

## QC Thresholds

| Metric | Acceptable | Warning | Fail |
|--------|------------|---------|------|
| Flow rate CV | <15% | 15-25% | >25% |
| Signal drift | <5% | 5-15% | >15% |
| Margin events | <1% | 1-5% | >5% |
| Dead cells | <10% | 10-30% | >30% |

## CyTOF-Specific Checks

| Check | Purpose | Typical Range |
|-------|---------|---------------|
| Event_length | Cell size proxy | 15-45 |
| DNA (Ir191/193) | Confirms nucleated cells | Bimodal peak |
| Gaussian Center/Width | Push quality | Per run baseline |

## Common Issues

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| High dead cells | Sample handling, temperature | Improve protocol |
| Flow rate instability | Clog, air bubble | Clean instrument |
| Signal drift | Laser warmup, temperature | Allow warmup, normalize |

## References
- flowAI: doi:10.1093/bioinformatics/btw191
- PeacoQC: doi:10.1002/cyto.a.24501
- CyTOF QC: doi:10.1002/cyto.a.22624
