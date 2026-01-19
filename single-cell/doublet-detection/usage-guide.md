# Doublet Detection Usage Guide

This guide covers detecting and removing doublets from single-cell RNA-seq data.

## What are Doublets?

Doublets occur when two or more cells are captured in the same droplet. They:
- Appear as intermediate cell populations
- Inflate cell counts
- Can lead to false biological conclusions

## When to Run Doublet Detection

Run doublet detection:
- After initial QC (gene/cell filtering)
- Before normalization and clustering
- On each sample separately (if pooled)

## Method Selection

### Scrublet (Python)
- Fast and simple
- Good for quick analysis
- Works well with Scanpy

### DoubletFinder (R)
- Most widely used in R
- Requires parameter optimization
- Works with Seurat

### scDblFinder (R)
- Fastest R method
- Machine learning based
- Often most accurate

## Expected Doublet Rates

Use the 10X formula: ~0.8% per 1,000 cells loaded

| Cells | Rate |
|-------|------|
| 5,000 | 4% |
| 10,000 | 8% |
| 15,000 | 12% |

## Requirements

```bash
# Python
pip install scrublet scanpy

# R
install.packages('Seurat')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
BiocManager::install('scDblFinder')
```

## Troubleshooting

### No doublets detected
- Check expected_doublet_rate
- Try lower threshold
- Verify data quality

### Too many doublets
- Lower expected rate
- Raise score threshold
- Check for batch effects

### Bimodal distribution unclear
- Use scDblFinder instead
- Set manual threshold
- Check data quality
