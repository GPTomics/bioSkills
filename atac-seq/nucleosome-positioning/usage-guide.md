# Nucleosome Positioning - Usage Guide

## Overview
Extract nucleosome positions from ATAC-seq fragment size patterns using ATACseqQC or NucleoATAC to understand chromatin structure at promoters and regulatory regions.

## Prerequisites
```r
BiocManager::install(c('ATACseqQC', 'TxDb.Hsapiens.UCSC.hg38.knownGene'))
```

```bash
pip install nucleoatac
```

## Quick Start
Tell your AI agent what you want to do:
- "Analyze nucleosome positioning from my ATAC-seq data"
- "Generate fragment size distribution to check nucleosome periodicity"

## Example Prompts
### Fragment Size Analysis
> "Plot the fragment size distribution from my ATAC-seq BAM to visualize nucleosome patterns"

### Nucleosome Calling
> "Call nucleosome positions from my ATAC-seq data using NucleoATAC"

### TSS Nucleosome Patterns
> "Analyze nucleosome positioning around TSS to visualize NFR and +1/-1 nucleosomes"

### NFR Analysis
> "Extract nucleosome-free region reads and analyze coverage at promoters"

## What the Agent Will Do
1. Calculate fragment size distribution to verify nucleosome periodicity
2. Separate reads by fragment size (NFR, mono-, di-, tri-nucleosomal)
3. Call nucleosome positions using NucleoATAC or ATACseqQC
4. Generate nucleosome occupancy signal tracks
5. Visualize nucleosome patterns around TSS

## Fragment Size Interpretation

| Fragment Size | Origin | Meaning |
|---------------|--------|---------|
| <100 bp | Nucleosome-free | Active promoter/enhancer |
| 180-247 bp | Mono-nucleosome | Well-positioned nucleosome |
| 315-473 bp | Di-nucleosome | Regular spacing |
| 558-615 bp | Tri-nucleosome | Chromatin structure |

## Nucleosome Pattern Interpretation

| Pattern | Meaning |
|---------|---------|
| Strong NFR peak | Active promoter/enhancer |
| Regular spacing | Well-positioned nucleosomes |
| Fuzzy positioning | Dynamic chromatin |

## Tips
- Paired-end data is required for accurate fragment sizing
- Use properly deduplicated BAM files
- Filter for high-quality fragments (MAPQ > 30)
- High mitochondrial contamination can skew fragment distributions
- NucleoATAC provides base-pair resolution nucleosome positions
- ATACseqQC is useful for V-plots around TSS showing nucleosome phasing
