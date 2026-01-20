# TF Footprinting - Usage Guide

## Overview
Identify transcription factor binding sites from ATAC-seq data by detecting characteristic "footprints" where DNA-bound proteins protect against Tn5 cutting.

## Prerequisites
```bash
pip install tobias pyBigWig numpy matplotlib
conda install -c bioconda samtools rgt
```

## Quick Start
Tell your AI agent what you want to do:
- "Run TF footprinting analysis on my ATAC-seq data"
- "Compare TF activity between treatment and control conditions"

## Example Prompts
### Basic Footprinting
> "Run TOBIAS footprinting analysis on my ATAC-seq BAM to identify bound transcription factors"

### Differential Footprinting
> "Compare TF binding between treatment and control ATAC-seq samples using TOBIAS"

### Visualize Footprints
> "Plot aggregate footprints for CTCF binding sites from my ATAC-seq data"

### Motif-Specific Analysis
> "Identify footprints for specific TF motifs from JASPAR in my accessible regions"

## What the Agent Will Do
1. Filter BAM for nucleosome-free reads (<100bp fragments)
2. Correct for Tn5 sequence bias using TOBIAS ATACorrect
3. Calculate footprint scores across accessible regions
4. Detect bound TF sites using motif analysis (BINDetect)
5. Generate per-TF binding reports and visualization

## Requirements for Footprinting
- **High depth**: >50M uniquely mapped reads
- **NFR reads**: Filter for fragments <100bp
- **Peak regions**: Accessible chromatin regions
- **TF motifs**: JASPAR, HOCOMOCO, or custom

## Interpreting Results

### bindetect_results.txt Columns
| Column | Description |
|--------|-------------|
| name | TF name |
| motif_id | JASPAR ID |
| n_detected | Number of binding sites |
| mean_score | Average footprint score |
| differential_score | Difference between conditions |
| pvalue | Statistical significance |

### Quality Indicators
- Good: Clear V-shaped dip at motif center, symmetric shoulders, many bound sites
- Poor: Flat or noisy signal, asymmetric pattern, few detected sites

## Tips
- Footprinting requires high sequencing depth (>50M reads)
- Use only nucleosome-free reads (<100bp) for best results
- Tn5 bias correction is critical for accurate footprinting
- Consider merging replicates if individual samples have low depth
- HINT-ATAC (RGT) is an alternative tool if TOBIAS is not available
