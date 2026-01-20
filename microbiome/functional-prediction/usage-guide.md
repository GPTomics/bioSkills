# Functional Prediction - Usage Guide

## Overview

PICRUSt2 predicts metagenome functional content from 16S/18S marker gene data by inferring gene content from phylogenetically related reference genomes.

## Prerequisites

```bash
# Install PICRUSt2
conda install -c bioconda picrust2

# Or via pip
pip install picrust2
```

## Quick Start

Tell your AI agent what you want to do:
- "Predict functional profiles from my 16S ASV data"
- "Run PICRUSt2 to get KEGG pathway abundances"

## Example Prompts

### Basic Prediction
> "Run PICRUSt2 on my ASV table and sequences"

> "Predict metagenome functional content from my 16S data"

### Quality Assessment
> "Check NSTI values to assess prediction confidence"

> "Filter ASVs with NSTI greater than 2 before functional prediction"

### Downstream Analysis
> "Run differential abundance on PICRUSt2 pathway output"

> "Compare predicted KEGG modules between treatment groups"

### Output Handling
> "Extract stratified output to see which taxa contribute to each pathway"

> "Convert PICRUSt2 output to a format for ALDEx2 analysis"

## What the Agent Will Do

1. Prepare ASV sequences and abundance table
2. Place sequences in reference tree
3. Predict gene content via hidden state prediction
4. Infer metagenome functional profiles
5. Reconstruct metabolic pathways
6. Assess prediction quality via NSTI
7. Generate KEGG ortholog and pathway tables

## Tips

- NSTI <0.5 indicates high confidence predictions
- Filter ASVs with NSTI >2 before inference
- Validate with shotgun metagenomics if possible
- Novel taxa have unreliable predictions
- Report NSTI values in methods section

## Output Files

| File | Description |
|------|-------------|
| `KO_metagenome_out/pred_metagenome_unstrat.tsv` | KEGG ortholog abundances |
| `EC_metagenome_out/pred_metagenome_unstrat.tsv` | EC number abundances |
| `pathways_out/path_abun_unstrat.tsv` | MetaCyc pathway abundances |
| `*_strat.tsv` | Stratified by contributing ASV |

## NSTI Quality Metrics

| NSTI Value | Confidence |
|------------|------------|
| <0.5 | High |
| 0.5-1.5 | Moderate |
| >2 | Low (consider filtering) |

## Limitations

- Assumes functional conservation with phylogeny
- Cannot detect horizontal gene transfer
- Novel taxa have unreliable predictions
- Lower resolution than shotgun metagenomics

## PICRUSt2 vs Shotgun

| Aspect | PICRUSt2 | Shotgun |
|--------|----------|---------|
| Cost | Low | High |
| Input | 16S amplicons | Whole DNA |
| Accuracy | Estimated | Direct |
| Novel functions | Cannot detect | Detectable |
