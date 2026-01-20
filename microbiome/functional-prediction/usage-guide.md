# Functional Prediction with PICRUSt2

## Overview

PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States) predicts metagenome functional content from 16S/18S marker gene data.

## How It Works

1. **Place sequences**: ASVs placed in reference tree
2. **Hidden state prediction**: Infer gene content from relatives
3. **Metagenome inference**: Weight predictions by abundance
4. **Pathway reconstruction**: Aggregate to pathways

## Input Requirements

- ASV sequences (FASTA)
- ASV abundance table (TSV)
- Sequences should be trimmed consistently

## Output Files

| File | Description |
|------|-------------|
| `KO_metagenome_out/pred_metagenome_unstrat.tsv` | KEGG ortholog abundances |
| `EC_metagenome_out/pred_metagenome_unstrat.tsv` | EC number abundances |
| `pathways_out/path_abun_unstrat.tsv` | MetaCyc pathway abundances |
| `*_strat.tsv` | Stratified by contributing ASV |

## Quality Metrics

### NSTI (Nearest Sequenced Taxon Index)
- Distance to nearest reference genome
- **<0.5**: High confidence
- **0.5-1.5**: Moderate confidence
- **>2**: Low confidence (consider filtering)

### Weighted NSTI
- NSTI weighted by ASV abundance
- Reported per sample
- Compare across groups

## Best Practices

1. **Filter high-NSTI ASVs**: Remove >2 NSTI before inference
2. **Use updated reference**: PICRUSt2 uses IMG database
3. **Validate with shotgun**: If possible, compare subset
4. **Report NSTI**: Include in methods

## Limitations

- Assumes functional conservation with phylogeny
- Cannot detect horizontal gene transfer
- Novel taxa have unreliable predictions
- Lower resolution than shotgun metagenomics

## Downstream Analysis

- **Differential abundance**: ALDEx2, ANCOM-BC on pathways
- **PCA/PCoA**: Ordination on functional profiles
- **Enrichment**: Against KEGG modules
- **Integration**: Combine with taxa-level analysis

## Comparison with Shotgun

| Aspect | PICRUSt2 | Shotgun |
|--------|----------|---------|
| Cost | Low | High |
| Input | 16S amplicons | Whole DNA |
| Accuracy | Estimated | Direct |
| Novel functions | Cannot detect | Detectable |
| Best for | Quick assessment | Comprehensive profiling |
