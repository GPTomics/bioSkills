# Taxonomy Assignment

## Overview

Taxonomic assignment classifies ASVs or OTUs to taxonomic ranks (Kingdom through Species) using reference databases.

## Classification Methods

### Naive Bayes Classifier
- **Used by**: DADA2, QIIME2 (sklearn)
- **Pros**: Fast, handles novel sequences
- **Cons**: May overclassify

### Exact Matching
- **Used by**: VSEARCH, BLAST
- **Pros**: High precision
- **Cons**: Misses novel taxa, slower

### Hybrid
- **Used by**: SINTAX, IDTAXA
- **Pros**: Balances speed and accuracy

## Reference Databases

### SILVA (16S/18S)
- Most comprehensive for general use
- Version 138.1 is current
- Includes bacteria, archaea, eukaryotes

### GTDB (16S)
- Genome-based taxonomy
- Better for novel/environmental samples
- More consistent naming

### UNITE (ITS)
- Gold standard for fungi
- Regular updates
- Species hypotheses for unknown taxa

### RDP (16S)
- Historical standard
- Less frequently updated
- 6-rank taxonomy

## Confidence Thresholds

Typical bootstrap thresholds:
- **Phylum**: 50-70%
- **Family**: 70-80%
- **Genus**: 80-90%
- **Species**: 95%+ (rarely achieved for 16S V4)

## Common Issues

### Unassigned taxa
- Check database coverage
- May be genuine novel taxa
- Try different database

### Inconsistent naming
- SILVA vs GTDB naming differs
- Stick to one database per study

### Low species-level assignment
- Normal for short amplicons
- V4 region has limited species resolution
- Consider longer reads or different region

## Output Formats

```
# Typical taxonomy table format
ASV    Kingdom   Phylum          Class          Order          Family         Genus
ASV1   Bacteria  Firmicutes      Clostridia     Clostridiales  Lachnospiraceae  Blautia
ASV2   Bacteria  Bacteroidota    Bacteroidia    Bacteroidales  Bacteroidaceae   Bacteroides
```
