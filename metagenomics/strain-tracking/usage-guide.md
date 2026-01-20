# Strain Tracking - Usage Guide

## Overview
Track bacterial strains at sub-species resolution using genome distance metrics (MASH, fastANI) and strain-level metagenomics (inStrain, sourmash).

## Prerequisites
```bash
# MASH - fast distance estimation
conda install -c bioconda mash

# sourmash - metagenome comparisons
pip install sourmash

# fastANI - average nucleotide identity
conda install -c bioconda fastani

# inStrain - strain-level metagenomics
pip install instrain
```

## Quick Start
Tell your AI agent what you want to do:
- "Calculate pairwise distances between my genome assemblies"
- "Track strains across timepoints in my longitudinal study"
- "Identify which reference genomes are in my metagenome"

## Example Prompts
### Genome Distance Calculations
> "Calculate MASH distances between all genomes in this directory"

> "Run fastANI to determine if these isolates are the same species"

### Outbreak Investigation
> "Cluster my outbreak isolates and identify closely related strains"

> "Find genomes with less than 0.001 MASH distance (same strain)"

### Metagenome Strain Analysis
> "Run sourmash gather to identify reference genomes in my metagenome"

> "Use inStrain to profile strain variation in my sample"

### Longitudinal Tracking
> "Track strain changes across my time-series samples using inStrain"

> "Compare strain populations between treatment and control groups"

## What the Agent Will Do
1. Create genome sketches or signatures for efficient comparison
2. Calculate pairwise distances or ANI values
3. Cluster strains based on distance thresholds
4. Profile within-sample variation for metagenomes
5. Compare strain profiles across samples or timepoints

## Tips
- MASH distance < 0.05 indicates same species (ANI > 95%)
- MASH distance < 0.001 suggests same strain
- sourmash uses MinHash sketches; compatible with large-scale comparisons
- inStrain requires BAM alignment to reference; provides SNV-level resolution
- fastANI is gold standard for species delineation

## Distance Interpretation

| MASH Distance | ANI | Interpretation |
|---------------|-----|----------------|
| 0.00 | 100% | Same strain |
| < 0.05 | > 95% | Same species |
| 0.05-0.10 | 90-95% | Related species |
| > 0.10 | < 90% | Different species |

## inStrain Key Metrics
- **popANI**: Population ANI across reads
- **conANI**: Consensus ANI
- **SNV density**: Variation within sample

## Resources
- [MASH documentation](https://mash.readthedocs.io/)
- [sourmash documentation](https://sourmash.readthedocs.io/)
- [inStrain tutorial](https://instrain.readthedocs.io/)
- [fastANI paper](https://doi.org/10.1038/s41467-018-07641-9)
