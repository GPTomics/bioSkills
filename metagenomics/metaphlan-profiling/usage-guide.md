# MetaPhlAn Profiling - Usage Guide

## Overview
MetaPhlAn uses clade-specific marker genes to profile shotgun metagenomes to species/SGB-level relative abundance. Its percentages are a cell fraction (genome-size-normalized taxonomic abundance), not a fraction of reads - which is why they must never be merged with Kraken/Bracken read fractions. MetaPhlAn 4 quantifies database-absent taxa as uSGBs, is high-precision and low-recall, and reports a separate unknown fraction whose default behavior changed in version 4.2.

## Prerequisites
```bash
conda install -c bioconda metaphlan

# Database downloads on first run (several GB); pin and place it explicitly:
metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 --db_dir /path/to/db
```

Conceptual prerequisites:
- A MetaPhlAn percentage is a cell fraction; a Kraken/Bracken percentage is a read fraction. Keep them apart.
- The atomic taxon is the SGB (species-level), not a strain; uSGBs have no Latin name but are real quantified taxa.
- Pin `--index` and keep one version across a study; the unknown-fraction default flip (4.2) is a hidden batch effect.
- MetaPhlAn is high precision, low recall - it cannot see clades whose markers are not in the database.

## Quick Start
Tell your AI agent what you want to do:
- "Profile species abundances for my gut metagenome with MetaPhlAn, pinning the database index"
- "Save the mapping file so I can re-profile at different levels without realigning"
- "Include the unknown fraction so my environmental sample abundances are honest"
- "Decide whether MetaPhlAn or mOTUs3 fits my under-characterized environment"

## Example Prompts

### Basic profiling
> "Run MetaPhlAn 4 on sample_R1.fastq.gz and sample_R2.fastq.gz at species level, pin the index, and cache the mapping file."

### Quantifying novel taxa
> "My samples are from rumen and a lot of the community is uncharacterized. Profile with MetaPhlAn 4 keeping uSGBs, and tell me how large the unknown fraction is."

### Marker-gene vs k-mer choice
> "I need maximum recall of novel species in a soil metagenome. Help me decide between MetaPhlAn, mOTUs3, and sourmash gather."

### Merging a study
> "Profile all my samples on one pinned index and merge them into one abundance table for compositional analysis."

## What the Agent Will Do
1. Confirm the MetaPhlAn version and pin the database `--index`.
2. Profile each sample, caching the read-to-marker mapping for cheap re-profiling.
3. Keep the unknown fraction consistent across the study and report its size for environmental samples.
4. Merge profiles built on the same index, and convert SGBs to GTDB names if needed.
5. State that the percentages are cell fractions and must not be merged with Kraken/Bracken read fractions.
6. Hand off compositional differential abundance (CLR/ANCOM-BC) to metagenome-visualization.

## Tips
- Pass paired reads as one comma-separated argument; MetaPhlAn treats them as two single-end files.
- Save `--mapout` once, then re-profile at different levels/settings from it for free.
- A low mapping rate is normal; a large unknown fraction means database-absent community, not a failure.
- SGBs are species-level. Do not call `--tax_lev t` "strain-level" - strains are StrainPhlAn.
- For recall-limited environments, reach for mOTUs3 or sourmash gather rather than lowering MetaPhlAn thresholds.

## MetaPhlAn vs Kraken2 (frame by precision/recall and abundance type, not "accuracy")

| Axis | MetaPhlAn 4 | Kraken2 + Bracken |
|------|-------------|-------------------|
| Method | clade-specific marker genes | k-mer LCA + Bayesian reestimation |
| Abundance type | cell fraction (taxonomic) | read fraction (sequence) |
| Precision / recall | high precision, lower recall | high recall, lower precision |
| Database-absent taxa | uSGBs (if binned) | assigned to nearest present relative |

## Common Issues

### No Database Found
```bash
metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403
```

### Low Mapping Rate
Normal - only marker genes are targeted. A large unknown fraction means uncharacterized community; a very low rate with low microbial yield suggests host contamination (see contamination-controls).

### Output All Zeros
- Check the input file is not empty and `--input_type` matches the format.
- The sample may be host-dominated or have very low microbial content.

## Related Skills

- kraken-classification - K-mer read classification; read fraction, not cell fraction
- abundance-estimation - Compositional handling and cross-tool comparison
- strain-tracking - StrainPhlAn strain resolution below the SGB
- functional-profiling - HUMAnN reuses the MetaPhlAn profile
- metagenome-visualization - Compositional stats and plotting
- workflows/metagenomics-pipeline - End-to-end shotgun profiling

## Resources
- [MetaPhlAn GitHub](https://github.com/biobakery/MetaPhlAn)
- [MetaPhlAn Wiki](https://github.com/biobakery/MetaPhlAn/wiki)
